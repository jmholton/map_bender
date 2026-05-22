#!/usr/bin/env ccp4-python
# -*- coding: utf-8 -*-
"""altalign.py - standalone alternate-indexing + origin search.

Finds the crystallographic (altindex x symop x origin) transformation that
best aligns a moving crystal onto a reference crystal, using an LSQ-based
solver:

  1. Kabsch rigid-body fit of matched CA atoms -> continuous R_lsq, t_lsq.
  2. Enumerate discrete (altindex x symop x origin) candidates via
     bendfinder's basis-change altindex enumeration (_get_altindex_ops)
     + the SG's own symops + an origin grid.
  3. Score every candidate two ways, both in honest cartesian RMSD (no
     per-atom fractional wrapping):
       - discrete:  apply R_cart + the exact crystallographic t.
       - comfit:    apply R_cart + the COM-optimal continuous t (= the
                    best rigid placement for that fixed rotation).
  4. Rank by rotational deviation from R_lsq; among the rotation-matched
     candidates, the comfit RMSD is the residual the bender must absorb,
     and the discrete RMSD tells you whether the exact crystallographic
     op alone is good enough to hand straight to refmac.

The point of keeping this standalone: iterate on the altindex/origin
search without running a full fitreso scan each time.

Usage:
    ccp4-python altalign.py <moving.pdb> <reference.pdb> \
        [out.pdb] [moving.mtz] [out.mtz]

If moving.mtz is given, the chosen discrete altindex op is applied to it
via the structure-factor transformation theorem (F'(h) = F(Rᵀh)·e^{2πi h·t})
so that the moving PDB and MTZ stay mutually consistent and refmac can be
run on the pair as usual.
"""
import os
import sys
import numpy as np
import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bendfinder as bf   # noqa: E402


def _rot_angle_deg(R):
    """Rotation angle (deg) of an orthogonal matrix R."""
    cos_t = max(-1.0, min(1.0, (np.trace(R) - 1.0) / 2.0))
    return float(np.degrees(np.arccos(cos_t)))


def search(mov_pdb, ref_pdb, max_M=1, det_max=1, verbose=True):
    """Run the LSQ-based altindex+origin search.

    Returns a dict with the LSQ baseline and a ranked candidate list.
    Each candidate is a dict: R_frac, t_frac, R_cart, drot, rmsd_discrete,
    rmsd_comfit, t_comfit_cart.
    """
    mov_st = gemmi.read_structure(mov_pdb)
    ref_st = gemmi.read_structure(ref_pdb)
    A, B = bf._collect_matched_ca(mov_st, ref_st)
    if len(A) < 3:
        raise SystemExit(f"altalign: <3 matched CA between {mov_pdb} and {ref_pdb}")

    rmsd_before = float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1))))
    R_lsq, t_lsq, rmsd_lsq = bf._kabsch(A, B)

    cell    = mov_st.cell
    sg      = mov_st.find_spacegroup()
    sg_name = mov_st.spacegroup_hm
    O, Finv = bf._ortho_pair(cell)
    G = O.T @ O

    if verbose:
        print(f"altalign: {mov_pdb}  ->  {ref_pdb}")
        print(f"  matched CA      : {len(A)}")
        print(f"  RMSD before fit : {rmsd_before:.2f} A")
        print(f"  LSQ rigid fit   : {rmsd_lsq:.2f} A   "
              f"(rotation {_rot_angle_deg(R_lsq):.1f} deg)")
        rc = (mov_st.cell.parameters, ref_st.cell.parameters)
        if not np.allclose(rc[0], rc[1], atol=1e-2):
            print(f"  NOTE cells differ: mov {tuple(round(x,2) for x in rc[0])}")
            print(f"                     ref {tuple(round(x,2) for x in rc[1])}")
        print(f"  enumerating altindex x symop x origin "
              f"(max_M={max_M}, det_max={det_max}) ...", flush=True)

    # Deduplicate (R_frac, t_frac mod 1) candidates from the enumerator.
    seen = set()
    cands = []
    B_com = B.mean(axis=0)
    for R_frac, t_frac in bf._enum_alt_rot_origin_candidates(cell, sg, sg_name):
        if not bf._is_metric_preserving(R_frac, G):
            continue
        R_int = np.round(R_frac).astype(int)
        t_mod = np.round(np.asarray(t_frac) % 1.0, 6)
        key = (tuple(R_int.flatten()), tuple(t_mod))
        if key in seen:
            continue
        seen.add(key)

        R_cart = O @ R_frac @ Finv
        A_rot  = (R_cart @ A.T).T

        # Discrete crystallographic translation.  The enumerated t_frac is
        # only defined mod 1; pick the whole-lattice image of the *whole*
        # moving model (one integer triplet for all atoms) that lands it on
        # the reference.  An integer-lattice shift is exact crystallography
        # — and invisible to the MTZ, since exp(2πi·h·integer)=1 — so we
        # fold it straight into the translation that gets applied.
        t_arr  = np.asarray(t_frac, dtype=float)
        fit0   = A_rot + O @ t_arr
        d_frac = Finv @ (fit0 - B).mean(axis=0)
        img    = np.round(d_frac)
        t_apply = t_arr - img                      # full fractional translation
        fit_disc = A_rot + O @ t_apply
        rmsd_disc = float(np.sqrt(np.mean(np.sum((fit_disc - B) ** 2, axis=1))))

        # COM-optimal continuous translation for this fixed rotation —
        # the theoretical floor (best placement for that rotation).
        t_comfit = B_com - A_rot.mean(axis=0)
        fit_com  = A_rot + t_comfit
        rmsd_com = float(np.sqrt(np.mean(np.sum((fit_com - B) ** 2, axis=1))))

        drot = bf._rot_deviation_deg(R_lsq, R_cart)
        cands.append(dict(R_frac=R_frac, t_frac=t_arr, t_apply=t_apply,
                          R_int=R_int, t_mod=t_mod, R_cart=R_cart,
                          drot=drot, rmsd_discrete=rmsd_disc,
                          rmsd_comfit=rmsd_com, t_comfit_cart=t_comfit))

    # Rank: the LSQ fit anchors which rotation is correct (drot), but the
    # quantity we actually need small is the honest *discrete* RMSD — the
    # residual after applying the exact crystallographic op.  Only an op
    # with a small discrete RMSD can be applied verbatim to both the PDB
    # and the MTZ (SF theorem) and then handed to refmac as usual.
    # comfit RMSD is the theoretical floor (best placement for that fixed
    # rotation) and is reported for context only.
    #
    # Sort primarily by discrete RMSD; a genuinely wrong rotation cannot
    # produce a small discrete RMSD (no single whole-lattice shift can
    # collapse a mis-rotated atom cloud), so the LSQ rotation is implicitly
    # enforced.  drot is carried for cross-checking against the LSQ answer.
    cands.sort(key=lambda c: (c['rmsd_discrete'], round(c['drot'], 3)))

    return dict(mov_pdb=mov_pdb, ref_pdb=ref_pdb, n_ca=len(A),
                rmsd_before=rmsd_before, rmsd_lsq=rmsd_lsq,
                R_lsq=R_lsq, t_lsq=t_lsq, cell=cell, sg_name=sg_name,
                O=O, Finv=Finv, candidates=cands)


def report(res, top_n=15):
    cands = res['candidates']
    print(f"\n  {len(cands)} unique (altindex x symop x origin) candidates")
    print(f"  top {top_n} ranked by honest discrete RMSD "
          f"(LSQ floor = {res['rmsd_lsq']:.2f} A):")
    print(f"  {'rank':>4} {'discrete':>9} {'drot':>7} {'comfit':>8}  "
          f"R (flat)                t (frac)")
    print("  " + "-" * 78)
    for i, c in enumerate(cands[:top_n], 1):
        rflat = ','.join(str(int(x)) for x in c['R_int'].flatten())
        tstr  = ','.join(f"{x:.3f}" for x in c['t_mod'])
        print(f"  {i:>4} {c['rmsd_discrete']:>9.2f} {c['drot']:>7.2f} "
              f"{c['rmsd_comfit']:>8.2f}  {rflat:<22s}  {tstr}")


def apply_and_write(res, choice, out_pdb, mov_mtz=None, out_mtz=None,
                    verbose=True):
    """Apply `choice` (a candidate dict) to the moving PDB (and MTZ).

    The PDB gets R_cart + the exact discrete crystallographic translation
    so it stays consistent with the altindexed MTZ; refmac can then be run
    on the pair as usual.
    """
    O = res['O']
    R_cart = choice['R_cart']
    # t_apply already folds in the whole-lattice image shift; an integer
    # part is invisible to the MTZ phases, so the same t goes to both.
    t_apply = choice['t_apply']
    bf._apply_cart_transform_to_pdb(res['mov_pdb'], R_cart, O @ t_apply,
                                    out_pdb, ref_cell=res['cell'])
    if verbose:
        print(f"\n  wrote {out_pdb}")
        print(f"    R_int = {choice['R_int'].flatten().tolist()}")
        print(f"    t     = {[round(float(x),4) for x in t_apply]}")
        print(f"    residual (honest discrete CA RMSD) = "
              f"{choice['rmsd_discrete']:.2f} A")

    if mov_mtz is not None and out_mtz is not None:
        bf._apply_op_to_mtz(mov_mtz, choice['R_int'], t_apply, out_mtz)
        if verbose:
            print(f"  wrote {out_mtz}  (SF theorem: F'(h)=F(Rᵀh)·e^{{2πi h·t}})")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    mov_pdb = sys.argv[1]
    ref_pdb = sys.argv[2]
    out_pdb = sys.argv[3] if len(sys.argv) > 3 else \
        os.path.splitext(os.path.basename(mov_pdb))[0] + '_altaligned.pdb'
    mov_mtz = sys.argv[4] if len(sys.argv) > 4 else None
    out_mtz = sys.argv[5] if len(sys.argv) > 5 else (
        os.path.splitext(os.path.basename(mov_mtz))[0] + '_altaligned.mtz'
        if mov_mtz else None)

    res = search(mov_pdb, ref_pdb)
    report(res)

    # Default choice: rank-1 by (drot, comfit).  This is the LSQ-based pick.
    choice = res['candidates'][0]
    apply_and_write(res, choice, out_pdb, mov_mtz, out_mtz)


if __name__ == '__main__':
    main()
