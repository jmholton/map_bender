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


def _normalizes_sg(R_int, sg):
    """Does R conjugate every SG symop back into the SG (R·g·R⁻¹ ∈ SG)?
    Non-normalizers land the reindexed crystal in a conjugate setting (e.g.
    obverse→reverse for rhombohedral H groups), which is the obverse/reverse
    problem the dual writer exists to handle.
    """
    Ri = np.round(np.linalg.inv(R_int.astype(float))).astype(int)
    sg_ops = []
    sg_set = set()
    for op in sg.operations():
        Rg = np.round(np.array(op.rot, dtype=float) / op.DEN).astype(int)
        tg = np.array(op.tran, dtype=float) / op.DEN
        sg_ops.append((Rg, tg))
        sg_set.add((tuple(Rg.flatten()), tuple(np.round(tg % 1.0, 4))))
    for Rg, tg in sg_ops:
        Rc = R_int @ Rg @ Ri
        tc = R_int @ tg
        key = (tuple(np.round(Rc).astype(int).flatten()),
               tuple(np.round(tc % 1.0, 4)))
        if key not in sg_set:
            return False
    return True


def _hex_to_R32R_setup(cell_hex):
    """Return (cell_rh, M_obv, sgR) — obverse-hex → R 3 2 :R basis change.

    M_obv is the linear part of gemmi's basisop for R 3 2 :R; it transforms
    fractional coords from obverse-hex into rhombohedral primitive.  No
    translation (basisop has tran=0).  The rhombohedral cell params come
    from the standard formula on hex (a,c).
    """
    a = cell_hex.a
    c = cell_hex.c
    a_rh = float(np.sqrt(3 * a * a + c * c) / 3.0)
    cos_al = (2 * c * c - 3 * a * a) / (2 * (c * c + 3 * a * a))
    al = float(np.degrees(np.arccos(cos_al)))
    cell_rh = gemmi.UnitCell(a_rh, a_rh, a_rh, al, al, al)
    sgR = gemmi.find_spacegroup_by_name("R 3 2:R")
    M_obv = np.array(sgR.basisop.rot, dtype=float) / sgR.basisop.DEN
    return cell_rh, M_obv, sgR


def _is_rhombohedral_in_hex(sg):
    """True for the R-lattice groups expressed in hexagonal axes (ext='H').
    These are the only SGs where the obverse/reverse problem can arise, and
    where R 3 2 :R-style rhombohedral primitive setting is meaningful.
    """
    return getattr(sg, 'ext', '') == 'H'


def apply_and_write(res, choice, out_pdb, mov_mtz=None, out_mtz=None,
                    verbose=True):
    """Apply `choice` to the moving PDB (and MTZ).

    For rhombohedral-in-hexagonal SGs (R3, R-3, R32, R3m, R3c, R-3m, R-3c
    in ext='H' setting) where the chosen op does NOT normalize the SG —
    i.e. it lands the moving crystal in the conjugate (reverse) setting —
    `out_pdb` is treated as a stem and TWO solutions are written:

      <stem>_H32.{pdb,mtz} — hexagonal axes preserved; atoms+data
        reindexed by (R, t).  gemmi cannot name reverse-H so the SYMM
        records gemmi writes for the MTZ stay obverse; the data is in the
        conjugate setting.  Caveat printed; downstream refmac needs the
        SYMM records swapped (e.g. via CCP4 mtzutils) for correct
        symmetry handling.  Reference workflow unchanged.

      <stem>_R32R.{pdb,mtz} — rhombohedral primitive R 3 2 :R via gemmi's
        basisop.  Setting-unambiguous, no centering, refmac-ready.
        Reference must also be converted to R 3 2 :R for downstream
        comparison (or bendfinder taught to bridge the settings).

    For all other cases (normalizing ops, non-rhombohedral SGs), a single
    `out_pdb` (+ `out_mtz`) is written as before.
    """
    sg = gemmi.find_spacegroup_by_name(res['sg_name'])
    is_rhomb_hex = _is_rhombohedral_in_hex(sg)
    is_norm = _normalizes_sg(choice['R_int'], sg)
    emit_dual = is_rhomb_hex and not is_norm

    O = res['O']; Finv = res['Finv']
    R_int = choice['R_int']
    t_apply = choice['t_apply']
    R_cart = choice['R_cart']

    if verbose:
        print(f"\n  chosen op:")
        print(f"    R_int = {R_int.flatten().tolist()}")
        print(f"    t     = {[round(float(x),4) for x in t_apply]}")
        print(f"    residual (honest discrete CA RMSD) = "
              f"{choice['rmsd_discrete']:.2f} A")
        print(f"    normalizes SG = {is_norm}    rhomb-in-hex SG = {is_rhomb_hex}")
        print(f"    emitting {'dual (H32+SYMM and R32:R)' if emit_dual else 'single (op is clean)'}")

    if not emit_dual:
        bf._apply_cart_transform_to_pdb(res['mov_pdb'], R_cart, O @ t_apply,
                                        out_pdb, ref_cell=res['cell'])
        if verbose:
            print(f"  wrote {out_pdb}")
        if mov_mtz is not None and out_mtz is not None:
            bf._apply_op_to_mtz(mov_mtz, R_int, t_apply, out_mtz)
            if verbose:
                print(f"  wrote {out_mtz}")
        return

    # ── dual emission ────────────────────────────────────────────────────
    stem, ext = os.path.splitext(out_pdb)
    H32_pdb  = f"{stem}_H32{ext}"
    R32R_pdb = f"{stem}_R32R{ext}"
    H32_mtz = R32R_mtz = None
    if out_mtz is not None:
        mstem, mext = os.path.splitext(out_mtz)
        H32_mtz  = f"{mstem}_H32{mext}"
        R32R_mtz = f"{mstem}_R32R{mext}"

    # Solution A — H 3 2 + SYMM (reverse SYMM is *not* what gemmi writes).
    # Same operation as the single-output case: cartesian reindex of the PDB,
    # SF-theorem reindex of the MTZ.  Cell and SG label stay H 3 2.
    bf._apply_cart_transform_to_pdb(res['mov_pdb'], R_cart, O @ t_apply,
                                    H32_pdb, ref_cell=res['cell'])
    if mov_mtz is not None and H32_mtz is not None:
        bf._apply_op_to_mtz(mov_mtz, R_int, t_apply, H32_mtz)

    # Solution B — R 3 2 :R rhombohedral primitive.
    # Composite fractional op original-hex → aligned R 3 2 :R:
    #   f_rh = M_rev @ (R @ f_hex + t_apply)
    # where M_rev = M_obv @ diag(-1,-1,1) is the reverse-hex → rhomb basis
    # change.  The diag(-1,-1,1) factor accounts for the obverse↔reverse
    # axis flip that the non-normalizing reindex introduces — M_obv alone
    # would mis-route reverse-setting fractional/HKLs.
    cell_hex = res['cell']
    cell_rh, M_obv, sgR = _hex_to_R32R_setup(cell_hex)
    M_rev = M_obv @ np.diag([-1.0, -1.0, 1.0])
    M_total = M_rev @ R_int.astype(float)
    t_total = M_rev @ np.asarray(t_apply, dtype=float)

    st = gemmi.read_structure(res['mov_pdb'])
    out = st.clone()
    out.cell = cell_rh
    out.spacegroup_hm = "R 3 2:R"
    for model in out:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    fpos = cell_hex.fractionalize(atom.pos)
                    f_hex = np.array([fpos.x, fpos.y, fpos.z])
                    f_rh = M_total @ f_hex + t_total
                    p = cell_rh.orthogonalize(gemmi.Fractional(*f_rh))
                    atom.pos = gemmi.Position(p.x, p.y, p.z)
    out.setup_entities()
    out.write_pdb(R32R_pdb)

    if mov_mtz is not None and R32R_mtz is not None:
        # HKL transform: h_rh = h_hex @ M_total^-1  (row-vector convention).
        # |det M_total| = 3 so 2/3 of input rows give non-integer h_rh —
        # these are exactly the H-centering systematic absences in the
        # original obverse-hex labeling (which become non-absent in the
        # reverse-hex labeling that the reindex produces — but in the
        # rhomb primitive description they collapse correctly).
        mtz_in = gemmi.read_mtz_file(mov_mtz)
        data = np.array(mtz_in.array, copy=True)
        H_hex = data[:, :3].astype(int)
        M_total_inv = np.linalg.inv(M_total)
        H_rh_f = H_hex.astype(float) @ M_total_inv
        H_rh = np.round(H_rh_f)
        is_int = np.all(np.abs(H_rh_f - H_rh) < 1e-3, axis=1)
        # Phase shift from the reindex step alone (R32:R basisop has no
        # translation).  The phase shift is 360° · h_revhex · t_apply, where
        # h_revhex = h_hex @ R_int^-1 (the HKL after just the reindex).
        R_inv = np.round(np.linalg.inv(R_int.astype(float))).astype(int)
        H_revhex = H_hex @ R_inv
        delta_phi_all = (360.0 * (H_revhex.astype(np.float64)
                                  @ np.asarray(t_apply, dtype=np.float64))
                         ).astype(np.float32)

        data_out = data[is_int].copy()
        data_out[:, :3] = H_rh[is_int].astype(np.float32)
        # Apply phase shifts to all PHI columns
        cols = list(mtz_in.columns)
        for i, c in enumerate(cols):
            if c.type == 'P':
                data_out[:, i] = (data_out[:, i] + delta_phi_all[is_int]) % 360.0

        mtz_in.cell = cell_rh
        mtz_in.spacegroup = sgR
        # CCP4 reads cells from the per-dataset records, not the file-level
        # cell; failing to update them gives refmac "Large differences
        # between cells from pdb and mtz" even when the file-level cell is
        # right.  Same for the file-level / per-dataset wavelength.
        for ds in mtz_in.datasets:
            ds.cell = cell_rh
        mtz_in.set_data(data_out)
        mtz_in.write_to_file(R32R_mtz)

    if verbose:
        print(f"\n  wrote H 3 2 (hex axes) solution:")
        print(f"    PDB: {H32_pdb}")
        if H32_mtz:
            print(f"    MTZ: {H32_mtz}")
            print(f"    NOTE: data is in the conjugate (reverse-H) setting; gemmi cannot")
            print(f"          label that, so the MTZ's SYMM records remain obverse.  For")
            print(f"          refmac on this output, swap SYMM with CCP4 mtzutils first.")
        print(f"  wrote R 3 2 :R (rhomb primitive) solution:")
        print(f"    PDB: {R32R_pdb}")
        if R32R_mtz:
            print(f"    MTZ: {R32R_mtz}")
            print(f"    Setting-unambiguous, refmac-ready.  Reference (3pou) must also be")
            print(f"          converted to R 3 2 :R for direct comparison, since the rhomb")
            print(f"          cell's gemmi-canonical cartesian frame differs from the hex frame.")


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
