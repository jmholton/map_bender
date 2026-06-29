#!/usr/bin/env ccp4-python
# -*- coding: utf-8 -*-
"""altalign.py - standalone alternate-indexing + origin search.

Finds the crystallographic (altindex x symop x origin) transformation that
best aligns a moving crystal onto a reference crystal.  Mirrors the same
pipeline bendfinder's internal `resolve_altindex` uses, so a standalone
altalign run picks the same op the full fitreso_scan would have picked:

  1. Cell-stretch pre-step.  If the moving and reference cells differ at
     all, re-orthogonalize the moving model into the reference cell with
     fractional coordinates preserved (and relabel the MTZ cell to match).
     This absorbs the isomorphous distortion as a uniform elastic stretch
     so the discrete enumeration sees a metric-matched pair.
  2. Kabsch rigid-body fit of matched CA atoms -> continuous R_lsq.
  3. Enumerate discrete (altindex x symop x origin) candidates via the
     basis-change altindex enumeration (_get_altindex_ops) + the SG's own
     symops + an origin grid.  Strict 1e-6 metric-preservation tolerance
     for same-cell pairs; loose 5% tolerance after a cell stretch to
     surface ops that hold in a nearby higher-symmetry holohedry the
     deposited cell only approximates.
  4. Score every candidate two ways, both in honest cartesian RMSD:
       - continuous: COM-optimal continuous t (the theoretical floor for
                     that rotation).  This is what bendfinder ranks by.
       - discrete:   the exact crystallographic translation.  Carried as
                     a diagnostic — if the discrete RMSD is close to the
                     continuous floor, the op is clean as-is; if not, the
                     remaining offset is sub-cell continuous and bendfinder
                     would absorb it via the COM-optimal t_cart_opt.
  5. Action gate on the rank-1 candidate (matches resolve_altindex):
       - none           : no candidate beats 0.7 x baseline
       - origin_only    : R = I, |t| > 0.01 (PDB translate + MTZ phase shift)
       - sg_op_origin   : R is in the SG (PDB+MTZ via SF theorem; no refmac)
       - altindex_refine: R != I and R not in SG (PDB cartesian transform
                          + MTZ reindex by R; --refine optionally chains
                          run_refinement; for non-normalizing R-lattice
                          ops also emits R32:R dual output)

Usage:
    ccp4-python altalign.py <moving.pdb> <reference.pdb> \
        [out.pdb] [moving.mtz] [out.mtz] [--refine] [--outdir DIR]

If moving.mtz is given, the chosen op is applied to it so the PDB and MTZ
stay mutually consistent.  --refine chains refmac5 after altindex_refine.
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


def search(mov_pdb, ref_pdb, mov_mtz=None, outdir='.',
           max_M=1, det_max=1, verbose=True):
    """Run the LSQ-based altindex+origin search.

    Mirrors bendfinder's resolve_altindex enumeration semantics: cells
    are matched via a stretch pre-step when they differ, and the metric
    tolerance is loosened from 1e-6 to 5% post-stretch so the
    same-habit-distorted-cell case (lipox-style) surfaces alt-cell ops
    that strict tolerance would reject.

    Returns a dict with the LSQ baseline and a ranked candidate list.
    Each candidate is a dict: R_frac, t_frac, t_apply, t_opt_frac,
    R_int, t_mod, R_cart, drot, rmsd_discrete, rmsd_continuous.

    If a cell stretch was applied, `mov_pdb` and `mov_mtz` in the
    returned dict point at the stretched copies in `outdir`; downstream
    apply steps consume those rather than the originals.
    """
    mov_st0 = gemmi.read_structure(mov_pdb)
    ref_st  = gemmi.read_structure(ref_pdb)

    # ── cell-stretch pre-step ─────────────────────────────────────────
    stretched = False
    if not bf._cells_match(mov_st0.cell, ref_st.cell):
        os.makedirs(outdir, exist_ok=True)
        s_pdb_stem = os.path.splitext(os.path.basename(mov_pdb))[0]
        s_pdb = os.path.join(outdir, f'{s_pdb_stem}_stretched.pdb')
        if verbose:
            print(f"altalign: cells differ — stretching moving into "
                  f"ref cell (fractions preserved)", flush=True)
            print(f"  moving cell : {tuple(round(x, 3) for x in mov_st0.cell.parameters)}")
            print(f"  ref    cell : {tuple(round(x, 3) for x in ref_st.cell.parameters)}",
                  flush=True)
        bf._stretch_pdb_to_cell(mov_pdb, ref_st.cell, s_pdb)
        mov_pdb = s_pdb
        if mov_mtz is not None:
            s_mtz_stem = os.path.splitext(os.path.basename(mov_mtz))[0]
            s_mtz = os.path.join(outdir, f'{s_mtz_stem}_stretched.mtz')
            bf._relabel_mtz_cell(mov_mtz, ref_st.cell, s_mtz)
            mov_mtz = s_mtz
        stretched = True

    mov_st = gemmi.read_structure(mov_pdb)
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
    A_frac = (Finv @ A.T).T
    B_frac = (Finv @ B.T).T
    metric_tol_rel = 0.05 if stretched else 1e-6

    if verbose:
        print(f"altalign: {mov_pdb}  ->  {ref_pdb}")
        print(f"  matched CA      : {len(A)}")
        print(f"  RMSD before fit : {rmsd_before:.2f} A")
        print(f"  LSQ rigid fit   : {rmsd_lsq:.2f} A   "
              f"(rotation {_rot_angle_deg(R_lsq):.1f} deg)")
        print(f"  metric tol      : {metric_tol_rel:g} "
              f"({'loose post-stretch' if stretched else 'strict same-cell'})")
        print(f"  enumerating altindex x symop x origin "
              f"(max_M={max_M}, det_max={det_max}) ...", flush=True)

    # Deduplicate (R_frac, t_frac mod 1) candidates from the enumerator.
    seen = set()
    cands = []
    B_com = B.mean(axis=0)
    for R_frac, t_frac in bf._enum_alt_rot_origin_candidates(
            cell, sg, sg_name, metric_tol_rel=metric_tol_rel):
        if not bf._is_metric_preserving(R_frac, G, atol=metric_tol_rel):
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

        # COM-optimal continuous translation in fractional space (the
        # quantity bendfinder ranks on and applies as t_cart_opt).  This
        # is the theoretical floor for the given rotation.
        A_rot_frac = (R_frac @ A_frac.T).T
        t_opt_frac = B_frac.mean(axis=0) - A_rot_frac.mean(axis=0)
        diff_frac  = (A_rot_frac + t_opt_frac) - B_frac
        cart_diff  = (O @ diff_frac.T).T
        rmsd_cont  = float(np.sqrt(np.mean(np.sum(cart_diff ** 2, axis=1))))

        drot = bf._rot_deviation_deg(R_lsq, R_cart)
        cands.append(dict(R_frac=R_frac, t_frac=t_arr, t_apply=t_apply,
                          t_opt_frac=t_opt_frac, R_int=R_int, t_mod=t_mod,
                          R_cart=R_cart, drot=drot,
                          rmsd_discrete=rmsd_disc,
                          rmsd_continuous=rmsd_cont))

    # Rank by the COM-optimal continuous RMSD — same metric bendfinder's
    # resolve_altindex picks on, so a standalone altalign run lands on
    # the same candidate the full fitreso_scan would.  Discrete RMSD is
    # carried as a diagnostic: if it's close to the continuous floor the
    # op is clean as-is; if not, the extra offset is sub-cell continuous
    # and gets folded into t_cart_opt at apply time.  drot vs the LSQ
    # rigid-body answer is a cross-check; an op with the wrong rotation
    # can't reach a small continuous RMSD anyway, so drot is mostly
    # informative not gating.
    cands.sort(key=lambda c: (c['rmsd_continuous'], round(c['drot'], 3)))

    return dict(mov_pdb=mov_pdb, mov_mtz=mov_mtz, ref_pdb=ref_pdb,
                stretched=stretched, n_ca=len(A),
                rmsd_before=rmsd_before, rmsd_lsq=rmsd_lsq,
                R_lsq=R_lsq, t_lsq=t_lsq, cell=cell, sg_name=sg_name,
                O=O, Finv=Finv, candidates=cands)


def report(res, top_n=15):
    cands = res['candidates']
    print(f"\n  {len(cands)} unique (altindex x symop x origin) candidates")
    print(f"  top {top_n} ranked by COM-optimal continuous RMSD "
          f"(LSQ floor = {res['rmsd_lsq']:.2f} A):")
    print(f"  {'rank':>4} {'cont':>8} {'discrete':>9} {'drot':>7}  "
          f"R (flat)                t (frac)")
    print("  " + "-" * 78)
    for i, c in enumerate(cands[:top_n], 1):
        rflat = ','.join(str(int(x)) for x in c['R_int'].flatten())
        tstr  = ','.join(f"{x:.3f}" for x in c['t_mod'])
        print(f"  {i:>4} {c['rmsd_continuous']:>8.2f} "
              f"{c['rmsd_discrete']:>9.2f} {c['drot']:>7.2f}  "
              f"{rflat:<22s}  {tstr}")


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


def _dual_write_R32R(res, choice, t_apply, out_pdb, mov_mtz, out_mtz, verbose):
    """Emit dual H32+R32R outputs for non-normalizing R-lattice altindex ops.

    See the obverse/reverse discussion in CLAUDE.md / README.  Stem-name
    convention: out_pdb -> <stem>_H32.{pdb} + <stem>_R32R.{pdb}; same for
    MTZ.  Both outputs apply (R, t_apply) by the SF theorem so PDB+MTZ
    stay mutually consistent without needing refmac.

    `t_apply` is supplied by the caller (typically the COM-optimal
    continuous t_frac picked by resolve()) — the SF theorem accepts any
    real t.  Post-stretch the discrete-grid entries from the enumerator
    are off by sub-cell amounts, so using the continuous t keeps the
    H32 and R32:R outputs self-consistent with the PDB.

    Caller is responsible for checking that the chosen op is in fact
    non-normalizing in an R-lattice (ext='H') SG before invoking this.
    """
    O = res['O']
    R_int = choice['R_int']
    R_cart = choice['R_cart']

    stem, ext = os.path.splitext(out_pdb)
    H32_pdb  = f"{stem}_H32{ext}"
    R32R_pdb = f"{stem}_R32R{ext}"
    H32_mtz = R32R_mtz = None
    if out_mtz is not None:
        mstem, mext = os.path.splitext(out_mtz)
        H32_mtz  = f"{mstem}_H32{mext}"
        R32R_mtz = f"{mstem}_R32R{mext}"

    # Solution A — H 3 2 + SYMM (reverse SYMM is *not* what gemmi writes).
    # Cartesian reindex of the PDB, SF-theorem reindex of the MTZ.  Cell
    # and SG label stay H 3 2.
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


def resolve(res, out_pdb, out_mtz=None, outdir='.', refine=False,
            refine_cycles=5, fill_asu=False, improve_threshold=0.7,
            verbose=True):
    """Apply the rank-1 candidate using bendfinder's 4-action gate.

    Returns a dict with action ('none' / 'origin_only' / 'sg_op_origin'
    / 'altindex_refine'), the chosen op, and the output file paths
    that were written (or None for action='none').

    For non-normalizing R-lattice (ext='H') altindex_refine ops, also
    emits dual H32+R32R outputs (file paths returned in `dual_outputs`).

    If `refine` is True and the action is altindex_refine, chain
    bendfinder.run_refinement on the reindexed pair.
    """
    mov_pdb = res['mov_pdb']
    mov_mtz = res['mov_mtz']
    choice  = res['candidates'][0]
    cell    = res['cell']
    sg      = gemmi.find_spacegroup_by_name(res['sg_name'])
    O       = res['O']
    R_frac  = choice['R_frac']
    R_int   = choice['R_int']
    R_cart  = choice['R_cart']
    t_frac  = choice['t_opt_frac']        # COM-optimal continuous t (bendfinder's pick)
    t_disc  = choice['t_apply']           # exact crystallographic t (carried for dual writer)
    rmsd_top  = choice['rmsd_continuous']
    rmsd_base = res['rmsd_before']

    # Action gate (same conditions as bendfinder.resolve_altindex).
    is_identity_R = np.allclose(R_frac, np.eye(3), atol=1e-6)
    is_zero_t     = float(np.linalg.norm(t_frac)) < 0.01
    R_in_sg = False
    if np.allclose(R_frac, R_int, atol=1e-4):
        for op in sg.operations():
            sg_R = (np.array(op.rot, dtype=int) // op.DEN)
            if np.array_equal(sg_R, R_int):
                R_in_sg = True
                break

    def _say(*a, **k):
        if verbose: print(*a, **k, flush=True)

    _say(f"\n  chosen op (rank 1):")
    _say(f"    R_int   = {R_int.flatten().tolist()}")
    _say(f"    t_cont  = {[round(float(x),4) for x in t_frac]}  (COM-optimal continuous)")
    _say(f"    t_disc  = {[round(float(x),4) for x in t_disc]}  (exact crystallographic)")
    _say(f"    rmsd    = {rmsd_top:.2f} A continuous  ({choice['rmsd_discrete']:.2f} A discrete)")
    _say(f"    baseline= {rmsd_base:.2f} A  (LSQ floor {res['rmsd_lsq']:.2f} A)")

    if rmsd_top >= improve_threshold * rmsd_base:
        _say(f"  action=none — rank-1 op fails to beat {improve_threshold:g}x baseline")
        _say(f"  (no output written; moving model already aligned within scope of discrete search)")
        return dict(action='none', choice=choice, out_pdb=None, out_mtz=None,
                    dual_outputs=None)

    if is_identity_R and is_zero_t:
        _say(f"  action=none — best op is identity, no transform needed")
        return dict(action='none', choice=choice, out_pdb=None, out_mtz=None,
                    dual_outputs=None)

    # ── sg_op_origin: R is already a symop of the SG ─────────────────────
    if R_in_sg and not is_identity_R:
        _say(f"  action=sg_op_origin  (R is an in-SG symop; full SF theorem on MTZ)")
        bf._apply_cart_transform_to_pdb(mov_pdb, R_cart, O @ t_frac,
                                        out_pdb, ref_cell=cell)
        _say(f"    wrote {out_pdb}")
        if mov_mtz is not None and out_mtz is not None:
            bf._apply_op_to_mtz(mov_mtz, R_int, t_frac, out_mtz)
            _say(f"    wrote {out_mtz}")
        return dict(action='sg_op_origin', choice=choice,
                    out_pdb=out_pdb, out_mtz=out_mtz, dual_outputs=None)

    # ── origin_only: R = I, translate PDB and phase-shift MTZ ────────────
    if is_identity_R:
        _say(f"  action=origin_only  (R = I; phase shift on MTZ, no reindex)")
        bf._apply_cart_transform_to_pdb(mov_pdb, R_cart, O @ t_frac,
                                        out_pdb, ref_cell=cell)
        _say(f"    wrote {out_pdb}")
        if mov_mtz is not None and out_mtz is not None:
            bf._shift_mtz_origin(mov_mtz, t_frac, out_mtz)
            _say(f"    wrote {out_mtz}")
        return dict(action='origin_only', choice=choice,
                    out_pdb=out_pdb, out_mtz=out_mtz, dual_outputs=None)

    # ── altindex_refine: R != I and R not in SG ──────────────────────────
    is_rhomb_hex = _is_rhombohedral_in_hex(sg)
    is_norm      = _normalizes_sg(R_int, sg)
    emit_dual    = is_rhomb_hex and not is_norm
    _say(f"  action=altindex_refine  (R is genuine altindex op)")
    _say(f"    normalizes SG = {is_norm}    rhomb-in-hex SG = {is_rhomb_hex}")

    # PDB: apply cartesian R + continuous t_cart_opt (matches bendfinder).
    # MTZ: full SF theorem with the same continuous t — keeps PDB+MTZ
    # mutually consistent without re-refinement.  bendfinder.resolve_altindex
    # uses reindex-only on the MTZ and relies on a subsequent refmac call
    # to regenerate map coefficients; altalign does not assume refmac will
    # follow so the phase shift is applied here.  Post-stretch the discrete
    # enumerator t entries are sub-cell off from the COM-optimal continuous
    # answer, so using the continuous t avoids landing the MTZ at a
    # different sub-cell origin than the PDB.  If --refine is on the
    # SF-theorem MTZ becomes the refmac input; rigid-body cycles absorb
    # any small residual.
    bf._apply_cart_transform_to_pdb(mov_pdb, R_cart, O @ t_frac,
                                    out_pdb, ref_cell=cell)
    if mov_mtz is not None and out_mtz is not None:
        bf._apply_op_to_mtz(mov_mtz, R_int, t_frac, out_mtz)

    dual_outputs = None
    if emit_dual:
        _say(f"    emitting dual (H32+SYMM and R32:R) — obverse/reverse case")
        _dual_write_R32R(res, choice, t_frac, out_pdb, mov_mtz, out_mtz,
                         verbose=verbose)
        stem, ext = os.path.splitext(out_pdb)
        dual_outputs = dict(H32_pdb=f"{stem}_H32{ext}",
                            R32R_pdb=f"{stem}_R32R{ext}")
        if out_mtz is not None:
            mstem, mext = os.path.splitext(out_mtz)
            dual_outputs['H32_mtz']  = f"{mstem}_H32{mext}"
            dual_outputs['R32R_mtz'] = f"{mstem}_R32R{mext}"
    else:
        _say(f"    wrote {out_pdb}")
        if out_mtz is not None:
            _say(f"    wrote {out_mtz}")

    refined = None
    if refine and out_mtz is not None:
        _say(f"  --refine: chaining refmac5 rigid-body ({refine_cycles} cycles)...")
        refined_pdb, refined_mtz = bf.run_refinement(
            out_pdb, out_mtz, outdir=outdir,
            n_cycles=refine_cycles, fill_asu=fill_asu)
        refined = dict(pdb=refined_pdb, mtz=refined_mtz)
        _say(f"    refined PDB: {refined_pdb}")
        _say(f"    refined MTZ: {refined_mtz}")

    return dict(action='altindex_refine', choice=choice,
                out_pdb=out_pdb, out_mtz=out_mtz,
                dual_outputs=dual_outputs, refined=refined)


def _parse_args(argv):
    """Minimal CLI parser: positional args (mov, ref, out_pdb, mov_mtz,
    out_mtz) + --refine, --outdir, --no-fill-asu flags."""
    pos = []
    refine = False
    outdir = '.'
    fill_asu = True
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--refine':
            refine = True
        elif a == '--no-refine':
            refine = False
        elif a == '--outdir':
            outdir = argv[i+1]; i += 1
        elif a == '--no-fill-asu':
            fill_asu = False
        elif a == '--fill-asu':
            fill_asu = True
        elif a in ('-h', '--help'):
            sys.exit(__doc__)
        else:
            pos.append(a)
        i += 1
    return pos, refine, outdir, fill_asu


def main():
    pos, refine, outdir, fill_asu = _parse_args(sys.argv[1:])
    if len(pos) < 2:
        sys.exit(__doc__)
    mov_pdb = pos[0]
    ref_pdb = pos[1]
    out_pdb = pos[2] if len(pos) > 2 else \
        os.path.splitext(os.path.basename(mov_pdb))[0] + '_altaligned.pdb'
    mov_mtz = pos[3] if len(pos) > 3 else None
    out_mtz = pos[4] if len(pos) > 4 else (
        os.path.splitext(os.path.basename(mov_mtz))[0] + '_altaligned.mtz'
        if mov_mtz else None)

    # search() may write stretched copies into outdir; out_pdb/out_mtz from
    # the user describe where the *final* aligned outputs go (resolve()
    # applies the chosen op to the post-stretch model and writes to those).
    res = search(mov_pdb, ref_pdb, mov_mtz=mov_mtz, outdir=outdir)
    report(res)
    resolve(res, out_pdb, out_mtz=out_mtz, outdir=outdir,
            refine=refine, fill_asu=fill_asu)


if __name__ == '__main__':
    main()
