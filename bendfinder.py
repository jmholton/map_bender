#!/programs/ccp4-8.0/bin/ccp4-python
"""bendfinder.py — fit a Fourier shift field between two crystal forms.

Usage (single fit):
    bendfinder.py moving.pdb reference.pdb moving.mtz reference.mtz [key=value ...]
    bendfinder.py moving.pdb reference.pdb [map.ccp4] [key=value ...]   # CCP4 map fallback

Usage (fitreso scan):
    bendfinder.py moving.pdb reference.pdb moving.mtz reference.mtz scan_dir=DIR [key=value ...]

MTZ inputs are expected to be refinement outputs with likelihood-weighted 2Fo-Fc
coefficients (FWT/PHWT from refmac or 2FOFCWT/PH2FOFCWT from phenix).

Parameters:
    nhkls=30          max Fourier terms (default 30); reports effective fitreso
    fitreso=X         use all HKLs with d ≥ X Å (e.g. fitreso=12); overrides nhkls
    drop_snr=1        post-fit SNR threshold (default 1; 0 = disable)
    frac=1            scale factor for applied shifts (default 1)
    geotest=false     run refmac5 geometry check (default false)
    use_symm=true     enforce space-group symmetry constraint (default true)
    dimensions=xyz    which shifts to fit (subset of xyz; default xyz)
    fitparams=file    load pre-computed psdvf.mtz and skip fitting
    deltamaps         write delta-x/y/z/r maps in addition to bent map
    labels=F,PHI      override 2Fo-Fc column detection (e.g. labels=FC,PHIC)
    run_refinement    run refmac5/phenix.refine to generate FWT/PHWT first
    refine_cycles=5   number of refinement cycles (default 5)
    sample_rate=0.0   FFT oversampling for MTZ → map (default 0.0 = gemmi auto)
    scan_dir=DIR      run full fitreso scan and write outputs to DIR
    subtract=ref      diff = bent − ref (default; positive peak = density in
                      bent absent from ref).  subtract=bent flips the sign.
    scan_all_fr=false when False (default) the fitreso scan stops early once
                      the combined `score` column has risen for early_stop_n
                      consecutive rows past its running argmin (early_stop_tol
                      relative tolerance).  scan_all_fr=true keeps all 8
                      fr-rows (fr20..fr5) — reference behavior for the gamut.
    early_stop_tol=0.01  early-stop relative tolerance around the running
                      argmin score (default 1% = 0.01).  A row counts as
                      "worse" only if `row.score > argmin.score * (1+tol)`.
    early_stop_n=2    minimum consecutive worse rows to trigger early-stop
                      (default 2 — one row's rise could be noise).

Public API:
    from bendfinder import bend_fit, bend_apply_pdb, bend_apply_map
    from bendfinder import fitreso_scan, detect_2fofc_cols, mtz_to_map_data
    result = bend_fit("moving.pdb", "reference.pdb", nhkls=30)
    print(result.rmsd)
"""

import sys
import os
import re
import time
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import lstsq as sp_lstsq, svd as sp_svd
from scipy.ndimage import map_coordinates
from scipy.special import erf as sp_erf, erfinv as sp_erfinv, gammainc as sp_gammainc, gammaincinv as sp_gammaincinv, betainc as sp_betainc, betaincinv as sp_betaincinv
import gemmi

TWO_PI = 2.0 * np.pi

_BACKBONE = frozenset({'N', 'CA', 'C', 'O'})


# ══════════════════════════════════════════════════════════════════════════════
# Unit cell math
# ══════════════════════════════════════════════════════════════════════════════

def _ortho_matrix(cell):
    """Orthogonalization matrix M: orth = M @ frac (CCP4 convention 1)."""
    a, b, c, al, be, ga = cell
    al, be, ga = np.radians([al, be, ga])
    cg, sg = np.cos(ga), np.sin(ga)
    cb, ca = np.cos(be), np.cos(al)
    v = np.sqrt(max(1 - ca**2 - cb**2 - cg**2 + 2*ca*cb*cg, 0.0))
    return np.array([
        [a,  b*cg,  c*cb          ],
        [0,  b*sg,  c*(ca-cb*cg)/sg],
        [0,  0,     c*v/sg         ],
    ])


def frac_to_orth(frac_xyz, cell):
    """Fractional → orthogonal Å.  frac_xyz shape (..., 3)."""
    return frac_xyz @ _ortho_matrix(cell).T


def orth_to_frac(orth_xyz, cell):
    return orth_xyz @ np.linalg.inv(_ortho_matrix(cell)).T


def _reciprocal_metric(cell):
    """G* (reciprocal metric tensor): 1/d² = H·G*·Hᵀ."""
    a, b, c, al, be, ga = cell
    al, be, ga = np.radians([al, be, ga])
    ca, cb, cg = np.cos(al), np.cos(be), np.cos(ga)
    sa, sb, sg = np.sin(al), np.sin(be), np.sin(ga)
    V = a*b*c * np.sqrt(max(1 - ca**2 - cb**2 - cg**2 + 2*ca*cb*cg, 0.0))
    a_s = b*c*sa/V
    b_s = a*c*sb/V
    c_s = a*b*sg/V
    ca_s = (cb*cg - ca) / (sb*sg)
    cb_s = (ca*cg - cb) / (sa*sg)
    cg_s = (ca*cb - cg) / (sa*sb)
    return np.array([
        [a_s**2,        a_s*b_s*cg_s,   a_s*c_s*cb_s],
        [a_s*b_s*cg_s,  b_s**2,         b_s*c_s*ca_s],
        [a_s*c_s*cb_s,  b_s*c_s*ca_s,   c_s**2      ],
    ])


# ══════════════════════════════════════════════════════════════════════════════
# PDB I/O and P1 expansion
# ══════════════════════════════════════════════════════════════════════════════

def _cell_tuple(gc):
    return (gc.a, gc.b, gc.c, gc.alpha, gc.beta, gc.gamma)


def _altloc_spread_tol(st, fallback=1.0):
    """Derive altloc tolerance from backbone atoms that have multiple conformers.

    For each backbone atom (N/CA/C/O) with ≥2 altloc positions, compute the
    max distance of any position from the centroid.  Return the median of
    those spreads if ≥5 such atoms are found, else return fallback (Å).
    """
    spreads = []
    for model in st:
        for chain in model:
            for res in chain:
                if res.name in ('HOH', 'WAT', 'H2O'):
                    continue
                by_name = {}
                for atom in res:
                    if atom.element == gemmi.Element('H'):
                        continue
                    aname = atom.name.strip()
                    if aname not in _BACKBONE:
                        continue
                    orth = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    by_name.setdefault(aname, []).append(orth)
                for positions in by_name.values():
                    if len(positions) < 2:
                        continue
                    pos_arr = np.array(positions)
                    centroid = pos_arr.mean(axis=0)
                    spread = float(np.max(np.linalg.norm(pos_arr - centroid, axis=1)))
                    spreads.append(spread)
    if len(spreads) >= 5:
        return float(np.median(spreads)), len(spreads)
    return fallback, len(spreads)


def reject_outliers(fitme, ca_mask, uids, cell, mad_sigma=3.0,
                    bfacs=None, b_sigma=None, verbose=True):
    """Remove atom pairs whose Å shift magnitude is an outlier (MAD filter).

    Optionally also reject atoms with outlier B-factors (max(B1,B2) > MAD threshold).

    Returns keep_mask (N,) bool alongside filtered arrays.
    """
    shift_orth = frac_to_orth(fitme[:, 3:6], cell)
    mags = np.linalg.norm(shift_orth, axis=1)
    med = np.median(mags)
    mad = np.median(np.abs(mags - med))
    sigma_eq = mad / 0.6745          # MAD → equivalent normal sigma
    tol = med + mad_sigma * sigma_eq
    keep = mags <= tol
    n_drop_shift = int((~keep).sum())
    if n_drop_shift and verbose:
        print(f"  shift outlier rejection: dropped {n_drop_shift}/{len(mags)} atoms "
              f"(|Δr|>{tol:.3f} Å; median={med:.3f}, σ_MAD={sigma_eq:.3f} Å)")

    if bfacs is not None and b_sigma is not None:
        bmax = np.maximum(bfacs[:, 0], bfacs[:, 1])
        b_med = np.median(bmax)
        b_mad = np.median(np.abs(bmax - b_med))
        b_sig_eq = b_mad / 0.6745
        b_tol = b_med + b_sigma * b_sig_eq
        b_keep = bmax <= b_tol
        n_drop_b = int((~b_keep & keep).sum())
        if n_drop_b and verbose:
            print(f"  B-factor outlier rejection: dropped additional {n_drop_b} atoms "
                  f"(B>{b_tol:.1f} Å²; median={b_med:.1f}, σ_MAD={b_sig_eq:.1f} Å²)")
        keep = keep & b_keep

    return (fitme[keep], ca_mask[keep],
            [u for u, k in zip(uids, keep) if k],
            bfacs[keep] if bfacs is not None else None,
            keep)


def expand_to_p1(pdb_path, altloc_filter=False, altloc_fallback=1.0, origin_shift=None):
    """Read PDB and expand ASU to P1 with all crystallographic symmetry ops.

    Parameters
    ----------
    altloc_filter : bool
        If True, collect all altloc positions for each atom and average them
        when the max spread from centroid is within the derived tolerance;
        atoms whose altlocs spread more than the tolerance are dropped.
        The tolerance is the median backbone-altloc spread if ≥5 such atoms
        are available, otherwise altloc_fallback Å.
        If False (default), only altloc 'A' (or no-altloc) atoms are kept.

    Returns
    -------
    atoms : list[dict]
        Keys: x y z bfac occ chain resnum icode resname atomname is_ca uid
        x,y,z are fractional in the structure's own cell, wrapped to [0,1).
        uid is unique per symmetry copy.
    cell : tuple (a,b,c,alpha,beta,gamma)
    sg   : str — space group name
    """
    st = gemmi.read_pdb(str(pdb_path))
    gc = st.cell
    cell = _cell_tuple(gc)
    M_inv = np.linalg.inv(_ortho_matrix(cell))

    sg_obj = st.find_spacegroup()
    if sg_obj is None:
        sg_obj = gemmi.find_spacegroup_by_name('P 1')
    sg_name = sg_obj.short_name()   # 'H32' not 'R 3 2' for hexagonal setting

    ops = list(sg_obj.operations())   # identity is ops[0]

    # Derive altloc tolerance once from backbone disorder, if requested
    if altloc_filter:
        altloc_tol, n_bb = _altloc_spread_tol(st, fallback=altloc_fallback)
        src = f'backbone median ({n_bb} atoms)' if n_bb >= 5 else f'fallback ({n_bb} backbone altlocs)'
        print(f"  altloc tol={altloc_tol:.3f} Å [{src}]")
    else:
        altloc_tol = None

    atoms = []
    for oi, op in enumerate(ops):
        for model in st:
            for chain in model:
                for res in chain:
                    if res.name in ('HOH', 'WAT', 'H2O'):
                        continue
                    rname = res.name.strip()
                    cname = chain.name.strip()
                    rnum  = res.seqid.num
                    icode = res.seqid.icode.strip()

                    if altloc_tol is not None:
                        # Collect all altloc positions per atom name, then filter
                        by_name = {}
                        for atom in res:
                            if atom.element == gemmi.Element('H'):
                                continue
                            aname = atom.name.strip()
                            orth = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                            d = by_name.setdefault(aname,
                                                   {'pos': [], 'bfac': [], 'occ': []})
                            d['pos'].append(orth)
                            d['bfac'].append(atom.b_iso)
                            d['occ'].append(atom.occ)
                        atom_iter = []
                        for aname, d in by_name.items():
                            pos_arr = np.array(d['pos'])
                            centroid = pos_arr.mean(axis=0)
                            spread = float(np.max(
                                np.linalg.norm(pos_arr - centroid, axis=1)))
                            if spread > altloc_tol:
                                continue   # incompatible conformers — skip
                            atom_iter.append((aname, centroid,
                                              float(np.mean(d['bfac'])),
                                              float(np.mean(d['occ']))))
                    else:
                        # Default: keep only no-altloc or altloc-A atoms
                        atom_iter = []
                        for atom in res:
                            if atom.element == gemmi.Element('H'):
                                continue
                            if atom.altloc not in ('\x00', 'A', ' ', '\0'):
                                continue
                            aname = atom.name.strip()
                            orth = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                            atom_iter.append((aname, orth, atom.b_iso, atom.occ))

                    for aname, orth, bfac, occ in atom_iter:
                        frac = M_inv @ orth
                        if origin_shift is not None:
                            frac = frac + origin_shift   # shift ASU before applying op
                        nf = op.apply_to_xyz(frac.tolist())
                        nf = [f % 1.0 for f in nf]
                        uid = f"{cname}_{rnum}{icode}_{rname}_{aname}_op{oi}"
                        atoms.append({
                            'x': nf[0], 'y': nf[1], 'z': nf[2],
                            'bfac': bfac, 'occ': occ,
                            'chain': cname, 'resnum': rnum, 'icode': icode,
                            'resname': rname, 'atomname': aname,
                            'is_ca': (aname == 'CA'),
                            'uid': uid,
                        })
    return atoms, cell, sg_name


def _maybe_remap_single_chain(atoms1, atoms2, seq_id_min=0.9, verbose=True):
    """If atoms1 and atoms2 are each single-chain but use different chain
    IDs (e.g. 'A' vs 'X' — common with older PDB-REDO outputs), return a
    chain-remapped *copy* of atoms2 with its chain renamed to match
    atoms1's, gated by a sequence-similarity sanity check.

    Sanity: among residues sharing (resnum, icode) in both inputs, at
    least `seq_id_min` fraction must have the same residue name.  This
    keeps the remap from silently aligning two different proteins that
    happen to overlap in residue numbering.

    Returns atoms2 unchanged when no remap is warranted (multi-chain
    inputs, matching chain IDs, or sequence-id below threshold).
    """
    chains1 = {a['chain'] for a in atoms1}
    chains2 = {a['chain'] for a in atoms2}
    if chains1 == chains2:
        return atoms2
    if len(chains1) != 1 or len(chains2) != 1:
        return atoms2   # multi-chain — too ambiguous, leave alone

    # Compare CA residue identities at each (resnum, icode)
    r1 = {(a['resnum'], a.get('icode', '')): a['resname']
          for a in atoms1 if a.get('is_ca')}
    r2 = {(a['resnum'], a.get('icode', '')): a['resname']
          for a in atoms2 if a.get('is_ca')}
    common = set(r1) & set(r2)
    if not common:
        return atoms2
    matches = sum(1 for k in common if r1[k] == r2[k])
    frac = matches / len(common)
    if frac < seq_id_min:
        if verbose:
            print(f"  match_atoms: chains differ ('{next(iter(chains1))}' vs "
                  f"'{next(iter(chains2))}') but sequence id only "
                  f"{100*frac:.0f}% — refusing to remap", flush=True)
        return atoms2

    src    = next(iter(chains2))
    target = next(iter(chains1))
    if verbose:
        print(f"  match_atoms: chain-ID remap '{src}' → '{target}' "
              f"({matches}/{len(common)} CA pairs same resname, "
              f"{100*frac:.0f}% sequence id)", flush=True)

    # Shallow-copy each atom; rewrite chain + the leading uid token
    remapped = []
    for a in atoms2:
        if a['chain'] != src:
            remapped.append(a)
            continue
        b = dict(a)
        b['chain'] = target
        parts = a['uid'].split('_')
        parts[0] = target
        b['uid'] = '_'.join(parts)
        remapped.append(b)
    return remapped


def match_atoms(atoms1, atoms2):
    """Match atoms between two P1 lists by uid.

    Returns
    -------
    fitme   : ndarray (N, 8) — xf yf zf  dxf dyf dzf  do  dB
    ca_mask : ndarray (N,) bool
    uids    : list[str]
    bfacs   : ndarray (N, 2) — [B1, B2] for each matched pair
    """
    # Auto-remap single-chain mismatches (e.g. 'X' vs 'A') if the residue
    # sequences agree.  No-op when both inputs already use the same chains
    # or the multi-chain case is too ambiguous to remap safely.
    atoms2 = _maybe_remap_single_chain(atoms1, atoms2)

    idx2 = {a['uid']: a for a in atoms2}
    rows, ca, uids, bfacs = [], [], [], []
    for a1 in atoms1:
        a2 = idx2.get(a1['uid'])
        if a2 is None:
            continue
        dx = a2['x'] - a1['x']
        dy = a2['y'] - a1['y']
        dz = a2['z'] - a1['z']
        # Minimum image convention — wrap shifts to (-0.5, 0.5]
        dx -= round(dx);  dy -= round(dy);  dz -= round(dz)
        rows.append([
            a1['x'], a1['y'], a1['z'],
            dx, dy, dz,
            a2['occ']  - a1['occ'],
            a2['bfac'] - a1['bfac'],
        ])
        ca.append(a1['is_ca'])
        uids.append(a1['uid'])
        bfacs.append([a1['bfac'], a2['bfac']])
    if not rows:
        raise ValueError("No atoms matched between the two structures — "
                         "check space groups match.")
    return (np.array(rows, dtype=float), np.array(ca, dtype=bool),
            uids, np.array(bfacs, dtype=float))


# ══════════════════════════════════════════════════════════════════════════════
# HKL generation
# ══════════════════════════════════════════════════════════════════════════════

def generate_hkls(cell, reso):
    """Friedel-unique (h,k,l) with d ≥ reso Å, sorted by d descending.

    (0,0,0) is prepended as the DC term.
    Returns ndarray (N, 3) int.
    """
    Gs = _reciprocal_metric(cell)
    a, b, c = cell[0], cell[1], cell[2]
    hmax = int(np.ceil(a / reso)) + 1
    kmax = int(np.ceil(b / reso)) + 1
    lmax = int(np.ceil(c / reso)) + 1

    hkls_d = []
    for h in range(-hmax, hmax+1):
        for k in range(-kmax, kmax+1):
            for l in range(-lmax, lmax+1):
                if h == 0 and k == 0 and l == 0:
                    continue
                hv = np.array([h, k, l], dtype=float)
                inv_d2 = float(hv @ Gs @ hv)
                if inv_d2 <= 0:
                    continue
                d = 1.0 / np.sqrt(inv_d2)
                if d < reso - 1e-9:
                    continue
                # Friedel dedup: flip so first nonzero index is positive
                for v in (h, k, l):
                    if v != 0:
                        if v < 0:
                            h, k, l = -h, -k, -l
                        break
                hkls_d.append((h, k, l, d))

    hkls_d.sort(key=lambda x: -x[3])   # low-res first

    seen = set()
    unique = []
    for h, k, l, d in hkls_d:
        key = (h, k, l)
        if key not in seen:
            seen.add(key)
            unique.append(key)

    return np.array([(0, 0, 0)] + unique, dtype=int)


# ══════════════════════════════════════════════════════════════════════════════
# Space-group symmetry for PSDVF fitting
# ══════════════════════════════════════════════════════════════════════════════

def get_proper_symops(sg_name):
    """Return list of (R, t) for proper (det=+1) rotation operators of sg_name.

    R : (3,3) float — rotation matrix in fractional coordinates
    t : (3,)  float — translation in fractional coordinates
    """
    sg = gemmi.find_spacegroup_by_name(sg_name)
    result = []
    for op in sg.operations():
        R = np.array(op.rot, dtype=float) / 24.0
        t = np.array(op.tran, dtype=float) / 24.0
        if round(np.linalg.det(R)) == 1:
            result.append((R, t))
    return result


def assign_canonical(hkls, proper_ops):
    """Assign each Friedel-unique HKL to its canonical orbit representative.

    The canonical HKL is the lexicographic minimum of {R_k^T H} over all
    proper operators k, reduced to Friedel-unique form.

    Returns
    -------
    canon_hkls   : (M, 3) int — one canonical HKL per orbit (M ≤ N)
    hkl_to_canon : (N,) int  — index into canon_hkls for each input HKL
    hkl_all_ops  : list of N lists of (k, flip) — ALL operators mapping
                   R_k^T canon → ±HKL[i]; flip=True means negation was needed
    """
    rot_mats = [np.round(R).astype(int) for R, _ in proper_ops]
    hkl_set = {tuple(int(x) for x in h) for h in hkls}

    def to_friedel(h):
        """Return (friedel_form, flipped).  flipped=True if negation was needed."""
        h = tuple(int(round(x)) for x in h)
        if h in hkl_set:
            return h, False
        neg = tuple(-x for x in h)
        if neg in hkl_set:
            return neg, True
        return None, False

    canon_map = {}
    canon_list = []
    hkl_to_canon = np.zeros(len(hkls), dtype=int)
    hkl_all_ops  = []          # list of lists; filled in second pass below

    for i, h in enumerate(hkls):
        h_t = tuple(int(x) for x in h)
        orbit = []
        for j, R in enumerate(rot_mats):
            fu, _ = to_friedel(R.T @ np.array(h))
            if fu is not None:
                orbit.append((fu, j))
        if not orbit:
            orbit = [(h_t, 0)]

        canon = min(set(fu for fu, _ in orbit))   # lexicographic min

        if canon not in canon_map:
            canon_map[canon] = len(canon_list)
            canon_list.append(canon)
        hkl_to_canon[i] = canon_map[canon]
        hkl_all_ops.append([])    # placeholder; filled below

    canon_hkls = np.array(canon_list, dtype=int)

    # Second pass: for each HKL i find ALL (k, flip) where R_k^T canon → ±H_i
    for i, h in enumerate(hkls):
        h_t = tuple(int(x) for x in h)
        m   = hkl_to_canon[i]
        canon = canon_hkls[m]
        ops_for_i = []
        for k, R in enumerate(rot_mats):
            fu, flipped = to_friedel(R.T @ canon)
            if fu == h_t:
                ops_for_i.append((k, flipped))
        hkl_all_ops[i] = ops_for_i

    return canon_hkls, hkl_to_canon, hkl_all_ops


def build_design_matrix_symm(frac_coords, canon_hkls, proper_ops,
                             hkl_to_canon, hkl_all_ops):
    """Build (3N, 6M) symmetry-adapted design matrix for PSDVF.

    Free parameters per canonical HKL m: A_x,A_y,A_z (cols 6m+0..2) and
    B_x,B_y,B_z (cols 6m+3..5).  Row 3n+α = atom n, dimension α.

    X[3n+α, 6m+β]   = Σ_k R_k^{-1}[α,β] sin(2π H_m·(R_k x_n + t_k))
    X[3n+α, 6m+3+β] = Σ_k R_k^{-1}[α,β] cos(2π H_m·(R_k x_n + t_k))

    NOTE: the vectorial factor is R_k^{-1}[α,β], NOT R_k^T[α,β].  These are
    equal when R_k is orthogonal (cubic/tetragonal/orthorhombic/monoclinic)
    but differ for trigonal/hexagonal fractional rotation matrices.

    Only operators k whose image R_k^T H_c is in the Friedel-unique HKL basis
    are included (via hkl_all_ops).  Operators mapping to HKLs outside the
    truncated basis are omitted, keeping the design matrix consistent with
    expand_ab_canon.
    """
    N, M = len(frac_coords), len(canon_hkls)
    X = np.zeros((3 * N, 6 * M))
    canon_f = canon_hkls.astype(float)

    # valid_mk[m, k] = True if operator k maps H_c_m to a basis HKL
    n_ops = len(proper_ops)
    valid_mk = np.zeros((M, n_ops), dtype=bool)
    for j in range(len(hkl_all_ops)):
        m = hkl_to_canon[j]
        for k, _flip in hkl_all_ops[j]:
            valid_mk[m, k] = True

    for ki, (R, t) in enumerate(proper_ops):
        R_inv = np.linalg.inv(R)                        # R^{-1} for vector transform
        x_t = frac_coords @ R.T + t                    # (N, 3): R_k x + t_k
        ph  = TWO_PI * (x_t @ canon_f.T)               # (N, M): 2π H_c·(R_k x+t_k)
        s = np.sin(ph)                                  # (N, M)
        c = np.cos(ph)
        # Zero columns for canonicals where this op maps outside the basis
        s[:, ~valid_mk[:, ki]] = 0.0
        c[:, ~valid_mk[:, ki]] = 0.0
        for alpha in range(3):
            for beta in range(3):
                r = R_inv[alpha, beta]                  # R_k^{-1}[alpha, beta]
                if r == 0.0:
                    continue
                X[alpha::3, beta::6]   += r * s
                X[alpha::3, 3+beta::6] += r * c
    return X


def _pnn_weight(snr, N):
    """Holton 'probability not noise' weight: P(snr is real | N HKLs tested).

    Pnn(snr, N) = erf(|snr|/√2)^N is the probability that an SNR this large
    would NOT be exceeded by random noise across N independent measurements
    (Bonferroni-style multiple-testing correction).  As N grows the SNR
    bar for "real signal" rises automatically — at N=500, snr<3 essentially
    zeroes, snr=4 ≈ 0.97 weight, snr=5 ≈ 1.

    Below the safe threshold (where Pnn would underflow to ~0), return 0
    exactly to avoid floating-point garbage.

    Note: half-normal is a slight statistical mis-match for the way snr
    is constructed in fit_lstsq_symm (joint |AB| / √Σσᵢ², which under H0
    is χ_k/√k for k=6 not a single half-normal — see _pnn_weight_chik).
    The half-normal form is retained as the empirical default because the
    SVD basis is correlated; the effective multiple-testing N is much
    smaller than the raw HKL count, and the looser half-normal CDF
    accidentally compensates.  Replacing with the matched chi-k Pnn
    over-suppresses moderate-snr HKLs at fine fitreso (8jee fr10 bondZ
    25 with chi-6 vs 6 with half-normal/softPnna) — leaving the SVD with
    too few effective basis vectors and producing ringing instead of a
    smooth field.  Keep half-normal as the workhorse; chi-k is available
    as `_pnn_weight_chik` for future experimentation.
    """
    snr_abs = np.abs(np.asarray(snr, dtype=float))
    s_thresh = np.sqrt(2.0) * sp_erfinv((1e-30) ** (1.0 / max(N, 1)))
    safe = snr_abs >= s_thresh
    w = np.zeros_like(snr_abs)
    if safe.any():
        w[safe] = sp_erf(snr_abs[safe] / np.sqrt(2.0)) ** N
    return w


def _pnn_weight_chik(snr, N, k=6):
    """Matched-null Pnn weight: F_χ_k(s·√k)^N via scipy gammainc.

    Theoretically correct for the joint-magnitude snr used by
    fit_lstsq_symm (k = 2 · N_dims).  Empirically too aggressive at fine
    fitreso because the SVD basis is correlated — effective N < raw N.
    Not used by default; retained for diagnostic comparison."""
    snr_abs = np.abs(np.asarray(snr, dtype=float))
    if N <= 0:
        return np.ones_like(snr_abs)
    F = sp_gammainc(k / 2.0, k * snr_abs * snr_abs / 2.0)
    logF = np.log(np.clip(F, 1e-300, 1.0))
    return np.exp(np.clip(N * logF, -700.0, 0.0))


def _soft_pnn_weight(snr, N):
    """Smoother polynomial approximation of Pnn (Holton's softPnna).

    Replaces the sharp erf-power Pnn transition with a generalised-logistic
    centred on a polynomial-fitted "50% threshold sigma" that varies with
    log10(N).  The exponent of the transition adapts with both N and snr
    so the curve is broader at low N and steeper at high N.

      softPnna(d, N) = 1 - 2^(-(d / σ₅₀(N))^p(d,N))

    where σ₅₀(N) ≈ a₀ + a₁ log₁₀N + a₂ (log₁₀N)² is the polynomial fit of
    the true Pnn 50% sigma, and p(d,N) = my·log₁₀N + mx·d + y₀.

    The result is a smoother tail than Pnn: low-snr HKLs get a small but
    nonzero weight (instead of hard-zero), high-snr HKLs ramp to 1 more
    gradually — better suited to fits where some signal at moderate snr
    is still useful.
    """
    Y0 = 1.0
    A2, A1, A0 = -0.0192266, 0.751694, 1.12482
    MX, MY = 0.21805, 0.736621

    snr_abs = np.abs(np.asarray(snr, dtype=float))
    logN = max(0.0, np.log10(max(N, 1)))

    asigma_50 = A2 * logN * logN + A1 * logN + A0
    if asigma_50 <= 0:
        asigma_50 = 1.0

    expo  = MY * logN + MX * snr_abs + Y0
    ratio = snr_abs / asigma_50

    with np.errstate(over='ignore', invalid='ignore'):
        val = -np.power(np.maximum(ratio, 1e-30), expo)
    val = np.clip(val, -1000.0, 0.0)
    return 1.0 - np.power(2.0, val)


def _pnn_weight_k(snr, N, k=1):
    """Exact k-th-from-top Pnn: probability that random noise across N tests
    would NOT produce ≥k deviates exceeding `snr`.

    For k=1 this is the standard Pnn (top order statistic): F^N.  For
    k≥2 the bar lowers — we only need fewer than k noise deviates to
    exceed snr.  Closed form via the order-statistic CDF:

        P( Y_(N-k+1) ≤ snr ) = I_F(snr) (N-k+1, k)

    where F = erf(snr/√2) and I is the regularized incomplete beta —
    scipy.special.betainc.  Numerically well-conditioned for all 1 ≤ k ≤ N.

    Use in step-down ranking: sort HKLs by snr descending and weight
    the i-th-ranked HKL by `_pnn_weight_k(snr_i, N, k=i)`.  The top HKL
    faces F^N (strict); the second faces F^N + N·F^(N-1)·(1-F); etc.
    Each successive 'spent' deviate widens the admit zone by one
    Bonferroni step — analogous to Holm-Bonferroni for sequential
    hypothesis testing.

    Diagnostic helper; not used by default.
    """
    snr_abs = np.abs(np.asarray(snr, dtype=float))
    if N <= 0 or k < 1 or k > N:
        return np.zeros_like(snr_abs)
    F = sp_erf(snr_abs / np.sqrt(2.0))
    return sp_betainc(N - k + 1, k, F)


def _soft_pnn_weight_k(snr, N, k=1):
    """Soft k-th-from-top Pnn: softPnna with effective N = N − k + 1.

    The step-down equivalent of softPnna.  Each of the (k−1) prior
    deviates is 'spent' on a higher-ranked test, so the multiple-
    testing budget shrinks by one per rank:

        softPnna_k(snr, N, k)  =  softPnna(snr, N − k + 1)

    This is the conditional (Holm) form, not the unconditional kth-
    order-statistic CDF — but matches the spirit of softPnna's
    polynomial-smoothed transition and reuses the same fitted
    constants.  For k=1 reduces to plain softPnna.

    Compared to the exact `_pnn_weight_k`:
      - exact form jumps sharply at σ_50(N, k);
      - soft form keeps moderate-snr HKLs with nonzero weight, which
        empirically holds up under SVD basis correlation (see chi-k
        note in `_pnn_weight`).

    Diagnostic helper; not used by default.
    """
    N_eff = max(1, int(N) - int(k) + 1)
    return _soft_pnn_weight(snr, N_eff)


def fit_lstsq_symm(X, b, drop_snr=0.0, max_rounds=5, bound_by_obs_frac=None,
                   pnn_mode=None):
    """Joint 3D lstsq with per-canonical-HKL Pnn weighting.

    X : (3N, 6M) symmetry-adapted design matrix
    b : (3N,)    interleaved shifts [Δx_0,Δy_0,Δz_0, Δx_1,...]

    Behaviour by `drop_snr`:
      - drop_snr == 0 (default): one SVD pass, then apply Pnn weight
        w_m = erf(|snr_m|/√2)^M to each 6-parameter canonical block,
        where M is the number of canonical HKLs (built-in multiple-
        testing correction).  Snr ≫ required threshold ⇒ w ≈ 1;
        snr at/below threshold ⇒ w ≈ 0.  No HKLs are dropped; all
        returned as active.  Stored snr is the RAW post-SVD SNR
        (pre-weight) so PSDVF.mtz diagnostics show what the weighting
        was based on.
      - drop_snr > 0: legacy hard-cut.  Iteratively drop HKLs with
        snr < drop_snr (up to max_rounds) and re-solve.  Surviving
        HKLs are returned at full (un-weighted) strength.

    `bound_by_obs_frac` (optional): max observed shift magnitude in
    FRACTIONAL units.  When set, the SVD pseudo-inverse `1/s` is
    replaced with the Tikhonov-filtered form `s/(s² + λ)` with
    `λ = σ_noise² · M / obs²` — calibrated so the prior on each of
    the 6M parameters has variance obs²/M (Parseval bound, total field
    energy ≤ obs²).  Shrinks small singular values; high-SNR
    components pass through unchanged.  Bounds the fitted field
    without per-HKL weighting.

    Returns
    -------
    params       : (6M,) float — fitted (A,B) per canonical HKL,
                                 Pnn-weighted when drop_snr == 0
    active_canon : (M,)  bool  — surviving canonical HKLs
    snr_canon    : (M,)  float — raw SNR per canonical HKL (pre-weight)
    """
    M = X.shape[1] // 6
    n = X.shape[0]
    active = np.ones(M, dtype=bool)
    params_full = np.zeros(6 * M)
    snr_full    = np.zeros(M)

    for _ in range(max_rounds):
        col_mask = np.repeat(active, 6)
        Xa = X[:, col_mask]
        ai = np.where(active)[0]

        try:
            U, s, Vt = sp_svd(Xa, full_matrices=False, check_finite=False)
        except np.linalg.LinAlgError:
            U, s, Vt = sp_svd(Xa, full_matrices=False, check_finite=False,
                               lapack_driver='gesvd')
        s_thresh = s.max() * max(Xa.shape) * np.finfo(float).eps * 100
        s_inv = np.zeros_like(s)
        nz = s > s_thresh
        s_inv[nz] = 1.0 / s[nz]

        # First-pass plain-inverse residual gives σ_noise² for the Tikhonov λ
        sol_lsq = Vt.T @ (s_inv * (U.T @ b))
        if bound_by_obs_frac is not None and bound_by_obs_frac > 0:
            resid_lsq = b - Xa @ sol_lsq
            dof_lsq   = max(n - int(nz.sum()), 1)
            sig2_lsq  = float(np.dot(resid_lsq, resid_lsq)) / dof_lsq
            M_act     = int(active.sum())
            # σ_prior² per param = obs_frac² / (6·M_act) — total prior energy = obs²
            sigma_prior2 = (bound_by_obs_frac ** 2) / max(6 * M_act, 1)
            lam = sig2_lsq / sigma_prior2
            # Tikhonov SVD filter: s/(s² + λ) replacing 1/s
            s_inv = np.where(nz, s / (s * s + lam), 0.0)

        sol   = Vt.T @ (s_inv * (U.T @ b))
        resid = b - Xa @ sol
        rank  = int(nz.sum())
        dof   = max(n - rank, 1)
        sig2  = float(np.dot(resid, resid)) / dof

        cov_diag = sig2 * np.sum((Vt * s_inv[:, None])**2, axis=0)   # (p,)

        # SNR per canonical HKL: 6-parameter block norm
        p_blk   = sol.reshape(-1, 6)        # (M_active, 6)
        cov_blk = cov_diag.reshape(-1, 6)
        a     = np.sqrt(np.sum(p_blk**2,   axis=1))
        sig_a = np.sqrt(np.sum(cov_blk,    axis=1))   # sqrt(sum of variances)
        with np.errstate(invalid='ignore', divide='ignore'):
            snr_a = np.where(sig_a > 0, a / sig_a, 0.0)

        snr_full[ai] = snr_a    # raw SNR — record BEFORE any weighting

        if drop_snr <= 0:
            # Default per-HKL weight: with ridge on, the Holm step-down
            # softPnna_kth wins the empirical gamut (lipox fr5 RMSD
            # 230 Å→0.59 Å, all other systems unchanged at best).
            # Without ridge, fall back to constant-N softPnna.
            mode = pnn_mode
            if mode is None:
                mode = 'softpnna_kth' if (bound_by_obs_frac is not None and
                                           bound_by_obs_frac > 0) else 'softpnna'

            N_act = len(snr_a)
            # Rank-conditional N_eff array for step-down variants
            ranks = np.empty(N_act, dtype=int)
            order_desc = np.argsort(-np.abs(snr_a))
            ranks[order_desc] = np.arange(1, N_act + 1)   # rank 1 = top |snr|
            N_eff_per = np.maximum(1, N_act - ranks + 1)  # Holm: N_eff for rank-k

            if mode == 'off':
                params_full[col_mask] = sol
            else:
                if mode == 'chik':
                    w = _pnn_weight_chik(snr_a, N_act, k=6)
                elif mode == 'softpnna_kth':
                    # softPnna with N_eff = N - rank + 1 per HKL (step-down)
                    w = np.zeros_like(snr_a)
                    for k_val in np.unique(N_eff_per):
                        m_k = (N_eff_per == k_val)
                        w[m_k] = _soft_pnn_weight(snr_a[m_k], int(k_val))
                elif mode == 'pnn_kth':
                    # Strict erf^N with N_eff = N - rank + 1 per HKL
                    w = np.zeros_like(snr_a)
                    for k_val in np.unique(N_eff_per):
                        m_k = (N_eff_per == k_val)
                        w[m_k] = _pnn_weight(snr_a[m_k], int(k_val))
                elif mode == 'chik_kth':
                    # chi-k with N_eff = N - rank + 1 per HKL
                    w = np.zeros_like(snr_a)
                    for k_val in np.unique(N_eff_per):
                        m_k = (N_eff_per == k_val)
                        w[m_k] = _pnn_weight_chik(snr_a[m_k], int(k_val), k=6)
                else:  # 'softpnna'
                    w = _soft_pnn_weight(snr_a, N_act)
                params_full[col_mask] = (p_blk * w[:, None]).ravel()
            break

        params_full[col_mask] = sol
        drop = snr_a < drop_snr
        if not drop.any():
            break
        active[ai[drop]] = False

    params_full[np.repeat(~active, 6)] = 0.0
    return params_full, active, snr_full


def expand_ab_canon(params, active_canon, canon_hkls, hkls,
                    hkl_to_canon, hkl_all_ops, proper_ops):
    """Expand canonical (A,B) params to full Friedel-unique AB (3, N_hkls, 2).

    Each Friedel-unique HKL j collects contributions from ALL operators k
    that map canon_m to ±H_j.  For operator k with flip flag f:
        A_j += (-1)^f * R_k^T (A_c cos φ_k - B_c sin φ_k)   [sin(-H·r)=-sin(H·r)]
        B_j +=           R_k^T (A_c sin φ_k + B_c cos φ_k)
    where φ_k = 2π H_c · t_k.
    """
    N_hkls = len(hkls)
    AB = np.zeros((3, N_hkls, 2))
    for j in range(N_hkls):
        m = hkl_to_canon[j]
        if not active_canon[m]:
            continue
        Ac  = params[6*m:6*m+3]
        Bc  = params[6*m+3:6*m+6]
        Hc  = canon_hkls[m].astype(float)
        A_j = np.zeros(3)
        B_j = np.zeros(3)
        for k, flip in hkl_all_ops[j]:
            R, t = proper_ops[k]
            R_inv = np.linalg.inv(R)                # R^{-1} (= R^T only for orthogonal R)
            phi  = TWO_PI * float(np.dot(Hc, t))
            cph, sph = np.cos(phi), np.sin(phi)
            dA = R_inv @ (Ac * cph - Bc * sph)
            dB = R_inv @ (Ac * sph + Bc * cph)
            A_j += -dA if flip else dA
            B_j += dB
        AB[:, j, 0] = A_j
        AB[:, j, 1] = B_j
    return AB


# ══════════════════════════════════════════════════════════════════════════════
# Design matrix and linear least-squares fit
# ══════════════════════════════════════════════════════════════════════════════

def build_design_matrix(frac_coords, hkls):
    """Build (N_atoms, 2*N_hkls) design matrix.

    Column 2j   = sin(2π H_j · r_i)
    Column 2j+1 = cos(2π H_j · r_i)
    """
    phase = TWO_PI * (frac_coords @ hkls.astype(float).T)   # (N, M)
    X = np.empty((phase.shape[0], 2 * phase.shape[1]))
    X[:, 0::2] = np.sin(phase)
    X[:, 1::2] = np.cos(phase)
    return X


def _snr_from_svd(U, s, Vt, s_inv, Xa, b, n):
    """Solve one dimension given a pre-computed thin SVD of Xa.

    Returns (AB array (M,2), snr array (M,)).
    The SVD (U, s, Vt, s_inv) is shared across dimensions to avoid recomputing.
    """
    sol   = Vt.T @ (s_inv * (U.T @ b))
    resid = b - Xa @ sol
    rank  = int((s_inv > 0).sum())
    dof   = max(n - rank, 1)
    sig2  = float(np.dot(resid, resid)) / dof

    # cov diagonal: sig2 * diag(V S⁻² Vᵀ)  — O(p²) not O(n·p)
    cov_diag = sig2 * np.sum((Vt * s_inv[:, None])**2, axis=0)   # (p,)

    AB     = sol.reshape(-1, 2)
    cov_AB = cov_diag.reshape(-1, 2)
    A, B   = AB[:, 0], AB[:, 1]
    varA, varB = cov_AB[:, 0], cov_AB[:, 1]

    a = np.sqrt(A**2 + B**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        sig_a = np.where(a > 0,
                         np.sqrt(np.abs(A**2 * varA + B**2 * varB)) / np.maximum(a, 1e-30),
                         0.0)
    snr = np.where(sig_a > 0, a / sig_a, 0.0)
    return AB, snr


def fit_lstsq(X, shifts, drop_snr=0.0, max_rounds=5):
    """Linear fit with per-HKL Pnn weighting (or legacy hard-cut).

    SVD of the design matrix is computed ONCE per round and shared across all
    dimensions — ~3× faster than per-dimension SVD.

    Behaviour by `drop_snr`: see fit_lstsq_symm — same semantics.  Default
    (drop_snr == 0) applies the Pnn weight w = erf(|snr|/√2)^N (Holton,
    multiple-test corrected) to each HKL's AB pair.  drop_snr > 0 is the
    legacy iterative hard-cut path.

    Parameters
    ----------
    X        : (N_atoms, 2*N_hkls) design matrix
    shifts   : (N_atoms, N_dims) fractional shifts
    drop_snr : float — hard-cut SNR threshold (0 = Pnn-weight instead)

    Returns
    -------
    AB     : (N_dims, N_hkls, 2) — Pnn-weighted (default) or raw
    active : (N_hkls,) bool — surviving HKLs
    snr    : (N_hkls,) float — raw mean SNR across dimensions (pre-weight)
    """
    N_hkls = X.shape[1] // 2
    N_dims  = shifts.shape[1]
    n       = X.shape[0]
    active  = np.ones(N_hkls, dtype=bool)
    AB_full  = np.zeros((N_dims, N_hkls, 2))
    snr_full = np.zeros(N_hkls)

    for _ in range(max_rounds):
        col_mask = np.repeat(active, 2)
        Xa = X[:, col_mask]
        ai = np.where(active)[0]

        # One SVD shared across all dimensions
        try:
            U, s, Vt = sp_svd(Xa, full_matrices=False, check_finite=False)
        except np.linalg.LinAlgError:
            U, s, Vt = sp_svd(Xa, full_matrices=False, check_finite=False,
                               lapack_driver='gesvd')
        s_thresh = s.max() * max(Xa.shape) * np.finfo(float).eps * 100
        s_inv = np.zeros_like(s)
        nz = s > s_thresh
        s_inv[nz] = 1.0 / s[nz]

        snr_sum = np.zeros(len(ai))
        for d in range(N_dims):
            AB_d, snr_d = _snr_from_svd(U, s, Vt, s_inv, Xa, shifts[:, d], n)
            AB_full[d, ai, :] = AB_d
            snr_sum += snr_d

        snr_active = snr_sum / N_dims
        snr_full[ai] = snr_active

        if drop_snr <= 0:
            # Pnn weight with N = number of HKLs being tested (multiple-
            # testing correction; see _pnn_weight)
            w = _soft_pnn_weight(snr_active, len(snr_active))      # (N_active,)
            AB_full[:, ai, :] *= w[None, :, None]
            break
        drop = snr_active < drop_snr
        if not drop.any():
            break
        active[ai[drop]] = False

    # Zero out dropped HKLs so they don't contribute to the shift field
    AB_full[:, ~active, :] = 0.0
    return AB_full, active, snr_full


# ══════════════════════════════════════════════════════════════════════════════
# Shift field evaluation and RMSD
# ══════════════════════════════════════════════════════════════════════════════

def eval_shift_field(frac_coords, hkls, AB):
    """Evaluate the shift field at fractional positions.

    AB : (N_dims, N_hkls, 2)
    Returns (N_pos, N_dims) fractional shifts.
    """
    phase  = TWO_PI * (frac_coords @ hkls.astype(float).T)   # (N, M)
    sin_ph = np.sin(phase)
    cos_ph = np.cos(phase)
    # (N_dims, M, 1) × (M, N) → (N_dims, N) → (N, N_dims)
    shifts = AB[:, :, 0] @ sin_ph.T + AB[:, :, 1] @ cos_ph.T
    return shifts.T


def _shift_magnitude_fourier(hkls, AB_xyz, cell):
    """Fourier amplitudes/phases of the displacement-magnitude map |Δr(x)|.

    The shift field has three independent vector components per voxel;
    |Δr|(x) = √(δx_orth² + δy_orth² + δz_orth²) is a derived real-space
    scalar map.  It is NOT a linear combination of the AB Fourier
    components — it's a non-linear function of them — so to express it
    as (amplitude, phase) per HKL we sample the field on a real-space
    grid, take voxelwise |Δr|, FFT, and look up at the input HKLs.

    Diagnostic in Coot: large |Δr| outside the protein body signals the
    field ringing where there are no constraints (a tell-tale of over-fit
    or wrong-direction bending).

    Parameters
    ----------
    hkls   : (N, 3) int   — same HKL set written to the rest of PSDVF.mtz
    AB_xyz : (3, N, 2)    — fractional-axis A, B coefficients
    cell   : 6-tuple      — unit cell (a,b,c,α,β,γ)

    Returns
    -------
    dR_amp   : (N,) float32 — |F(h)| of the |Δr| map (Å³, gemmi convention)
    PHIR_deg : (N,) float32 — phase in degrees
    """
    if len(hkls) == 0 or AB_xyz.shape[0] != 3:
        return (np.zeros(len(hkls), dtype=np.float32),
                np.zeros(len(hkls), dtype=np.float32))

    # Grid sized to comfortably resolve the highest miller index in the set
    hmax = max(int(np.abs(hkls[:, 0]).max()), 1)
    kmax = max(int(np.abs(hkls[:, 1]).max()), 1)
    lmax = max(int(np.abs(hkls[:, 2]).max()), 1)
    nx, ny, nz = 4 * (hmax + 1), 4 * (kmax + 1), 4 * (lmax + 1)

    # Fractional coordinate grid (NX × NY × NZ)
    xs = np.arange(nx, dtype=np.float64) / nx
    ys = np.arange(ny, dtype=np.float64) / ny
    zs = np.arange(nz, dtype=np.float64) / nz
    Xg, Yg, Zg = np.meshgrid(xs, ys, zs, indexing='ij')
    frac_pts = np.stack([Xg.ravel(), Yg.ravel(), Zg.ravel()], axis=1)

    # Field at each voxel → orthogonal Å → scalar magnitude
    delta_frac = eval_shift_field(frac_pts, hkls, AB_xyz)        # (N_voxels, 3)
    delta_orth = frac_to_orth(delta_frac, cell)                  # (N_voxels, 3)
    mag = np.linalg.norm(delta_orth, axis=1).reshape((nx, ny, nz))

    # FFT and look up at each HKL using the +2πi crystallographic
    # convention (gemmi/refmac), matching the existing _map2mtz path.
    F = np.fft.fftn(mag)
    H = hkls[:, 0].astype(int)
    K = hkls[:, 1].astype(int)
    L = hkls[:, 2].astype(int)
    cell_vol = gemmi.UnitCell(*cell).volume
    F_h = F[(-H) % nx, (-K) % ny, (-L) % nz] * (cell_vol / (nx * ny * nz))

    dR_amp = np.abs(F_h).astype(np.float32)
    PHIR_deg = np.degrees(np.angle(F_h)).astype(np.float32)
    return dR_amp, PHIR_deg


def eval_divergence(frac_coords, hkls, AB):
    """Analytic divergence div(δ) = ∂δx/∂x + ∂δy/∂y + ∂δz/∂z in fractional coords.

    Equal to det(J) - 1 to first order, where J = I + ∇δ is the Jacobian of the
    forward transformation.  Multiply resampled density by (1 - div(δ)) to
    account for the change in voxel volume under the transformation.

    AB : (3, N_hkls, 2) — must have exactly 3 dimensions (x, y, z).
    Returns (N_pos,) divergence values.
    """
    phase  = TWO_PI * (frac_coords @ hkls.astype(float).T)   # (N, M)
    sin_ph = np.sin(phase)
    cos_ph = np.cos(phase)
    h, k, l = hkls[:, 0].astype(float), hkls[:, 1].astype(float), hkls[:, 2].astype(float)
    # ∂δx/∂x = Σ_j 2π h_j (Ax_j cos - Bx_j sin)
    ddx = TWO_PI * ((AB[0, :, 0] * cos_ph - AB[0, :, 1] * sin_ph) @ h)
    ddy = TWO_PI * ((AB[1, :, 0] * cos_ph - AB[1, :, 1] * sin_ph) @ k)
    ddz = TWO_PI * ((AB[2, :, 0] * cos_ph - AB[2, :, 1] * sin_ph) @ l)
    return ddx + ddy + ddz


def rmsd_ca(fitme, ca_mask, hkls, AB_xyz, cell2, frac=1.0):
    """RMSD of CA atoms in Å between predicted bent pdb1 and pdb2.

    AB_xyz : (3, N_hkls, 2) — x, y, z dimensions.
    cell2  : tuple — reference cell for fractional → Å conversion.
    """
    fc   = fitme[:, :3]          # fractional coords (N, 3)
    true = fitme[:, 3:6]         # true fractional shifts
    pred = eval_shift_field(fc, hkls, AB_xyz)   # (N, 3)
    resid_frac = true - frac * pred
    resid_orth = frac_to_orth(resid_frac[ca_mask], cell2)
    return float(np.sqrt(np.mean(np.sum(resid_orth**2, axis=1))))


# ══════════════════════════════════════════════════════════════════════════════
# Output PDB
# ══════════════════════════════════════════════════════════════════════════════

def write_bent_pdb(pdb1_path, pdb2_path, hkls, AB_xyz, outpath, frac=1.0,
                   origin_shift=None, R_alt=None):
    """Apply shift field to pdb1, convert to pdb2 orthogonal frame, write PDB.

    AB_xyz : (3, N_hkls, 2) — x, y, z shift coefficients.
    origin_shift : (3,) fractional shift added to pdb1 frac coords before
                   evaluating the shift field (matches expand_to_p1 behaviour).
    R_alt        : (3,3) altindex rotation that was applied to pdb2 in fractional
                   coords during alignment.  After origin+delta brings atoms1 into
                   atoms2_processed (= R_alt @ atoms2_orig), apply R_alt^-1 so the
                   output bent atoms sit in pdb2's original frame.
    """
    st1 = gemmi.read_pdb(str(pdb1_path))
    st2 = gemmi.read_pdb(str(pdb2_path))

    cell1 = _cell_tuple(st1.cell)
    cell2 = _cell_tuple(st2.cell)
    M1_inv = np.linalg.inv(_ortho_matrix(cell1))
    M2     = _ortho_matrix(cell2)

    orth_list = []
    atom_refs = []
    for model in st1:
        for chain in model:
            for res in chain:
                for atom in res:
                    p = atom.pos
                    orth_list.append([p.x, p.y, p.z])
                    atom_refs.append(atom)

    orth_arr  = np.array(orth_list)                          # (N, 3)
    frac_arr  = orth_arr @ M1_inv.T                          # (N, 3)
    if origin_shift is not None:
        frac_arr = frac_arr + np.asarray(origin_shift, dtype=float)
    d_frac    = eval_shift_field(frac_arr, hkls, AB_xyz)     # (N, 3)
    bent_frac = frac_arr + frac * d_frac                     # (N, 3)  in atoms1 frame
    if R_alt is not None:
        # Rotate from atoms1's (R_alt-rotated) frame back to ref's original frame
        bent_frac = bent_frac @ np.linalg.inv(np.asarray(R_alt, dtype=float)).T
    bent_orth = bent_frac @ M2.T                             # (N, 3)

    for atom, bpos in zip(atom_refs, bent_orth):
        atom.pos = gemmi.Position(float(bpos[0]), float(bpos[1]), float(bpos[2]))

    st1.cell = st2.cell
    st1.spacegroup_hm = st2.spacegroup_hm
    st1.write_pdb(str(outpath))


_MONLIB_CACHE = {}

def _load_monlib(resnames):
    """Cached MonLib loader.  Reading the CCP4 monomer lib takes ~0.1 s for a
    typical protein; per-scan-point geometry checks would otherwise re-read it.
    Cache key is the frozenset of residue names that must be covered."""
    key = frozenset(resnames)
    if key in _MONLIB_CACHE:
        return _MONLIB_CACHE[key]
    monlib_dir = os.environ.get('CLIBD_MON')
    if not monlib_dir or not os.path.isdir(monlib_dir):
        return None
    ml = gemmi.MonLib()
    try:
        ml.read_monomer_lib(monlib_dir, sorted(resnames))
    except Exception:
        return None
    _MONLIB_CACHE[key] = ml
    return ml


def check_geometry(pdb_path, outlier_z=5.0):
    """Bond + angle RMSZ via gemmi.prepare_topology against the CCP4 monomer lib.

    The PSDVF is smooth and bandlimited at the chosen fitreso.  At fine fitreso
    the field gradient outpaces the inter-atom spacing and bonds get stretched
    or compressed by amounts the underlying smooth-deformation model is happy
    with but real chemistry is not.  This check reports the z-score deviation
    from monomer-library ideals — refmac's standard geometry metric.

    Returns dict with bond/angle RMSZ + max z + outlier count + total counts.
    Returns None if monlib loading or topology prep fails.

    Empirically: bond RMSZ ~0.5–1.0 is refined-quality; ~2–4 is a loose
    deposit; >5 is starting to break; ~100 is catastrophic.
    """
    try:
        st = gemmi.read_structure(str(pdb_path))
        st.setup_entities()
    except Exception:
        return None
    # gemmi.prepare_topology raises RuntimeError on the FIRST atom whose name
    # doesn't match its monomer-lib definition (e.g. Y/PYR 225/C1 vs C in
    # pyruvoyl-dependent aspartate decarboxylase).  One mismatched atom kills
    # topology prep for the entire structure, so bondZ silently drops to None
    # for every scan point.  Strip the offending resname and retry — this
    # keeps the protein body's geometry check while sacrificing the cofactor.
    # Iterate: some structures have multiple problem cofactors.
    dropped = []
    for _ in range(8):
        resnames = {res.name for model in st for chain in model for res in chain}
        ml = _load_monlib(resnames)
        if ml is None:
            return None
        try:
            topo = gemmi.prepare_topology(st, ml, model_index=0)
            break
        except RuntimeError as e:
            msg = str(e)
            m = re.search(r'\b([A-Za-z0-9]{1,4})\s+\d+\s*/\s*[A-Za-z0-9]', msg)
            bad = m.group(1) if m else None
            if not bad or bad in dropped:
                return None
            dropped.append(bad)
            for model in st:
                for chain in model:
                    for i in range(len(chain) - 1, -1, -1):
                        if chain[i].name == bad:
                            del chain[i]
    else:
        return None
    bonds  = topo.bonds
    angles = topo.angles
    if not bonds:
        return None
    # Some bonds/angles can lack an ideal target (e.g. linkage to a residue
    # whose monomer-lib entry doesn't fully constrain it).  Their
    # calculate_z() returns NaN and would propagate to the RMSZ.  Filter
    # out non-finite z's before reducing.
    bz_all = np.array([b.calculate_z() for b in bonds])
    az_all = np.array([a.calculate_z() for a in angles]) if angles else np.zeros(0)
    bz = bz_all[np.isfinite(bz_all)]
    az = az_all[np.isfinite(az_all)]
    if len(bz) == 0:
        return None     # truly no usable bonds left
    return dict(n_bonds=int(len(bz)),
                bond_rmsz=float(np.sqrt(np.mean(bz*bz))),
                bond_max_z=float(np.max(np.abs(bz))),
                n_bond_outliers=int(np.sum(np.abs(bz) > outlier_z)),
                n_angles=int(len(az)),
                angle_rmsz=float(np.sqrt(np.mean(az*az))) if len(az) else 0.0,
                angle_max_z=float(np.max(np.abs(az))) if len(az) else 0.0,
                n_angle_outliers=int(np.sum(np.abs(az) > outlier_z)) if len(az) else 0,
                outlier_z=outlier_z)


# ══════════════════════════════════════════════════════════════════════════════
# CCP4 map I/O (pure binary — no external dependencies)
# ══════════════════════════════════════════════════════════════════════════════

def read_ccp4(path):
    """Read a CCP4/MRC map.  Returns (data float32 (ns,nr,nc), hdr dict)."""
    with open(path, 'rb') as f:
        raw = f.read()

    i4 = lambda off: struct.unpack_from('<i', raw, off)[0]
    f4 = lambda off: struct.unpack_from('<f', raw, off)[0]

    nc, nr, ns = i4(0), i4(4), i4(8)
    mode       = i4(12)
    ncstart    = i4(16);  nrstart = i4(20);  nsstart = i4(24)
    nx, ny, nz = i4(28),  i4(32),  i4(36)
    cell_a     = (f4(40), f4(44), f4(48), f4(52), f4(56), f4(60))
    mapc, mapr, maps = i4(64), i4(68), i4(72)
    ispg       = i4(88)
    nsymbt     = i4(92)

    hdr_size = 1024 + nsymbt
    dtypes   = {0: np.int8, 1: np.int16, 2: np.float32}
    if mode not in dtypes:
        raise ValueError(f"Unsupported CCP4 map mode {mode}")
    data = np.frombuffer(raw[hdr_size:], dtype=dtypes[mode]).copy()
    data = data.reshape(ns, nr, nc).astype(np.float32)

    hdr = {
        'nc': nc, 'nr': nr, 'ns': ns, 'mode': mode,
        'ncstart': ncstart, 'nrstart': nrstart, 'nsstart': nsstart,
        'nx': nx, 'ny': ny, 'nz': nz,
        'cell': cell_a,
        'mapc': mapc, 'mapr': mapr, 'maps': maps,
        'spacegroup': ispg,
        'nsymbt': nsymbt,
        'raw_header': raw[:hdr_size],
    }
    return data, hdr


def read_ccp4_fullcell(path):
    """Read a CCP4 map and expand to the full unit cell using gemmi symmetry.

    Uses gemmi's setup() to fill the entire unit cell from the asymmetric unit,
    then rearranges axes to match the original map's mapc/mapr/maps convention.
    Returns (data float32 (ns_full, nr_full, nc_full), hdr dict with starts=0).
    """
    import gemmi
    _, hdr = read_ccp4(path)
    ccp4 = gemmi.read_ccp4_map(path, setup=True)
    grid = ccp4.grid
    # gemmi grid.array shape is (nu, nv, nw) = (NX, NY, NZ) after setup
    garr = np.array(grid.array, dtype=np.float32)
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']
    # Rearrange to (ns, nr, nc) = (maps-axis, mapr-axis, mapc-axis)
    data = np.transpose(garr, (maps - 1, mapr - 1, mapc - 1))
    Ng = {1: grid.nu, 2: grid.nv, 3: grid.nw}
    hdr_full = dict(hdr)
    hdr_full.update(nc=Ng[mapc], nr=Ng[mapr], ns=Ng[maps],
                    ncstart=0, nrstart=0, nsstart=0)
    return data, hdr_full


def _symop_block_for_sg(sg_number):
    """Return CCP4 symop string block (80-char ASCII records, no terminator).

    Each record is "X1, X2, X3" in CCP4 convention, padded to 80 chars.
    Falls back to identity (P1) if the spacegroup can't be resolved.
    """
    import gemmi as _gm
    sg = _gm.find_spacegroup_by_number(int(sg_number)) if sg_number else None
    if sg is None:
        ops = [_gm.Op('x,y,z')]
    else:
        ops = list(sg.operations())
    records = b''
    for op in ops:
        txt = op.triplet().upper().replace(' ', '').replace(',', ', ')
        records += txt.ljust(80).encode('ascii', 'replace')[:80]
    return records


def write_ccp4(path, data, hdr_template, cell_override=None):
    """Write a CCP4 map, borrowing the header from hdr_template.

    Overwrites DMIN/DMAX/DMEAN (words 20-22) and RMS (word 54) with
    actual stats.  ISPG (word 23, byte 88) is preserved from the template.
    Symop strings (after the 1024-byte header) are written from the
    spacegroup when synthesizing a fresh header — older CCP4 tools
    (mapdump, sfall) require these to be present, falling back to P1
    and corrupting downstream FFT/scaling when they're missing.
    cell_override: optional (a,b,c,al,be,ga) to replace cell in header.

    If hdr_template lacks 'raw_header' (e.g. MTZ-derived header), synthesises
    a 1024-byte CCP4 header + symop block from the dict fields.
    """
    sg_num = hdr_template.get('spacegroup')
    if 'raw_header' in hdr_template:
        raw = bytearray(hdr_template['raw_header'])
        # If the template header lacks symops (NSYMBT=0) but does name
        # a non-trivial spacegroup, regenerate them so the output map
        # is mapdump-readable.
        nsymbt_existing = struct.unpack_from('<i', raw, 92)[0]
        ispg_existing   = struct.unpack_from('<i', raw, 88)[0]
        if nsymbt_existing == 0 and ispg_existing > 1:
            symops = _symop_block_for_sg(ispg_existing)
            if symops:
                raw = bytearray(raw[:1024]) + bytearray(symops)
                struct.pack_into('<i', raw, 92, len(symops))
        symops_block = b''  # already inlined into raw above
    else:
        raw = bytearray(1024)
        struct.pack_into('<i', raw, 12, 2)       # MODE = float32
        nx = hdr_template.get('nx', hdr_template['nc'])
        ny = hdr_template.get('ny', hdr_template['nr'])
        nz = hdr_template.get('nz', hdr_template['ns'])
        for off, v in ((28, nx), (32, ny), (36, nz)):
            struct.pack_into('<i', raw, off, int(v))
        cell = hdr_template['cell']
        for i, v in enumerate(cell):
            struct.pack_into('<f', raw, 40 + i*4, float(v))
        struct.pack_into('<i', raw, 64, hdr_template.get('mapc', 1))
        struct.pack_into('<i', raw, 68, hdr_template.get('mapr', 2))
        struct.pack_into('<i', raw, 72, hdr_template.get('maps', 3))
        struct.pack_into('<i', raw, 88, int(sg_num or 1))
        raw[208:212] = b'MAP '
        raw[212:216] = b'\x44\x44\x00\x00'   # little-endian machine stamp
        symops_block = _symop_block_for_sg(sg_num or 1)
        struct.pack_into('<i', raw, 92, len(symops_block))   # NSYMBT

    mn, mx = float(data.min()), float(data.max())
    mean   = float(data.mean())
    rms    = float(np.sqrt(np.mean(data**2)))
    for off, val in zip((76, 80, 84, 216), (mn, mx, mean, rms)):
        struct.pack_into('<f', raw, off, val)

    # Update grid dimensions and starts from the header dict (may differ from template)
    for off, key in ((0,'nc'),(4,'nr'),(8,'ns'),(16,'ncstart'),(20,'nrstart'),(24,'nsstart')):
        struct.pack_into('<i', raw, off, int(hdr_template[key]))

    if cell_override is not None:
        for i, val in enumerate(cell_override):
            struct.pack_into('<f', raw, 40 + i*4, val)

    with open(path, 'wb') as f:
        f.write(raw)
        if symops_block:
            f.write(symops_block)
        f.write(data.astype(np.float32).tobytes())


def _frac_to_grid_indices(frac_xyz, hdr):
    """Convert fractional coords (N,3) → grid indices [sec, row, col] (3, N).

    ncstart/nrstart/nsstart belong to the MAPC/MAPR/MAPS axes respectively,
    not necessarily to the X/Y/Z crystallographic axes — so we build a
    per-axis-number lookup before indexing.
    """
    nx, ny, nz = hdr['nx'], hdr['ny'], hdr['nz']
    ncstart, nrstart, nsstart = hdr['ncstart'], hdr['nrstart'], hdr['nsstart']
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']

    # start[ax] = origin offset for crystallographic axis ax (1=X, 2=Y, 3=Z)
    start = {mapc: ncstart, mapr: nrstart, maps: nsstart}
    Nxyz  = {1: nx, 2: ny, 3: nz}

    gxyz = {ax: frac_xyz[:, ax-1] * Nxyz[ax] - start[ax] for ax in (1, 2, 3)}
    return np.array([gxyz[maps], gxyz[mapr], gxyz[mapc]])


def interpolate_map(data, hdr, probe_frac, pad_mode='reflect'):
    """Tricubic interpolation of map data at fractional positions.

    probe_frac : (N, 3) fractional coords to sample.
    pad_mode   : 'reflect' (default, for ASU maps — boundaries are not periodic)
                 'wrap'    (for full-cell maps — x=0 and x=1 are the same point)
    Returns (N,) float32 density values.
    Pads the array with 5 voxels before computing b-spline coefficients.
    The IIR prefilter boundary influence decays as 0.268^k, so 5 voxels
    gives < 0.15% contamination at any interior query point.
    """
    idx = _frac_to_grid_indices(probe_frac, hdr)   # (3, N)
    pad = 5
    padded = np.pad(data, pad, mode=pad_mode)
    return map_coordinates(padded, idx + pad, order=3, mode='nearest').astype(np.float32)


# ══════════════════════════════════════════════════════════════════════════════
# Bent map (replace mapman + floatgen pipeline)
# ══════════════════════════════════════════════════════════════════════════════

def _make_bent_map_data(data, hdr, hkls, AB_xyz, outpath, cell2=None,
                        use_jacobian=False):
    """Resample in-memory map (data, hdr) through the shift field; write outpath."""
    nc, nr, ns = hdr['nc'], hdr['nr'], hdr['ns']
    nx, ny, nz = hdr['nx'], hdr['ny'], hdr['nz']
    ncstart, nrstart, nsstart = hdr['ncstart'], hdr['nrstart'], hdr['nsstart']
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']

    sec_idx = np.arange(ns)
    row_idx = np.arange(nr)
    col_idx = np.arange(nc)
    sec_g, row_g, col_g = np.meshgrid(sec_idx, row_idx, col_idx, indexing='ij')

    axis_to_frac = {mapc: col_g, mapr: row_g, maps: sec_g}
    start  = {mapc: ncstart, mapr: nrstart, maps: nsstart}
    Nxyz   = {1: nx, 2: ny, 3: nz}

    frac_x = (axis_to_frac[1].ravel() + start[1]) / Nxyz[1]
    frac_y = (axis_to_frac[2].ravel() + start[2]) / Nxyz[2]
    frac_z = (axis_to_frac[3].ravel() + start[3]) / Nxyz[3]
    frac_pts = np.stack([frac_x, frac_y, frac_z], axis=1)

    delta    = eval_shift_field(frac_pts, hkls, AB_xyz)
    new_vals = interpolate_map(data, hdr, frac_pts - delta)

    if use_jacobian:
        div = eval_divergence(frac_pts, hkls, AB_xyz)
        new_vals = new_vals * (1.0 - div)

    new_data = new_vals.reshape(ns, nr, nc)
    write_ccp4(outpath, new_data, hdr, cell_override=cell2)
    print(f"bent map written to {outpath}")
    return outpath


def make_bent_map(map_path, hkls, AB_xyz, outpath, cell2=None, use_jacobian=False):
    """Resample map_path through the shift field and write outpath.

    For each grid point g in the OUTPUT map, we find where it came from
    in the INPUT map: sample at (g - δ(g)).  The negative sign is because
    the shift field maps the moving frame → reference frame; to resample
    we need the inverse.

    cell2: reference-crystal unit cell (a,b,c,α,β,γ).  The output map
    is in the reference frame, so its header should carry cell2.
    use_jacobian: if True, multiply each voxel by (1 - div(δ)) to account
    for the change in voxel volume under the coordinate transformation.
    """
    data, hdr = read_ccp4(map_path)
    return _make_bent_map_data(data, hdr, hkls, AB_xyz, outpath,
                               cell2=cell2, use_jacobian=use_jacobian)


def make_delta_maps(map_path, hkls, AB_xyz, nhkls):
    """Write delta-x, delta-y, delta-z, delta-r maps over the full P1 unit cell.

    Maps are always written as full-cell P1 (starts=0, dimensions=Nxyz) so that
    map viewers never need 'extend xtal' — the shift field is smooth everywhere
    and vector components are not scalars under space-group symmetry.
    """
    _, hdr = read_ccp4(map_path)
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']
    Nxyz = {1: hdr['nx'], 2: hdr['ny'], 3: hdr['nz']}

    # Full unit cell grid: indices 0..Nxyz[axis]-1, frac = idx/Nxyz[axis]
    nc_full = Nxyz[mapc]
    nr_full = Nxyz[mapr]
    ns_full = Nxyz[maps]
    sec_g, row_g, col_g = np.meshgrid(np.arange(ns_full), np.arange(nr_full),
                                       np.arange(nc_full), indexing='ij')
    axis_to_grid = {mapc: col_g, mapr: row_g, maps: sec_g}
    frac_pts = np.stack([axis_to_grid[1].ravel() / Nxyz[1],
                         axis_to_grid[2].ravel() / Nxyz[2],
                         axis_to_grid[3].ravel() / Nxyz[3]], axis=1)

    delta = eval_shift_field(frac_pts, hkls, AB_xyz)  # (N, 3) fractional
    delta_orth = frac_to_orth(delta, hdr['cell'])      # (N, 3) Å

    hdr_p1 = dict(hdr)
    hdr_p1.update(nc=nc_full, nr=nr_full, ns=ns_full,
                  ncstart=0, nrstart=0, nsstart=0)

    for i, label in enumerate('xyz'):
        d = delta_orth[:, i].reshape(ns_full, nr_full, nc_full).astype(np.float32)
        write_ccp4(f"delta_{label}{nhkls}.map", d, hdr_p1)
    dr = np.sqrt(np.sum(delta_orth**2, axis=1)).reshape(ns_full, nr_full, nc_full).astype(np.float32)
    write_ccp4(f"delta_r{nhkls}.map", dr, hdr_p1)
    print(f"delta maps written: delta_x/y/z/r{nhkls}.map")


def _make_delta_maps_hdr(hdr, hkls, AB_xyz, nhkls):
    """Write delta maps using a header dict (no density data needed)."""
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']
    Nxyz = {1: hdr['nx'], 2: hdr['ny'], 3: hdr['nz']}
    nc_full = Nxyz[mapc]
    nr_full = Nxyz[mapr]
    ns_full = Nxyz[maps]
    sec_g, row_g, col_g = np.meshgrid(np.arange(ns_full), np.arange(nr_full),
                                       np.arange(nc_full), indexing='ij')
    axis_to_grid = {mapc: col_g, mapr: row_g, maps: sec_g}
    frac_pts = np.stack([axis_to_grid[1].ravel() / Nxyz[1],
                         axis_to_grid[2].ravel() / Nxyz[2],
                         axis_to_grid[3].ravel() / Nxyz[3]], axis=1)
    delta      = eval_shift_field(frac_pts, hkls, AB_xyz)
    delta_orth = frac_to_orth(delta, hdr['cell'])
    hdr_p1 = dict(hdr)
    hdr_p1.update(nc=nc_full, nr=nr_full, ns=ns_full,
                  ncstart=0, nrstart=0, nsstart=0)
    for i, label in enumerate('xyz'):
        d = delta_orth[:, i].reshape(ns_full, nr_full, nc_full).astype(np.float32)
        write_ccp4(f"delta_{label}{nhkls}.map", d, hdr_p1)
    dr = np.sqrt(np.sum(delta_orth**2, axis=1)).reshape(
        ns_full, nr_full, nc_full).astype(np.float32)
    write_ccp4(f"delta_r{nhkls}.map", dr, hdr_p1)
    print(f"delta maps written: delta_x/y/z/r{nhkls}.map")


# ══════════════════════════════════════════════════════════════════════════════
# fitparams save / load (MTZ format)
# ══════════════════════════════════════════════════════════════════════════════

def save_fitparams(path, hkls, AB, active, snr, cell1, cell2, dimensions, rmsd):
    path = str(path)
    for suffix in ('.npz', '.mtz'):
        if path.endswith(suffix):
            path = path[:-len(suffix)]
            break
    path = path + '.mtz'

    n_hkls = len(hkls)
    n_dims = len(dimensions)

    mtz = gemmi.Mtz()
    mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz.cell = gemmi.UnitCell(*cell1)
    mtz.title = 'bendfinder shift-field coefficients'
    mtz.history = [
        'BENDFINDER cell1 ' + ' '.join(f'{x:.6g}' for x in cell1),
        'BENDFINDER cell2 ' + ' '.join(f'{x:.6g}' for x in cell2),
        f'BENDFINDER dimensions {dimensions}',
        f'BENDFINDER rmsd {rmsd:.6f}',
    ]

    # Whether to compute the derived |Δr| Fourier amplitudes (dR/PHIR).
    # Only meaningful when all three axes are fitted — otherwise skip.
    write_dR = (dimensions == 'xyz' or set(dimensions) == set('xyz'))

    mtz.add_dataset('bendfinder')
    mtz.add_column('H', 'H')
    mtz.add_column('K', 'H')
    mtz.add_column('L', 'H')
    for dim in dimensions:
        mtz.add_column(f'd{dim.upper()}', 'F')   # amplitude of d{x,y,z} shift (Å)
        mtz.add_column(f'PH{dim.upper()}', 'P')  # phase (degrees)
    if write_dR:
        # Derived |Δr(x)| = √(δx² + δy² + δz²) map in Å, FFT'd to (amp, phase)
        # per HKL.  Diagnostic in Coot — large |Δr| outside the protein body
        # signals over-fit / ringing.  Not a linear function of AB so this
        # is computed via real-space sample + FFT, not by combining AB.
        mtz.add_column('dR',   'F')
        mtz.add_column('PHIR', 'P')
    mtz.add_column('SNR', 'R')
    mtz.add_column('ACTIVE', 'R')

    # Convert (A, B) fractional shift coefficients to (amplitude, phase).
    # Shift field: d_frac = Σ a*cos(2π(hx+ky+lz) + phi)
    # where A = -a*sin(phi_rad), B = a*cos(phi_rad)
    # so a = sqrt(A²+B²), phi_rad = atan2(-A, B).
    # Then convert fractional amplitudes to Å by multiplying by the cell
    # length along that dimension's axis (a/b/c for x/y/z) — much more
    # interpretable for diagnostics.  load_fitparams reads the
    # `amplitude_units` history tag and back-converts to fractional for
    # all downstream consumers (eval_shift_field, write_bent_pdb, etc.)
    # so this is purely a display-units change.
    _AXIS_LEN = {'x': cell1[0], 'y': cell1[1], 'z': cell1[2]}
    mtz.history = list(mtz.history) + ['BENDFINDER amplitude_units angstrom']
    n_extra = 2 if write_dR else 0
    n_cols = 3 + 2 * n_dims + n_extra + 2
    data = np.zeros((n_hkls, n_cols), dtype=np.float32)
    data[:, 0] = hkls[:, 0]
    data[:, 1] = hkls[:, 1]
    data[:, 2] = hkls[:, 2]
    for d, dim in enumerate(dimensions):
        A, B = AB[d, :, 0], AB[d, :, 1]
        amp_frac = np.sqrt(A**2 + B**2)
        data[:, 3 + 2*d]     = amp_frac * _AXIS_LEN[dim]
        data[:, 3 + 2*d + 1] = np.degrees(np.arctan2(-A, B))
    if write_dR:
        # Need AB ordered as (x,y,z); reorder if `dimensions` was permuted.
        order = [dimensions.index(d) for d in 'xyz']
        AB_xyz = AB[order, :, :]
        dR_amp, PHIR_deg = _shift_magnitude_fourier(hkls, AB_xyz, cell1)
        data[:, 3 + 2*n_dims]     = dR_amp
        data[:, 3 + 2*n_dims + 1] = PHIR_deg
    data[:, 3 + 2*n_dims + n_extra]     = snr
    data[:, 3 + 2*n_dims + n_extra + 1] = active.astype(np.float32)

    mtz.set_data(data)
    mtz.write_to_file(path)
    return path


def load_fitparams(path):
    path = str(path)
    if not (path.endswith('.mtz') or path.endswith('.npz')):
        if os.path.exists(path + '.mtz'):
            path = path + '.mtz'
        elif os.path.exists(path + '.npz'):
            path = path + '.npz'

    # Backward compat: read old numpy format
    if path.endswith('.npz'):
        d = np.load(path, allow_pickle=False)
        return (d['hkls'], d['AB'], d['active'].astype(bool), d['snr'],
                tuple(d['cell1']), tuple(d['cell2']),
                ''.join(d['dimensions']), float(d['rmsd']))

    mtz = gemmi.read_mtz_file(path)

    cell1 = tuple(mtz.cell.parameters)
    cell2 = cell1
    dimensions = 'xyz'
    rmsd = 0.0
    amp_units = 'fractional'   # default for files written before the unit tag
    for line in mtz.history:
        parts = line.split()
        if len(parts) < 2 or parts[0] != 'BENDFINDER':
            continue
        key = parts[1]
        if key == 'cell1':
            cell1 = tuple(float(x) for x in parts[2:8])
        elif key == 'cell2':
            cell2 = tuple(float(x) for x in parts[2:8])
        elif key == 'dimensions':
            dimensions = parts[2]
        elif key == 'rmsd':
            rmsd = float(parts[2])
        elif key == 'amplitude_units' and len(parts) >= 3:
            amp_units = parts[2].lower()

    data = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]

    H = data[:, labels.index('H')].astype(int)
    K = data[:, labels.index('K')].astype(int)
    L = data[:, labels.index('L')].astype(int)
    hkls = np.stack([H, K, L], axis=1)

    n_hkls = len(hkls)
    n_dims = len(dimensions)
    AB = np.zeros((n_dims, n_hkls, 2))
    # Back-convert Å amplitudes to fractional for downstream consumers
    # (eval_shift_field expects fractional AB).  Files without the
    # `amplitude_units` history tag are read as fractional unchanged.
    _AXIS_LEN = {'x': cell1[0], 'y': cell1[1], 'z': cell1[2]}
    amp_scale = (1.0 / np.array([_AXIS_LEN[d] for d in dimensions])
                 if amp_units == 'angstrom' else np.ones(n_dims))
    for d, dim in enumerate(dimensions):
        DIM = dim.upper()
        amp_col = (f'd{DIM}' if f'd{DIM}' in labels else
                   f'D{DIM}' if f'D{DIM}' in labels else
                   f'F{DIM}' if f'F{DIM}' in labels else None)
        if amp_col is not None:
            # Amplitude + phase in degrees
            a       = data[:, labels.index(amp_col)] * amp_scale[d]
            phi_rad = np.radians(data[:, labels.index(f'PH{DIM}')])
            AB[d, :, 0] = -a * np.sin(phi_rad)
            AB[d, :, 1] =  a * np.cos(phi_rad)
        else:
            # Legacy Cartesian A, B coefficients (always fractional)
            AB[d, :, 0] = data[:, labels.index(f'A{DIM}')]
            AB[d, :, 1] = data[:, labels.index(f'B{DIM}')]

    snr    = data[:, labels.index('SNR')]
    active = data[:, labels.index('ACTIVE')].astype(bool)

    return hkls, AB, active, snr, cell1, cell2, dimensions, rmsd


# ══════════════════════════════════════════════════════════════════════════════
# Origin search (handles symmetry-equivalent origins in polar/centred SGs)
# ══════════════════════════════════════════════════════════════════════════════

# Allowed fractional origin shifts per space group (from CCP4 convention).
# 'x','y','z' = polar axis (continuous free parameter).
# 'x=y=z'     = R3 special (all three equal).
_ORIGINS_TABLE = {
    # P1 — fully free
    'P1': [('x', 'y', 'z')],
    # Monoclinic — polar y
    **{sg: [(0, 'y', 0), (0, 'y', .5), (.5, 'y', 0), (.5, 'y', .5)]
       for sg in ('P2', 'P21', 'C2', 'B2', 'A2', 'I2')},
    # Orthorhombic / some cubic — all discrete
    **{sg: [(0,0,0),(0,0,.5),(0,.5,0),(0,.5,.5),(.5,0,0),(.5,0,.5),(.5,.5,0),(.5,.5,.5)]
       for sg in ('P222','P2221','P2212','P2122','P21212','P21221','P22121','P212121',
                  'C2221','C222','I222','I212121','F432','F4132')},
    # F-centred — 16 origins
    **{sg: [(0,0,0),(0,0,.5),(0,.5,0),(0,.5,.5),(.5,0,0),(.5,0,.5),(.5,.5,0),(.5,.5,.5),
            (.25,.25,.25),(.25,.25,.75),(.25,.75,.25),(.25,.75,.75),
            (.75,.25,.25),(.75,.25,.75),(.75,.75,.25),(.75,.75,.75)]
       for sg in ('F222', 'F23')},
    # Tetragonal — polar z (2 xy positions)
    **{sg: [(0, 0, 'z'), (.5, .5, 'z')]
       for sg in ('P4','P41','P42','P43','I4','I41')},
    # Tetragonal 422 — 4 discrete origins
    **{sg: [(0,0,0),(0,0,.5),(.5,.5,0),(.5,.5,.5)]
       for sg in ('P422','P4212','P4122','P41212','P4222','P42212','P4322','P43212',
                  'I422','I4122')},
    # Trigonal P — polar z (3 xy positions)
    **{sg: [(0, 0, 'z'), (1/3, 2/3, 'z'), (2/3, 1/3, 'z')]
       for sg in ('P3', 'P31', 'P32')},
    # R3 — special: (t, t, t) diagonal
    'R3': [('x=y=z',)],
    # H3 — polar z (3 xy positions)
    **{sg: [(0, 0, 'z'), (1/3, 2/3, 'z'), (2/3, 1/3, 'z')]
       for sg in ('H3',)},
    # P312 family
    **{sg: [(0,0,0),(0,0,.5),(1/3,2/3,0),(2/3,1/3,0),(1/3,2/3,.5),(2/3,1/3,.5)]
       for sg in ('P312','P3112','P3212')},
    # P321 / P622 family
    **{sg: [(0,0,0),(0,0,.5)]
       for sg in ('P321','P3121','P3221','P622','P6122','P6522','P6222','P6422','P6322')},
    # R32 — 2 origins
    'R32': [(0,0,0), (.5,.5,.5)],
    # H32 — 6 origins
    'H32': [(0,0,0),(0,0,.5),(1/3,2/3,2/3),(2/3,1/3,1/3),(1/3,2/3,1/6),(2/3,1/3,5/6)],
    # Hexagonal polar z
    **{sg: [(0, 0, 'z')]
       for sg in ('P6','P61','P65','P62','P64','P63')},
    # Cubic with 2 origins
    **{sg: [(0,0,0), (.5,.5,.5)]
       for sg in ('P23','P213','P432','P4232','I432','P4332','P4132','I4132')},
}


def _sg_origins_key(sg_name):
    """Map gemmi space group name to _ORIGINS_TABLE key (strip spaces)."""
    k = sg_name.replace(' ', '').replace('_', '')
    # Handle H/R notation variants gemmi may produce
    alt = {'H3': 'H3', 'H32': 'H32', 'R3': 'R3', 'R32': 'R32'}
    return alt.get(k, k)


def _expand_origin_entries(entries, n_polar=12):
    """Expand polar placeholders to concrete (dx,dy,dz) float tuples."""
    polar_vals = [i / n_polar for i in range(n_polar)]
    result = []
    for entry in entries:
        if len(entry) == 1 and entry[0] == 'x=y=z':
            for t in polar_vals:
                result.append((t, t, t))
            continue
        coords = list(entry)
        def _vals(c):
            return polar_vals if isinstance(c, str) else [float(c)]
        for dx in _vals(coords[0]):
            for dy in _vals(coords[1]):
                for dz in _vals(coords[2]):
                    result.append((dx, dy, dz))
    return result


def _op_idx(uid):
    """Return the integer operator index from a uid ending in _opN."""
    i = uid.rfind('_op')
    return int(uid[i + 3:]) if i >= 0 else 0


def _base_uid(uid):
    """Strip _opN suffix from a uid."""
    i = uid.rfind('_op')
    return uid[:i] if i >= 0 else uid


def _cross_op_ca_shifts(atoms1_op0, atoms2, k):
    """CA shifts: atoms2 op_k minus atoms1 op_0, matched by base uid.

    Returns ndarray (N_ca, 3) of fractional shifts (minimum-image convention).
    """
    op_tag = f'_op{k}'
    a2_map = {_base_uid(a['uid']): a
              for a in atoms2 if a['is_ca'] and a['uid'].endswith(op_tag)}
    shifts = []
    for a1 in atoms1_op0:
        if not a1['is_ca']:
            continue
        a2 = a2_map.get(_base_uid(a1['uid']))
        if a2 is None:
            continue
        dx = a2['x'] - a1['x'];  dx -= round(dx)
        dy = a2['y'] - a1['y'];  dy -= round(dy)
        dz = a2['z'] - a1['z'];  dz -= round(dz)
        shifts.append([dx, dy, dz])
    return np.array(shifts, dtype=float) if shifts else np.zeros((0, 3))


def _relabel_ops(atoms, k_offset, n_ops):
    """Permute operator labels: old _opJ → new _op_{(J - k_offset) % n_ops}."""
    if k_offset == 0:
        return atoms
    result = []
    for a in atoms:
        uid = a['uid']
        i = uid.rfind('_op')
        if i >= 0:
            j = int(uid[i + 3:])
            a = {**a, 'uid': uid[:i] + f'_op{(j - k_offset) % n_ops}'}
        result.append(a)
    return result


# Alternative-indexing candidates by crystal system (fractional coords, det=+1).
# These are proper rotations in the Laue-group holohedry that may be absent
# from lower-symmetry space groups.
def _metric_tensor(cell):
    """Lattice metric tensor G = O^T·O in cartesian, where O is the
    fractional→cartesian basis matrix (columns are a, b, c vectors).
    Accepts a (a,b,c,α,β,γ) tuple or a gemmi.UnitCell."""
    if not isinstance(cell, gemmi.UnitCell):
        cell = gemmi.UnitCell(*cell)
    O = np.zeros((3, 3))
    for i, b in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
        p = cell.orthogonalize(gemmi.Fractional(*b))
        O[:, i] = [p.x, p.y, p.z]
    return O.T @ O


_ALT_SYSTEM_GROUPS = {
    'triclinic':    ('triclinic',),
    'monoclinic':   ('triclinic', 'monoclinic'),
    'orthorhombic': ('triclinic', 'monoclinic', 'orthorhombic'),
    'tetragonal':   ('triclinic', 'monoclinic', 'orthorhombic', 'tetragonal'),
    'trig_hex':     ('triclinic', 'monoclinic', 'orthorhombic',
                     'trigonal', 'hexagonal'),
    'cubic':        ('triclinic', 'monoclinic', 'orthorhombic',
                     'tetragonal', 'cubic'),
}


def _cell_params_from_G(G):
    a  = float(np.sqrt(G[0, 0]))
    b  = float(np.sqrt(G[1, 1]))
    c  = float(np.sqrt(G[2, 2]))
    al = float(np.degrees(np.arccos(np.clip(G[1, 2] / (b * c), -1, 1))))
    be = float(np.degrees(np.arccos(np.clip(G[0, 2] / (a * c), -1, 1))))
    ga = float(np.degrees(np.arccos(np.clip(G[0, 1] / (a * b), -1, 1))))
    return (a, b, c, al, be, ga)


def _crystal_system_of_params(params, deg=1.5, atol_a=1.0):
    a, b, c, al, be, ga = params
    is_90  = lambda x: abs(x - 90)  < deg
    is_120 = lambda x: abs(x - 120) < deg
    eq     = lambda x, y: abs(x - y) < atol_a
    if is_90(al) and is_90(be) and is_90(ga):
        if eq(a, b) and eq(b, c):    return 'cubic'
        if eq(a, b) or eq(b, c) or eq(a, c): return 'tetragonal'
        return 'orthorhombic'
    if is_90(al) and is_90(be) and is_120(ga) and eq(a, b):
        return 'trig_hex'
    if sum(1 for x in (al, be, ga) if not is_90(x)) == 1:
        return 'monoclinic'
    return 'triclinic'


import functools as _functools


@_functools.lru_cache(maxsize=64)
def _get_altindex_ops_cached(sg_name, cell_tuple, max_M, det_max, metric_tol_rel):
    """LRU-cached worker for _get_altindex_ops.  Cell is a tuple for
    hashability; cell=None is encoded as None.  ~minute-long enumeration
    runs only once per (sg, cell, max_M, det_max) within a process."""
    return _get_altindex_ops_impl(sg_name, cell_tuple, max_M, det_max,
                                   metric_tol_rel)


def _get_altindex_ops(sg_name, cell=None, max_M=1, det_max=1,
                       metric_tol_rel=1e-6):
    """metric_tol_rel : relative tolerance for metric preservation,
    expressed as a fraction of max(|G|).  Default 1e-6 admits only ops
    that EXACTLY preserve the lattice metric (the rigorous altindex set).
    A larger value (e.g. 0.05) admits ops that preserve a NEARBY
    higher-symmetry metric — useful for pseudo-symmetric cells where
    the deposited cell is slightly distorted from an actual orthorhombic /
    tetragonal / etc. lattice and the natural alt-cell ops fall just
    outside strict tolerance.  resolve_altindex uses 0.05 in the
    cross-cell (post-stretch) branch so a 180°-about-z rotation in a
    nearly-orthorhombic monoclinic, for example, is found and can be
    applied as a discrete reindexing (preserving experimental Fobs)
    rather than forcing a Fcalc-only fallback."""
    if cell is None:
        cell_tuple = None
    elif hasattr(cell, 'a'):       # gemmi.UnitCell
        cell_tuple = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    else:
        cell_tuple = tuple(cell)
    return _get_altindex_ops_cached(sg_name, cell_tuple, max_M, det_max,
                                     metric_tol_rel)


@_functools.lru_cache(maxsize=1)
def _catalog_ops_by_cs():
    """Pre-extract (R_int, t_frac) arrays for every catalog SG, grouped
    by crystal system.  Saves ~100× over per-call iteration of
    gemmi.spacegroup_table() + .operations() + per-op array construction.
    Returns dict: crystal_system → (R_all (n,3,3) float, t_all (n,3) float).
    """
    by_cs = {}
    for cand in gemmi.spacegroup_table():
        Rs, ts = [], []
        for op in cand.operations():
            Rs.append(np.array(op.rot, dtype=float) / op.DEN)
            ts.append(np.array(op.tran, dtype=float) / op.DEN)
        cs = cand.crystal_system_str()
        prev = by_cs.get(cs)
        if prev is None:
            by_cs[cs] = ([Rs], [ts])
        else:
            prev[0].append(Rs); prev[1].append(ts)
    # Concatenate per-cs into single arrays
    out = {}
    for cs, (Rs_list, ts_list) in by_cs.items():
        R_flat = np.concatenate([np.array(r) for r in Rs_list], axis=0)
        t_flat = np.concatenate([np.array(t) for t in ts_list], axis=0)
        out[cs] = (R_flat, t_flat)
    return out


def _get_altindex_ops_impl(sg_name, cell=None, max_M=1, det_max=1,
                             metric_tol_rel=1e-6):
    """All proper-rotation (R, t) ops that preserve the cell's lattice metric
    and are NOT already in this SG.  These are the "altindex" candidates:
    rotations + translations that map the lattice to itself but aren't
    already symmetries of the SG.

    Uses a basis-change enumeration ported from ~/projects/origins/claude/
    gemmi_altindex.py: integer M matrices with entries in {-max_M..max_M}
    and |det(M)| ≤ det_max (defaults 2 and 4); for each M, build the alt-cell
    metric G_alt = M G Mᵀ, identify its crystal system, look up catalog SGs
    matching that system, then conjugate each catalog symop back to the
    original frame as R_old = Mᵀ R_alt M⁻ᵀ, t_old = Mᵀ t_alt.  Keep ops
    where R_old is integer with entries in {-1,0,1} and t_old has
    denominator in {1,2,3} (centring fractions of standard SG settings).

    Returns a list of (R_float, t_float) tuples.  Empty if the SG already
    spans its lattice holohedry.  Identity is excluded.

    A `cell` is needed because metric-preservation depends on the lattice
    parameters.  If cell is None we fall back to an isotropic generic
    metric — imperfect but harmless when the SG strongly constrains the
    cell anyway.
    """
    import itertools
    from fractions import Fraction
    sg = gemmi.find_spacegroup_by_name(sg_name)
    cell_for_G = cell if cell is not None else (10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
    G = _metric_tensor(cell_for_G)

    # SG's existing proper rotations (drop these from the altindex set)
    sg_R = set()
    for op in sg.operations():
        R = (np.array(op.rot, dtype=int) // op.DEN)
        if int(round(np.linalg.det(R))) == 1:
            sg_R.add(tuple(R.flatten()))

    catalog = _catalog_ops_by_cs()
    tol_G = float(metric_tol_rel) * float(np.max(np.abs(G)))
    all_ops = set()

    # Pre-filter M matrices by det in [-det_max..det_max] excluding 0.
    # At max_M=1 det_max=1 this drops 19683 → ~3000 (only |det|=1 survive).
    M_candidates = []
    for M_flat in itertools.product(range(-max_M, max_M + 1), repeat=9):
        M = np.array(M_flat, dtype=float).reshape(3, 3)
        d = int(round(np.linalg.det(M)))
        if d == 0 or abs(d) > det_max:
            continue
        M_candidates.append(M)

    for M in M_candidates:
        G_alt = M @ G @ M.T
        try:
            params = _cell_params_from_G(G_alt)
        except (ValueError, FloatingPointError):
            continue
        if any(p < 1.0 for p in params[:3]):
            continue
        try:
            M_inv_T = np.linalg.inv(M.T)
        except np.linalg.LinAlgError:
            continue
        cs = _crystal_system_of_params(params)

        # Vectorized per-M: build R_all/t_all by concatenating per-CS arrays
        R_all_list, t_all_list = [], []
        for sub_cs in _ALT_SYSTEM_GROUPS[cs]:
            entry = catalog.get(sub_cs)
            if entry is None: continue
            R_all_list.append(entry[0])
            t_all_list.append(entry[1])
        if not R_all_list:
            continue
        R_alt = np.concatenate(R_all_list, axis=0)   # (N, 3, 3)
        t_alt = np.concatenate(t_all_list, axis=0)   # (N, 3)

        # Bulk conjugate: R_old[n] = M.T @ R_alt[n] @ M_inv_T
        R_old = np.einsum('ij,njk,kl->nil', M.T, R_alt, M_inv_T)
        # Filter integer R with |entries| ≤ 1
        R_round = np.round(R_old)
        ok_int = np.all(np.abs(R_old - R_round) < 1e-4, axis=(1, 2))
        ok_small = np.all(np.abs(R_round) <= 1.0, axis=(1, 2))
        keep = ok_int & ok_small
        if not keep.any():
            continue
        R_int = R_round[keep].astype(int)            # (M, 3, 3)
        t_old = (M.T @ t_alt[keep].T).T              # (M, 3)
        det_R = np.round(np.linalg.det(R_int))
        keep2 = det_R == 1
        if not keep2.any():
            continue
        R_int = R_int[keep2]; t_old = t_old[keep2]

        # Metric-preservation: R_int.T @ G @ R_int ≈ G  (fast manual check)
        R_int_f = R_int.astype(float)
        prod = np.einsum('nji,jk,nkl->nil', R_int_f, G, R_int_f)
        keep3 = (np.max(np.abs(prod - G).reshape(len(R_int), -1), axis=1)
                 < tol_G)
        if not keep3.any():
            continue
        R_int = R_int[keep3]; t_old = t_old[keep3]

        # Translation denominator ∈ {1, 2, 3}: equivalent to 2t OR 3t integer.
        t_mod = t_old - np.floor(t_old + 1e-9)        # wrap to [0,1)
        t2 = np.abs(t_mod * 2.0 - np.round(t_mod * 2.0)) < 1e-3
        t3 = np.abs(t_mod * 3.0 - np.round(t_mod * 3.0)) < 1e-3
        per_comp_ok = t2 | t3
        keep4 = np.all(per_comp_ok, axis=1)
        if not keep4.any():
            continue
        R_int = R_int[keep4]; t_mod = t_mod[keep4]
        # Snap t_mod to its nearest /6 grid value to canonicalize for dedup
        t_snap = np.round(t_mod * 6.0) / 6.0
        # Wrap any 1.0 back to 0.0
        t_snap = np.round(t_snap, 6) % 1.0

        for i in range(len(R_int)):
            R_tup = tuple(int(v) for v in R_int[i].flatten())
            t_tup = tuple(float(round(v, 6)) for v in t_snap[i])
            if R_tup in sg_R and all(abs(v) < 1e-6 for v in t_tup):
                continue                              # plain SG identity-translation op
            all_ops.add((R_tup, t_tup))

    return [(np.array(R_tup, dtype=float).reshape(3, 3),
             np.array(t_tup, dtype=float))
            for R_tup, t_tup in sorted(all_ops)]


def _apply_rotation_to_atoms(atoms, R):
    """Apply 3×3 rotation R (fractional coords) to atom positions; uid unchanged.

    NOTE: rotates every P1 atom uniformly.  For non-trivial R that is not in the
    SG, the resulting non-_op0 atoms are NOT at symop_k(rotated_asu); use
    _reexpand_atoms_with_rotation when symmetry-consistent atoms are needed.
    """
    result = []
    for a in atoms:
        xyz = R @ np.array([a['x'], a['y'], a['z']])
        result.append({**a, 'x': xyz[0], 'y': xyz[1], 'z': xyz[2]})
    return result


def _reexpand_atoms_with_rotation(atoms_full, R_alt, sg_name):
    """Rotate the ASU (uid endswith _op0) atoms by R_alt in fractional coords,
    then re-expand to P1 using sg's proper symops.

    Returns a new atoms list whose _op_k positions are symop_k(R_alt @ asu),
    matching the convention of expand_to_p1 — so atoms1 (which is built that
    way) and the re-expanded atoms2 are directly comparable by uid.
    """
    asu = [a for a in atoms_full if a['uid'].endswith('_op0')]
    proper_ops = get_proper_symops(sg_name)
    out = []
    R = np.asarray(R_alt, dtype=float)
    for k, (Rop, top) in enumerate(proper_ops):
        for a in asu:
            f      = R @ np.array([a['x'], a['y'], a['z']])
            f_op   = Rop @ f + top
            uid    = _base_uid(a['uid']) + f'_op{k}'
            out.append({**a, 'x': f_op[0], 'y': f_op[1], 'z': f_op[2],
                        'uid': uid})
    return out


def _find_best_origin(atoms1_op0, atoms2, sg_name, cell=None, n_polar=12):
    """Find best (origin_shift, symop_idx, altindex_R, improved).

    Searches (allowed_shift, operator_k) combinations.  If no improvement,
    tries each altindex rotation (proper ops in the Laue holohedry absent from
    the SG) on atoms2 and repeats the search.

    Returns
    -------
    best_d      : ndarray (3,) fractional shift to apply to atoms1
    best_k      : int         operator index in atoms2 that matches atoms1 op0
    best_R_alt  : ndarray (3,3) or None — altindex rotation applied to atoms2
    improved    : bool
    """
    key = _sg_origins_key(sg_name)
    entries = _ORIGINS_TABLE.get(key)
    if entries is None:
        candidates = np.array([(dx, dy, dz)
                                for dx in (0, .5) for dy in (0, .5) for dz in (0, .5)],
                               dtype=float)
    else:
        candidates = np.array(_expand_origin_entries(entries, n_polar), dtype=float)

    n_ops = 1 + max((_op_idx(a['uid']) for a in atoms2), default=0)

    sh0 = _cross_op_ca_shifts(atoms1_op0, atoms2, 0)
    if len(sh0) < 3:
        return np.zeros(3), 0, None, False
    ref_score = float(np.median(np.sqrt(np.sum(sh0**2, axis=1))))

    def _search(atoms2_try):
        """Inner search over all (k, candidate) pairs; returns (score, d, k)."""
        best_s = ref_score
        best_d = np.zeros(3)
        best_k = 0
        for k in range(n_ops):
            ca_sh = _cross_op_ca_shifts(atoms1_op0, atoms2_try, k)
            if len(ca_sh) < 3:
                continue
            new_sh  = ca_sh[None, :, :] - candidates[:, None, :]
            new_sh -= np.round(new_sh)
            medians = np.median(np.sqrt(np.sum(new_sh**2, axis=2)), axis=1)
            i = int(np.argmin(medians))
            if medians[i] < best_s:
                best_s = medians[i]; best_d = candidates[i]; best_k = k
        return best_s, best_d, best_k

    best_score, best_d, best_k = _search(atoms2)
    best_R_alt = None

    if not (best_score < ref_score * 0.9):
        for R_alt, _t_alt in _get_altindex_ops(sg_name, cell=cell):
            # _find_best_origin's caller applies altindex via
            # _reexpand_atoms_with_rotation which takes R only; the
            # translation is absorbed by the existing origin candidates
            # loop above.
            atoms2_alt = _reexpand_atoms_with_rotation(atoms2, R_alt, sg_name)
            s, d, k = _search(atoms2_alt)
            if s < best_score:
                best_score = s; best_d = d; best_k = k; best_R_alt = R_alt

    improved = bool(best_score < ref_score * 0.9)
    return best_d, best_k, best_R_alt, improved


def find_origin_alignment(pdb1_path, pdb2_path,
                          altloc_filter=False, altloc_fallback=1.0):
    """Run the origin / altindex / symop search for aligning pdb1 onto pdb2.

    Returns (shift, best_k, best_R_alt, improved).  When improved=False the
    other values are zero/None and the caller should skip alignment.
    """
    atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    atoms2, cell2, sg2 = expand_to_p1(pdb2_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)
    op0_ca_sh = fitme[[u.endswith('_op0') and m
                        for u, m in zip(uids, ca_mask)], 3:6]
    if not len(op0_ca_sh) or np.max(np.sqrt(np.sum(op0_ca_sh**2, axis=1))) <= 0.1:
        return np.zeros(3), 0, None, False
    shift, best_k, best_R_alt, improved = _find_best_origin(
        [a for a in atoms1 if a['uid'].endswith('_op0')], atoms2, sg1,
        cell=cell1)
    if improved and (np.any(np.abs(shift) > 1e-6) or best_k != 0
                     or best_R_alt is not None):
        return shift, best_k, best_R_alt, True
    return np.zeros(3), 0, None, False


# ══════════════════════════════════════════════════════════════════════════════
# Public API
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class BendResult:
    hkls:       object    # ndarray (N_hkls, 3) int
    AB:         object    # ndarray (N_dims, N_hkls, 2)
    active:     object    # ndarray (N_hkls,) bool
    snr:        object    # ndarray (N_hkls,) float
    fitme:      object    # ndarray (N_atoms, 8) — xf yf zf dxf dyf dzf do dB
    ca_mask:    object    # ndarray (N_atoms,) bool
    rmsd:       float
    cell1:      tuple     # moving crystal cell (a,b,c,alpha,beta,gamma)
    cell2:      tuple     # reference crystal cell
    dimensions: str       # e.g. 'xyz'


def bend_fit(pdb1_path, pdb2_path, nhkls=30, fitreso=None, drop_snr=1.0,
             dimensions='xyz', geotest=False, frac=1.0, verbose=True,
             use_symm=True, atom_sel='auto', outlier_sigma=3.0, b_sigma=None,
             altloc_filter=False, altloc_fallback=1.0):
    """Fit a Fourier shift field between pdb1 (moving) and pdb2 (reference).

    HKL selection: pass either fitreso (d-spacing cutoff in Å) or nhkls (count).
    When nhkls is given the effective fitreso is printed for reference.

    Returns BendResult with fitted coefficients and RMSD.
    """
    t0 = time.time()

    if verbose:
        print(f"expanding {pdb1_path} from", end=' ', flush=True)
    atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    if verbose:
        print(f"{sg1} to P1")
        print(f"expanding {pdb2_path} from", end=' ', flush=True)
    atoms2, cell2, sg2 = expand_to_p1(pdb2_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    if verbose:
        print(f"{sg2} to P1")

    fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)

    # ── Origin search — only when ASU atoms show large shifts
    op0_ca_sh = fitme[[u.endswith('_op0') and m
                        for u, m in zip(uids, ca_mask)], 3:6]
    if len(op0_ca_sh) and np.max(np.sqrt(np.sum(op0_ca_sh**2, axis=1))) > 0.1:
        shift, best_k, best_R_alt, improved = _find_best_origin(
            [a for a in atoms1 if a['uid'].endswith('_op0')], atoms2, sg1,
            cell=cell1)
        if improved and (np.any(np.abs(shift) > 1e-6) or best_k != 0
                         or best_R_alt is not None):
            msg = f"  origin shift: ({shift[0]:.4f}, {shift[1]:.4f}, {shift[2]:.4f}) for {sg1}"
            if best_R_alt is not None:
                atoms2 = _reexpand_atoms_with_rotation(atoms2, best_R_alt, sg1)
                msg += f"  +altindex"
            if best_k != 0:
                n_ops = 1 + max(_op_idx(a['uid']) for a in atoms2)
                atoms2 = _relabel_ops(atoms2, best_k, n_ops)
                msg += f"  +symop_idx={best_k}"
            print(msg, flush=True)
            atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                               altloc_filter=altloc_filter,
                                               altloc_fallback=altloc_fallback,
                                               origin_shift=shift)
            fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)

    # ── Atom selection ────────────────────────────────────────────────────────
    # 'auto'     : all atoms, MAD outlier rejection applied
    # 'all'      : all atoms, no filtering
    # 'backbone' : N/CA/C/O only
    # 'ca'       : CA only
    if atom_sel == 'auto':
        fitme, ca_mask, uids, bfacs, _ = reject_outliers(
            fitme, ca_mask, uids, cell1, mad_sigma=outlier_sigma,
            bfacs=bfacs, b_sigma=b_sigma)
        fit_mask = np.ones(len(fitme), dtype=bool)
    elif atom_sel == 'ca':
        fit_mask = ca_mask.copy()
    elif atom_sel == 'backbone':
        fit_mask = np.array([uid.split('_')[-2] in _BACKBONE for uid in uids])
    else:   # 'all'
        fit_mask = np.ones(len(fitme), dtype=bool)

    if verbose and atom_sel != 'all':
        print(f"  atom_sel={atom_sel!r}: {fit_mask.sum()} atoms used for fitting")

    fitme_fit   = fitme[fit_mask]
    ca_mask_fit = ca_mask[fit_mask]
    uids_fit    = [u for u, k in zip(uids, fit_mask) if k]

    # Report largest CA fractional shift and RMS CA shift in Å after origin fix
    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    ca_orth   = frac_to_orth(ca_shifts, cell1)
    ca_rmsd   = float(np.sqrt(np.mean(np.sum(ca_orth**2, axis=1))))
    print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}")
    print(f"CA RMSD after origin fix: {ca_rmsd:.3f} Å  ({ca_mask.sum()} CA pairs)")
    if ca_rmsd > 10.0:
        raise ValueError(f"CA RMSD after origin fix {ca_rmsd:.3f} Å > 10.0 Å. "
                         "Check structure alignment.")
    if ca_mags[max_i] > 0.35:
        raise ValueError(f"Largest fractional CA shift {ca_mags[max_i]:.4f} > 0.35 cells. "
                         "Check structure alignment.")

    _FINE_RESO = 1.0   # internal pool limit for nhkls mode

    if fitreso is not None:
        # fitreso mode: generate only HKLs with d >= fitreso
        all_hkls  = generate_hkls(cell1, fitreso)
        hkls      = all_hkls
        nhkls_use = len(hkls)
        if verbose:
            print(f"fitreso={fitreso} A → {nhkls_use - 1} HKLs "
                  f"for each dimension: {dimensions}")
    else:
        # nhkls mode: generate fine pool, take first nhkls, report effective fitreso
        all_hkls  = generate_hkls(cell1, _FINE_RESO)
        nhkls_use = min(nhkls, len(all_hkls))
        hkls      = all_hkls[:nhkls_use]
        if len(hkls) > 1 and verbose:
            Gs     = _reciprocal_metric(cell1)
            last   = hkls[-1].astype(float)
            inv_d2 = float(last @ Gs @ last)
            eff_fr = 1.0 / np.sqrt(inv_d2) if inv_d2 > 0 else np.inf
            print(f"nhkls={nhkls_use} → fitreso={eff_fr:.2f} A "
                  f"for each dimension: {dimensions}")

    # Dimension mapping: x=col3, y=col4, z=col5, o=col6, B=col7
    dim_col = {'x': 3, 'y': 4, 'z': 5, 'o': 6, 'B': 7}
    dim_list = [d for d in dimensions if d in dim_col]
    dim_indices = [dim_col[d] for d in dim_list]
    shifts      = fitme_fit[:, dim_indices]   # (N_fit, N_dims)
    frac_coords = fitme_fit[:, :3]

    # Use symmetry-constrained fitting when all xyz dims are present and sg > P1
    proper_ops = get_proper_symops(sg1)
    do_symm = (use_symm and len(proper_ops) > 1
               and all(d in dim_list for d in 'xyz'))

    dt = time.time() - t0
    print(f"\n=================  Starting nhkls {nhkls_use} fit at {dt:.0f} s\n")

    if do_symm:
        canon_hkls, hkl_to_canon, hkl_all_ops = assign_canonical(hkls, proper_ops)
        n_canon = len(canon_hkls)
        print(f"{n_canon} canonical HKLs ({len(proper_ops)} proper symops, "
              f"reduction {nhkls_use}/{n_canon}={nhkls_use/max(n_canon,1):.1f}x)")
        # Use only ASU atoms (op0): build_design_matrix_symm applies all ops
        # internally, so feeding all P1 copies would double-count the data.
        asu_mask    = np.array([uid.endswith('_op0') for uid in uids_fit])
        frac_asu    = fitme_fit[asu_mask, :3]
        X_symm = build_design_matrix_symm(frac_asu, canon_hkls, proper_ops,
                                           hkl_to_canon, hkl_all_ops)
        xyz_col = [dim_list.index(d) for d in 'xyz']
        shifts_xyz  = shifts[asu_mask][:, xyz_col]   # (N_asu, 3)
        shifts_flat = shifts_xyz.ravel()              # (3*N_asu,) interleaved
        params, active_c, snr_c = fit_lstsq_symm(
            X_symm, shifts_flat, drop_snr=drop_snr)
        AB_xyz = expand_ab_canon(params, active_c, canon_hkls, hkls,
                                 hkl_to_canon, hkl_all_ops, proper_ops)
        # Build full AB (N_dims, N_hkls, 2) — xyz from symm fit
        AB = np.zeros((len(dim_list), nhkls_use, 2))
        for new_i, d in enumerate(dim_list):
            if d in 'xyz':
                AB[new_i] = AB_xyz['xyz'.index(d)]
        # Propagate active/snr from canonical → Friedel-unique
        active = np.array([active_c[hkl_to_canon[j]] for j in range(nhkls_use)])
        snr    = np.array([snr_c[hkl_to_canon[j]]    for j in range(nhkls_use)])
    else:
        X = build_design_matrix(frac_coords, hkls)
        AB, active, snr = fit_lstsq(X, shifts, drop_snr=drop_snr)

    n_dropped = nhkls_use - int(active.sum())
    if n_dropped > 0:
        print(f"SNR-pruned {n_dropped} params")

    # Extract x,y,z AB for RMSD and PDB output
    xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
    AB_xyz_all = AB[xyz_dims]   # (≤3, N_hkls, 2)
    # Pad to (3, N_hkls, 2) if fewer than 3 dims were fitted
    if len(xyz_dims) < 3:
        AB_xyz = np.zeros((3, nhkls_use, 2))
        for new_i, old_i in enumerate(xyz_dims):
            AB_xyz[new_i] = AB_xyz_all[new_i]
    else:
        AB_xyz = AB_xyz_all

    rmsd = rmsd_ca(fitme, ca_mask, hkls, AB_xyz, cell2, frac)
    n_ca = int(ca_mask.sum())
    n_active = int(active.sum())
    dt2 = time.time() - t0
    print(f"RMSD(CA) = {rmsd:.3f} A  ({n_ca} CA atoms, {n_active} active HKLs, {dt2:.0f} s)")

    return BendResult(
        hkls=hkls, AB=AB, active=active, snr=snr,
        fitme=fitme, ca_mask=ca_mask, rmsd=rmsd,
        cell1=cell1, cell2=cell2, dimensions=''.join(dim_list),
    )


def bend_fit_progressive(pdb1_path, pdb2_path,
                         fitreso_start=20.0, fitreso_end=7.0,
                         drop_snr=0.0, batch_hkls=100, od_margin=1.5,
                         outlier_sigma=2.5, b_sigma=3.0,
                         use_symm=True, dimensions='xyz', atom_sel='all',
                         frac=1.0, verbose=True,
                         altloc_filter=False, altloc_fallback=1.0,
                         iter_callback=None, max_canon=None,
                         precomputed_origin=None,
                         iteration_schedule=None,
                         bound_by_obs=False,
                         pnn_mode=None,
                         holdout_frac=0.0, holdout_seed=0,
                         cv_callback=None):
    """Fit a Fourier shift field progressively, admitting HKLs coarsest-first.

    HKLs are admitted in batches from fitreso_start (coarsest, large d) down to
    fitreso_end (finest, small d).  After each batch the fit residuals are
    MAD-filtered to adaptively update the active atom set for the next iteration.
    Stops when the overdetermination ratio (N_asu_active × 3) / (N_canon × 6)
    drops below od_margin, or fitreso_end is reached.

    Defaults optimised on DHFR 1rx1→1rx2 (P2₁2₁2₁, 636 CA pairs):
      batch_hkls=100, fitreso_end=7.0, drop_snr=0.0, outlier_sigma=2.5, b_sigma=3.0
      → CA RMSD 0.071 Å, lig diff-map 5.76 σ in ~190 s.
    The OD limit naturally stops around 7 Å for typical protein sizes, so
    fitreso_end=7.0 avoids wasteful HKL pool generation beyond the stopping point.
    drop_snr=0 skips per-iteration SVD pruning; residual-MAD atom filtering
    handles outlier rejection instead, which is faster and equally effective.

    Parameters
    ----------
    fitreso_start : float — d-spacing (Å) of the initial coarse batch (default 20 Å)
    fitreso_end   : float — finest d-spacing to pool (default 7 Å; OD limit typically
                            stops here anyway — go finer only if cell is very large)
    batch_hkls    : int   — canonical HKLs to add per iteration (default 100)
    od_margin     : float — stop when OD ratio drops below this (default 1.5)
    drop_snr      : float — SVD SNR pruning threshold per iteration; 0 disables
                            (default 0 — residual-MAD handles outliers instead)
    outlier_sigma : float — MAD sigma for initial shift rejection and per-iter
                            residual rejection (default 2.5)
    b_sigma       : float — B-factor outlier rejection threshold in σ (default 3.0;
                            None disables B-factor filtering)

    Returns a BendResult identical to bend_fit.
    """
    t0 = time.time()
    pool_reso = fitreso_end if fitreso_end is not None else 2.0

    if verbose:
        print(f"expanding {pdb1_path} from", end=' ', flush=True)
    atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    if verbose:
        print(f"{sg1} to P1", flush=True)
        print(f"expanding {pdb2_path} from", end=' ', flush=True)
    atoms2, cell2, sg2 = expand_to_p1(pdb2_path,
                                       altloc_filter=altloc_filter,
                                       altloc_fallback=altloc_fallback)
    if verbose:
        print(f"{sg2} to P1", flush=True)

    # If atoms1 and atoms2 use different single-chain IDs but the same
    # residue sequence (e.g. 'A' vs 'X' in older PDB-REDO outputs),
    # remap atoms2's chain to atoms1's BEFORE the origin search — which
    # also matches by uid and would otherwise hit zero pairs.
    atoms2 = _maybe_remap_single_chain(atoms1, atoms2, verbose=verbose)

    if verbose:
        print("matching atoms...", end=' ', flush=True)
    fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)
    if verbose:
        print(f"{ca_mask.sum()} CA pairs", flush=True)

    # ── Origin alignment — precomputed (skip search) or auto-search ──────────
    if precomputed_origin is not None:
        shift, best_k, best_R_alt = precomputed_origin
        shift = np.asarray(shift, dtype=float)
        applied = (np.any(np.abs(shift) > 1e-6) or best_k != 0
                   or best_R_alt is not None)
        if applied:
            msg = (f"  origin (precomputed): "
                   f"({shift[0]:.4f}, {shift[1]:.4f}, {shift[2]:.4f}) for {sg1}")
            if best_R_alt is not None:
                atoms2 = _reexpand_atoms_with_rotation(atoms2, best_R_alt, sg1)
                msg += "  +altindex"
            if best_k != 0:
                n_ops = 1 + max(_op_idx(a['uid']) for a in atoms2)
                atoms2 = _relabel_ops(atoms2, best_k, n_ops)
                msg += f"  +symop_idx={best_k}"
            if verbose:
                print(msg, flush=True)
            if np.any(np.abs(shift) > 1e-6):
                atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                                   altloc_filter=altloc_filter,
                                                   altloc_fallback=altloc_fallback,
                                                   origin_shift=shift)
            fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)
            if verbose:
                print(f"  after origin fix: {ca_mask.sum()} CA pairs", flush=True)
    else:
        op0_ca_sh = fitme[[u.endswith('_op0') and m
                            for u, m in zip(uids, ca_mask)], 3:6]
        if len(op0_ca_sh) and np.max(np.sqrt(np.sum(op0_ca_sh**2, axis=1))) > 0.1:
            shift, best_k, best_R_alt, improved = _find_best_origin(
                [a for a in atoms1 if a['uid'].endswith('_op0')], atoms2, sg1,
                cell=cell1)
            if improved and (np.any(np.abs(shift) > 1e-6) or best_k != 0
                             or best_R_alt is not None):
                msg = f"  origin shift: ({shift[0]:.4f}, {shift[1]:.4f}, {shift[2]:.4f}) for {sg1}"
                if best_R_alt is not None:
                    atoms2 = _reexpand_atoms_with_rotation(atoms2, best_R_alt, sg1)
                    msg += f"  +altindex"
                if best_k != 0:
                    n_ops = 1 + max(_op_idx(a['uid']) for a in atoms2)
                    atoms2 = _relabel_ops(atoms2, best_k, n_ops)
                    msg += f"  +symop_idx={best_k}"
                print(msg, flush=True)
                atoms1, cell1, sg1 = expand_to_p1(pdb1_path,
                                                   altloc_filter=altloc_filter,
                                                   altloc_fallback=altloc_fallback,
                                                   origin_shift=shift)
                fitme, ca_mask, uids, bfacs = match_atoms(atoms1, atoms2)
                if verbose:
                    print(f"  after origin fix: {ca_mask.sum()} CA pairs", flush=True)

    # Atom selection — restrict the fitting population BEFORE outlier MAD.
    # Default 'all' = every matched atom (backbone + side chain + waters +
    # ligands).  'backbone' = main-chain {N, CA, C, O} only — drops the
    # rotamer-noisy side chains and any non-protein content.  'ca' = CAs only.
    if atom_sel == 'all':
        sel_mask = np.ones(len(fitme), dtype=bool)
    elif atom_sel == 'backbone':
        sel_mask = np.array([u.split('_')[-2] in _BACKBONE for u in uids])
    elif atom_sel == 'ca':
        sel_mask = ca_mask.copy()
    else:
        raise ValueError(f"atom_sel must be 'all'|'backbone'|'ca', got {atom_sel!r}")
    if not sel_mask.all():
        if verbose:
            print(f"atom_sel={atom_sel!r}: {int(sel_mask.sum())}/{len(fitme)} "
                  f"atoms retained", flush=True)
        fitme   = fitme[sel_mask]
        ca_mask = ca_mask[sel_mask]
        uids    = [u for u, k in zip(uids, sel_mask) if k]
        bfacs   = bfacs[sel_mask]

    # Permanent initial outlier rejection (B-factor + shift MAD)
    if verbose:
        print("rejecting outliers...", end=' ', flush=True)
    fitme, ca_mask, uids, bfacs, _ = reject_outliers(
        fitme, ca_mask, uids, cell1, mad_sigma=outlier_sigma,
        bfacs=bfacs, b_sigma=b_sigma, verbose=verbose)
    if verbose:
        print(f"{ca_mask.sum()} CA remaining", flush=True)

    # ── Cross-validation holdout ────────────────────────────────────────
    # holdout_frac > 0: reserve a random subset of atoms that are NEVER
    # used to constrain the fit.  They stay in `fitme` and their positions
    # are still evaluated by eval_shift_field each iteration so we can
    # compare predicted vs observed shift on unseen atoms — a proper
    # error estimate that doesn't rely on the SVD covariance assuming
    # iid noise across atoms.
    holdout_mask = np.zeros(len(fitme), dtype=bool)
    if holdout_frac and holdout_frac > 0.0:
        _rng = np.random.RandomState(int(holdout_seed))
        holdout_mask = _rng.rand(len(fitme)) < float(holdout_frac)
        if verbose:
            n_hold = int(holdout_mask.sum())
            n_hold_ca = int((holdout_mask & ca_mask).sum())
            print(f"holdout: {n_hold}/{len(fitme)} atoms reserved "
                  f"({n_hold_ca} CA) — fit sees {len(fitme)-n_hold}", flush=True)
    train_mask = ~holdout_mask

    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    ca_orth = frac_to_orth(ca_shifts, cell1)
    ca_rmsd = float(np.sqrt(np.mean(np.sum(ca_orth**2, axis=1))))
    if verbose:
        print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}", flush=True)
        print(f"CA RMSD after origin fix: {ca_rmsd:.3f} Å  ({ca_mask.sum()} CA pairs)", flush=True)
    if ca_rmsd > 10.0:
        raise ValueError(f"CA RMSD after origin fix {ca_rmsd:.3f} Å > 10.0 Å. "
                         "Check structure alignment.")
    if ca_mags[max_i] > 0.35:
        raise ValueError(f"Largest fractional CA shift {ca_mags[max_i]:.4f} > 0.35 cells. "
                         "Check structure alignment.")

    # Dimension and symmetry setup
    dim_col  = {'x': 3, 'y': 4, 'z': 5, 'o': 6, 'B': 7}
    dim_list = [d for d in dimensions if d in dim_col]
    dim_indices = [dim_col[d] for d in dim_list]

    if verbose:
        print("getting symops...", end=' ', flush=True)
    proper_ops = get_proper_symops(sg1)
    do_symm = (use_symm and len(proper_ops) > 1
               and all(d in dim_list for d in 'xyz'))
    asu_mask = np.array([uid.endswith('_op0') for uid in uids])
    if verbose:
        print(f"{len(proper_ops)} proper ops, do_symm={do_symm}", flush=True)

    # Build full HKL pool down to pool_reso, sorted coarse→fine
    if verbose:
        print(f"generating HKLs to {pool_reso} Å...", end=' ', flush=True)
    all_hkls = generate_hkls(cell1, pool_reso)   # (N+1, 3), DC at index 0
    Gs     = _reciprocal_metric(cell1)
    hkl_ndc = all_hkls[1:]                        # non-DC HKLs
    hv_f    = hkl_ndc.astype(float)
    inv_d2  = np.sum((hv_f @ Gs) * hv_f, axis=1)
    d_all   = np.where(inv_d2 > 0, 1.0 / np.sqrt(inv_d2), np.inf)  # (N,)
    if verbose:
        print(f"{len(hkl_ndc)} HKLs in pool", flush=True)

    # Initial batch: all HKLs with d >= fitreso_start.  If an explicit
    # iteration_schedule was passed, ignore batch_hkls / fitreso_start /
    # fitreso_end pool sizing for the advance — use the schedule's n_used
    # values directly.  Used by fitreso_scan to land iterations exactly on
    # fr-row thresholds so one combined call produces the same per-row
    # snapshots as N independent per-fr-row calls (without redoing the
    # coarser work).
    if iteration_schedule is not None:
        sched = sorted({int(n) for n in iteration_schedule
                        if int(n) >= 2 and int(n) <= len(all_hkls)})
        if not sched:
            raise ValueError("iteration_schedule contained no usable n_used "
                             "values (need 2 ≤ n ≤ len(all_hkls))")
        n_used = sched[0]
        _sched_idx = 0
    else:
        n_initial = int(np.sum(d_all >= fitreso_start - 1e-9))
        n_initial = max(n_initial, batch_hkls)         # at least one batch
        n_initial = min(n_initial, len(hkl_ndc))
        n_used = n_initial + 1                         # +1 for DC at index 0
        sched = None

    # fit_active tracks which atoms to use; updated by residual MAD each iter
    fit_active = np.ones(len(fitme), dtype=bool)

    result_AB     = None
    result_active = None
    result_snr    = None
    result_hkls   = None
    result_rmsd   = float('inf')
    iter_count    = 0

    if verbose:
        print(f"\nProgressive fit: fitreso {fitreso_start}→{pool_reso} Å, "
              f"batch={batch_hkls}, od_margin={od_margin}, drop_snr={drop_snr}", flush=True)

    while True:
        hkls_now = all_hkls[:n_used]
        n_hkls   = len(hkls_now)

        # Symmetry reduction for current HKL set
        if do_symm:
            canon_hkls, hkl_to_canon, hkl_all_ops = assign_canonical(
                hkls_now, proper_ops)
            n_canon = len(canon_hkls)
        else:
            n_canon = n_hkls

        # Overdetermination ratio — stop before fitting if too low
        n_asu_active = int((asu_mask & fit_active).sum())
        od_ratio = (n_asu_active * 3) / max(n_canon * 6, 1)
        eff_reso  = float(d_all[n_used - 2])

        if od_ratio < od_margin:
            if verbose:
                print(f"  [iter {iter_count:2d}] fitreso={eff_reso:.2f}Å "
                      f"OD={od_ratio:.2f} < {od_margin} — stopping", flush=True)
            break

        # ── Fit ──────────────────────────────────────────────────────────────
        if do_symm:
            asu_fit_mask = asu_mask & fit_active & train_mask
            frac_asu     = fitme[asu_fit_mask, :3]
            shifts_xyz   = fitme[asu_fit_mask][:, [dim_col[d] for d in 'xyz']]
            shifts_flat  = shifts_xyz.ravel()

            obs_max_frac = None
            if bound_by_obs:
                # Per-atom fractional-shift magnitude, max across the fit
                # population.  Used as the L∞ field-magnitude prior in the
                # SVD ridge (λ = σ_noise² · 6M / obs_max²).  Recomputed per
                # iter because the active atom set changes via outlier rej.
                mags = np.sqrt(np.sum(shifts_xyz**2, axis=1))
                obs_max_frac = float(mags.max()) if len(mags) else None

            X_symm = build_design_matrix_symm(frac_asu, canon_hkls, proper_ops,
                                               hkl_to_canon, hkl_all_ops)
            params, active_c, snr_c = fit_lstsq_symm(
                X_symm, shifts_flat, drop_snr=drop_snr,
                bound_by_obs_frac=obs_max_frac, pnn_mode=pnn_mode)
            AB_xyz = expand_ab_canon(params, active_c, canon_hkls, hkls_now,
                                     hkl_to_canon, hkl_all_ops, proper_ops)

            AB = np.zeros((len(dim_list), n_hkls, 2))
            for new_i, d in enumerate(dim_list):
                if d in 'xyz':
                    AB[new_i] = AB_xyz['xyz'.index(d)]
            active = np.array([active_c[hkl_to_canon[j]] for j in range(n_hkls)])
            snr    = np.array([snr_c[hkl_to_canon[j]]    for j in range(n_hkls)])
            AB_xyz_eval = AB_xyz   # (3, n_hkls, 2)

        else:
            fit_row_mask = fit_active & train_mask
            frac_fit = fitme[fit_row_mask, :3]
            shifts   = fitme[fit_row_mask][:, dim_indices]
            X        = build_design_matrix(frac_fit, hkls_now)
            AB, active, snr = fit_lstsq(X, shifts, drop_snr=drop_snr)

            xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
            if len(xyz_dims) < 3:
                AB_xyz_eval = np.zeros((3, n_hkls, 2))
                for new_i, old_i in enumerate(xyz_dims):
                    AB_xyz_eval[new_i] = AB[old_i]
            else:
                AB_xyz_eval = AB[xyz_dims]

        # ── Residual MAD → update fit_active for next iteration ───────────
        pred       = eval_shift_field(fitme[:, :3], hkls_now, AB_xyz_eval)
        resid_frac = fitme[:, 3:6] - frac * pred
        resid_orth = frac_to_orth(resid_frac, cell2)
        resid_mag  = np.linalg.norm(resid_orth, axis=1)

        # ── CV metrics: how well does the fit predict held-out atoms? ─────
        # Cartesian RMSD of (obs - frac*pred) over train subset vs holdout
        # subset, plus the pre-fit "no-field" RMSD of holdout obs as the
        # baseline that the field is trying to beat.
        cv_train_rmsd = cv_hold_rmsd = cv_hold_prefit = None
        cv_hold_ca_rmsd = cv_hold_ca_prefit = None
        if holdout_mask.any():
            t_orth = frac_to_orth(resid_frac[train_mask], cell2)
            h_orth = frac_to_orth(resid_frac[holdout_mask], cell2)
            h_pre  = frac_to_orth(fitme[holdout_mask, 3:6], cell2)
            cv_train_rmsd  = float(np.sqrt(np.mean(np.sum(t_orth**2, axis=1))))
            cv_hold_rmsd   = float(np.sqrt(np.mean(np.sum(h_orth**2, axis=1))))
            cv_hold_prefit = float(np.sqrt(np.mean(np.sum(h_pre**2, axis=1))))
            # CA-only metrics — useful because atom_sel='all' mixes CAs with
            # side-chain / water atoms of very different shift magnitudes.
            hca = holdout_mask & ca_mask
            if hca.any():
                hca_orth = frac_to_orth(resid_frac[hca], cell2)
                hca_pre  = frac_to_orth(fitme[hca, 3:6], cell2)
                cv_hold_ca_rmsd   = float(np.sqrt(np.mean(np.sum(hca_orth**2, axis=1))))
                cv_hold_ca_prefit = float(np.sqrt(np.mean(np.sum(hca_pre**2, axis=1))))

        med_r = np.median(resid_mag)
        mad_r = np.median(np.abs(resid_mag - med_r))
        sig_r = mad_r / 0.6745
        tol_r = med_r + outlier_sigma * sig_r
        fit_active_new   = resid_mag <= tol_r
        n_dropped_resid  = int((~fit_active_new & fit_active).sum())
        fit_active       = fit_active_new

        # ── Progress report ───────────────────────────────────────────────
        rmsd     = rmsd_ca(fitme, ca_mask, hkls_now, AB_xyz_eval, cell2, frac)
        n_active = int(active.sum())
        dt       = time.time() - t0
        if verbose:
            drop_str = f"  -{n_dropped_resid} resid" if n_dropped_resid else ""
            cv_str = ""
            if cv_hold_rmsd is not None:
                cv_str = (f"  [CV] train={cv_train_rmsd:.3f} "
                          f"hold={cv_hold_rmsd:.3f}Å "
                          f"(pre={cv_hold_prefit:.3f})")
                if cv_hold_ca_rmsd is not None:
                    cv_str += (f"  holdCA={cv_hold_ca_rmsd:.3f} "
                               f"(pre={cv_hold_ca_prefit:.3f})")
            print(f"  [iter {iter_count:2d}] fitreso={eff_reso:.2f}Å  "
                  f"nhkls={n_hkls-1}  canon={n_canon}  OD={od_ratio:.2f}  "
                  f"active={n_active}  RMSD={rmsd:.3f}Å{drop_str}{cv_str}  {dt:.0f}s", flush=True)

        result_AB     = AB
        result_active = active
        result_snr    = snr
        result_hkls   = hkls_now
        result_rmsd   = rmsd
        iter_count   += 1

        if iter_callback is not None:
            iter_callback(iter_count, n_hkls - 1, n_canon, rmsd,
                          hkls_now, AB_xyz_eval, active, snr)

        if cv_callback is not None and cv_hold_rmsd is not None:
            cv_callback(iter_count, eff_reso, n_hkls - 1, n_canon,
                        rmsd, cv_train_rmsd, cv_hold_rmsd,
                        cv_hold_prefit, cv_hold_ca_rmsd,
                        cv_hold_ca_prefit, int(holdout_mask.sum()))

        if max_canon is not None and n_canon >= max_canon:
            break

        # Advance: explicit schedule (one n_used per iter) or fixed-batch step.
        if sched is not None:
            _sched_idx += 1
            if _sched_idx >= len(sched):
                break
            n_used = sched[_sched_idx]
        else:
            n_next = n_used + batch_hkls
            if n_next > len(all_hkls):
                break
            n_used = n_next

    if result_hkls is None:
        raise ValueError("No batches fitted — adjust fitreso_start/fitreso_end")

    n_ca     = int(ca_mask.sum())
    n_active = int(result_active.sum())
    dt_total = time.time() - t0
    print(f"\nProgressive fit done: {iter_count} iterations  "
          f"RMSD(CA)={result_rmsd:.3f} Å  "
          f"({n_ca} CA, {n_active} active HKLs, {dt_total:.0f} s)", flush=True)

    return BendResult(
        hkls=result_hkls, AB=result_AB, active=result_active, snr=result_snr,
        fitme=fitme, ca_mask=ca_mask, rmsd=result_rmsd,
        cell1=cell1, cell2=cell2, dimensions=''.join(dim_list),
    )


def bend_apply_pdb(pdb1_path, pdb2_path, result, frac=1.0, outpath=None):
    """Apply BendResult to pdb1, write a bent PDB aligned with pdb2."""
    nhkls = int(result.active.sum())
    if outpath is None:
        outpath = f"bent{nhkls}.pdb"

    # Build AB_xyz (3, N_hkls, 2) from result
    dim_list = list(result.dimensions)
    xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
    AB_xyz = np.zeros((3, len(result.hkls), 2))
    for new_i, old_i in enumerate(xyz_dims):
        AB_xyz[new_i] = result.AB[old_i]

    write_bent_pdb(pdb1_path, pdb2_path, result.hkls, AB_xyz, outpath, frac)
    print(f"bent PDB written to {outpath}")
    return outpath


def bend_apply_map(map_path, result, frac=1.0, outpath=None, delta=False):
    """Apply BendResult to a CCP4 map, write bent map."""
    nhkls = int(result.active.sum())
    if outpath is None:
        outpath = f"bent{nhkls}.map"

    dim_list = list(result.dimensions)
    xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
    AB_xyz = np.zeros((3, len(result.hkls), 2))
    for new_i, old_i in enumerate(xyz_dims):
        AB_xyz[new_i] = result.AB[old_i]

    make_bent_map(map_path, result.hkls, AB_xyz, outpath, cell2=result.cell2)
    if delta:
        make_delta_maps(map_path, result.hkls, AB_xyz, nhkls)
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def _parse_args(argv):
    """Parse bendfinder.com-style key=value arguments.

    Returns (pdb1, pdb2, mapfile, mtz1, mtz2, params).  When two CCP4 maps are
    given, mapfile is a tuple (mov_map, ref_map); when one is given, mapfile is
    the path string (single-fit fallback)."""
    pdb1 = pdb2 = mtz1 = mtz2 = None
    map1 = map2 = None
    params = {
        'nhkls':          30,
        'fitreso':        None,
        'drop_snr':       1.0,
        'frac':           1.0,
        'geotest':        False,
        'use_symm':       True,
        'dimensions':     'xyz',
        'fitparams':      None,
        'deltamaps':      False,
        'f_col':          None,
        'phi_col':        None,
        'run_refinement': False,
        'refine_cycles':  5,
        'fill_asu':     False,
        'fill_method':  'nan',
        'sample_rate':    0.0,
        'scan_dir':       None,
        'subtract':       'ref',
        'scan_all_fr':    False,
        'early_stop_tol': 0.01,
        'early_stop_n':   2,
    }
    for arg in argv:
        if arg.lower().endswith('.pdb'):
            if pdb1 is None:
                pdb1 = arg
            elif pdb2 is None:
                pdb2 = arg
        elif arg.lower().endswith(('.map', '.ccp4', '.mrc')):
            if map1 is None:
                map1 = arg
            elif map2 is None:
                map2 = arg
        elif arg.lower().endswith('.mtz'):
            if mtz1 is None:
                mtz1 = arg
            elif mtz2 is None:
                mtz2 = arg
        elif '=' in arg:
            key, _, val = arg.partition('=')
            key = key.strip()
            val = val.strip()
            if key in ('nhkls', 'starthkls', 'maxhkls'):
                params['nhkls'] = int(val)
            elif key == 'fitreso':
                params['fitreso'] = float(val)
            elif key == 'drop_snr':
                params['drop_snr'] = float(val)
            elif key == 'frac':
                params['frac'] = float(val)
            elif key == 'geotest':
                params['geotest'] = val.lower() in ('true', '1', 'yes')
            elif key == 'use_symm':
                params['use_symm'] = val.lower() not in ('false', '0', 'no')
            elif key == 'dimensions':
                params['dimensions'] = val
            elif key == 'fitparams':
                params['fitparams'] = val
            elif key in ('deltamaps', 'delta'):
                params['deltamaps'] = True
            elif key == 'labels':
                parts = val.split(',')
                if len(parts) == 2:
                    params['f_col']   = parts[0].strip()
                    params['phi_col'] = parts[1].strip()
            elif key in ('run_refinement', 'refine'):
                params['run_refinement'] = val.lower() not in ('false', '0', 'no')
            elif key == 'refine_cycles':
                params['refine_cycles'] = int(val)
            elif key in ('fill_asu', 'fill-asu'):
                params['fill_asu'] = val.lower() not in ('false', '0', 'no')
            elif key in ('fill_method', 'fill-method'):
                v = val.strip().lower()
                if v not in ('nan', 'sigmaa'):
                    print(f"ERROR: fill_method must be 'nan' or 'sigmaa', "
                          f"got {val!r}", file=sys.stderr)
                    sys.exit(2)
                params['fill_method'] = v
                # Choosing sigmaa implies opting into the fill itself.
                if v == 'sigmaa':
                    params['fill_asu'] = True
            elif key == 'sample_rate':
                params['sample_rate'] = float(val)
            elif key == 'scan_dir':
                params['scan_dir'] = val
            elif key == 'subtract':
                v = val.strip().lower()
                if v not in ('ref', 'bent'):
                    print(f"ERROR: subtract must be 'ref' or 'bent', got "
                          f"{val!r}", file=sys.stderr)
                    sys.exit(2)
                params['subtract'] = v
            elif key in ('scan_all_fr', 'scan-all-fr'):
                params['scan_all_fr'] = val.lower() not in ('false', '0', 'no')
            elif key in ('early_stop_tol', 'early-stop-tol'):
                params['early_stop_tol'] = float(val)
            elif key in ('early_stop_n', 'early-stop-n'):
                params['early_stop_n'] = int(val)
        elif arg in ('deltamaps', 'delta', 'nofit', 'run_refinement', 'refine',
                     'fill_asu', 'fill-asu', '--fill-asu',
                     'scan_all_fr', 'scan-all-fr', '--scan-all-fr'):
            if arg in ('deltamaps', 'delta'):
                params['deltamaps'] = True
            elif arg in ('run_refinement', 'refine'):
                params['run_refinement'] = True
            elif arg in ('fill_asu', 'fill-asu', '--fill-asu'):
                params['fill_asu'] = True
            elif arg in ('scan_all_fr', 'scan-all-fr', '--scan-all-fr'):
                params['scan_all_fr'] = True
            # nofit: ignored — lstsq always fits
    # mapfile: tuple (mov, ref) if both maps given; bare path if one; None if none
    if map2 is not None:
        mapfile = (map1, map2)
    elif map1 is not None:
        mapfile = map1
    else:
        mapfile = None
    return pdb1, pdb2, mapfile, mtz1, mtz2, params


# ══════════════════════════════════════════════════════════════════════════════
# fitreso_scan — consolidated scan: hkl00..hkl10 + fr20..fr5
# ══════════════════════════════════════════════════════════════════════════════

import math as _math
import contextlib as _contextlib
import shutil as _shutil
import tempfile as _tempfile
from scipy import ndimage as _ndimage


def _scan_cell_matrix(hdr):
    a, b, c, al, be, ga = hdr['cell']
    al, be, ga = _math.radians(al), _math.radians(be), _math.radians(ga)
    cg, sg_ = _math.cos(ga), _math.sin(ga)
    ca, cb  = _math.cos(al), _math.cos(be)
    v = _math.sqrt(1 - ca**2 - cb**2 - cg**2 + 2*ca*cb*cg)
    return np.array([
        [a,  b*cg,  c*cb],
        [0,  b*sg_, c*(ca - cb*cg)/sg_],
        [0,  0,     c*v/sg_],
    ])


def _scan_read_pdb_atoms_p1(pdb_path, src='ref', cell_override=None):
    """Read PDB, expand to P1, return list of atom dicts with frac coords.

    src           : tag stored in each atom (e.g. 'ref', 'mov') so peak
                    labels can show which structure an atom came from.
    cell_override : if not None, fractionalize Cartesian atom positions with
                    this gemmi.UnitCell instead of the PDB's own cell.  Use
                    when loading mov atoms for peak hunting against a
                    diff_map on ref's grid — keeps all fractional coords in
                    the same frame even if mov/ref cells differ slightly."""
    st       = gemmi.read_pdb(pdb_path)
    sg       = st.find_spacegroup()
    ops      = sg.operations() if sg else gemmi.SpaceGroup('P1').operations()
    frac_cell = cell_override if cell_override is not None else st.cell
    out = []
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.has_altloc() and atom.altloc != 'A':
                        continue
                    f0 = frac_cell.fractionalize(atom.pos)
                    f0 = np.array([f0.x, f0.y, f0.z])
                    for op in ops:
                        R = np.array([[op.rot[i][j] for j in range(3)]
                                      for i in range(3)], dtype=float) / 24.
                        t = np.array([op.tran[i] for i in range(3)],
                                     dtype=float) / 24.
                        out.append({'chain': chain.name,
                                    'resseq': res.seqid.num,
                                    'resname': res.name,
                                    'atomname': atom.name,
                                    'src': src,
                                    'frac': (R @ f0 + t) % 1.0})
    return out


def _scan_voxel_to_frac(idx_3d, hdr):
    nx, ny, nz = hdr['nx'], hdr['ny'], hdr['nz']
    ss, sr, sc = hdr['nsstart'], hdr['nrstart'], hdr['ncstart']
    ms, mr, mc = hdr['maps'], hdr['mapr'], hdr['mapc']
    is_, ir, ic = idx_3d
    grid  = {mc: ic,  mr: ir,  ms: is_}
    start = {mc: sc,  mr: sr,  ms: ss}
    N     = {1: nx, 2: ny, 3: nz}
    return np.array([(grid[ax] + start[ax]) / N[ax] for ax in [1, 2, 3]])


def _scan_nearest_atom(frac, atoms, M):
    arr  = np.array([a['frac'] for a in atoms])
    d    = arr - frac[np.newaxis, :]
    d   -= np.round(d)
    i    = np.argmin(np.linalg.norm(d @ M.T, axis=1))
    return atoms[i], np.linalg.norm((d @ M.T)[i])


def _scan_find_peaks(data, hdr, atoms, M, fp=5):
    """Return {'pos': (sigma, atom, dist, idx), 'neg': (sigma, atom, dist, idx)}
    where idx is the (i,j,k) voxel index of the peak in `data`.  Callers that
    used 3-tuple unpacking should be updated to ignore or use the idx."""
    sigma = data.std()
    fp3   = np.ones((fp,)*3)
    out   = {}
    _dummy = {'chain': '?', 'resseq': 0, 'resname': '?', 'atomname': '?'}
    for sign, lbl in [(1, 'pos'), (-1, 'neg')]:
        d     = sign * data
        lmax  = d == _ndimage.maximum_filter(d, footprint=fp3)
        labeled, nf = _ndimage.label(lmax & (d > 3 * sigma))
        peaks = []
        for k in range(1, nf + 1):
            mask = labeled == k
            idx  = np.unravel_index(np.argmax(d * mask), d.shape)
            peaks.append((sign * d[idx] / sigma, np.array(idx)))
        peaks.sort(key=lambda x: -abs(x[0]))
        if peaks:
            val, idx = peaks[0]
            a, dist  = _scan_nearest_atom(_scan_voxel_to_frac(idx, hdr), atoms, M)
            out[lbl] = (val, a, dist, tuple(int(i) for i in idx))
        else:
            out[lbl] = (0., _dummy, 0., None)
    return out


def _atom_label(a):
    tag = a.get('src')
    suffix = f"({tag[0]})" if tag else ''
    return f"{a['chain']}/{a['resseq']}{a['resname']}/{a['atomname']}{suffix}"


def _grid_xyz_fullcell(data, hdr):
    """Return full-cell density as (NX, NY, NZ) float64, in canonical
    crystallographic axis order (axis 0 = a, axis 1 = b, axis 2 = c).

    `data` is the (ns, nr, nc) CCP4-shaped array described by `hdr`.  If the
    header indicates an ASU subset (nonzero starts or dims < unit-cell grid),
    the map is round-tripped through `gemmi.read_ccp4_map(..., setup=True)`
    so the spacegroup operators expand the ASU to the full cell.
    """
    NX, NY, NZ = hdr['nx'], hdr['ny'], hdr['nz']
    mc, mr, ms = hdr['mapc'], hdr['mapr'], hdr['maps']
    N_by_cryst = {mc: hdr['nc'], mr: hdr['nr'], ms: hdr['ns']}
    is_fullcell = (N_by_cryst.get(1) == NX
                   and N_by_cryst.get(2) == NY
                   and N_by_cryst.get(3) == NZ
                   and hdr['ncstart'] == 0
                   and hdr['nrstart'] == 0
                   and hdr['nsstart'] == 0)
    if is_fullcell:
        # data[i_s, i_r, i_c] where i_s is the maps-axis index, etc.
        # Permute so new axis 0 = crystallographic axis 1 (a), etc.
        data_axis_for_cryst = {ms: 0, mr: 1, mc: 2}
        perm = (data_axis_for_cryst[1],
                data_axis_for_cryst[2],
                data_axis_for_cryst[3])
        return np.ascontiguousarray(data.transpose(perm).astype(np.float64))

    # ASU subset: round-trip via gemmi to apply spacegroup expansion.
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix='.map', delete=False).name
    try:
        write_ccp4(tmp, data, hdr)
        ccp4 = gemmi.read_ccp4_map(tmp)
        ccp4.setup(0.0)
        # gemmi grid.array shape is (nu, nv, nw) = (NX, NY, NZ) after setup
        return np.array(ccp4.grid.array, dtype=np.float64)
    finally:
        try: os.unlink(tmp)
        except OSError: pass


def _fft_lookup_at_hkls(grid_xyz, hkls, cell_volume):
    """Direct numpy FFT + crystallographic-convention F lookup.

    grid_xyz : (NX, NY, NZ) float density.
    hkls     : (M, 3) int Miller indices.
    cell_volume : Å³, used to match gemmi's F-vs-density normalization.

    Returns complex128 array of shape (M,) with F values matching the
    convention of `gemmi.transform_map_to_f_phi` (so amplitudes line up
    with a reference MTZ written by gemmi for the same density).
    """
    NX, NY, NZ = grid_xyz.shape
    N_grid     = NX * NY * NZ
    # numpy: A[k] = Σ a[n] exp(-2πi k·n/N); crystallographic F has +2πi.
    # ⇒ F(h,k,l) = A[(-h)%NX, (-k)%NY, (-l)%NZ]
    # Numpy F is unnormalized — multiply by V/N to match gemmi's convention.
    F_grid = np.fft.fftn(grid_xyz)
    H = hkls[:, 0].astype(int)
    K = hkls[:, 1].astype(int)
    L = hkls[:, 2].astype(int)
    F = F_grid[(-H) % NX, (-K) % NY, (-L) % NZ]
    return F * (cell_volume / N_grid)


def _fft_box_hkls(NX, NY, NZ):
    """Enumerate Friedel-unique (h,k,l) inside the FFT box of an
    (NX,NY,NZ) grid.  Excludes (0,0,0).  Used by `_map2mtz` and
    `_write_scan_mtz` to write complete MTZs without going through
    the broken `prepare_asu_data` path.
    """
    h_max = (NX - 1) // 2
    k_max = (NY - 1) // 2
    l_max = (NZ - 1) // 2
    hs = np.arange(-h_max, h_max + 1)
    ks = np.arange(-k_max, k_max + 1)
    ls = np.arange(-l_max, l_max + 1)
    H, K, L = np.meshgrid(hs, ks, ls, indexing='ij')
    H = H.ravel(); K = K.ravel(); L = L.ravel()
    # Friedel dedup: keep (l>0) OR (l==0 & k>0) OR (l==0 & k==0 & h>0)
    keep = (L > 0) | ((L == 0) & ((K > 0) | ((K == 0) & (H > 0))))
    return np.stack([H[keep], K[keep], L[keep]], axis=1).astype(np.int32)


def _canonical_hkl(h, proper_ops_R_int):
    """Return the SG-canonical representative of HKL h.

    Convention matches `assign_canonical`: lexicographic minimum of the orbit
    {R_k^T h} reduced to Friedel-unique form (first nonzero index positive).
    `proper_ops_R_int` is a list of integer 3×3 rotation matrices (the
    fractional rotation parts of the SG's proper operators).
    """
    h_arr = np.asarray(h, dtype=int)
    orbit = set()
    for R in proper_ops_R_int:
        rh = R.T @ h_arr
        t  = (int(rh[0]), int(rh[1]), int(rh[2]))
        for v in t:
            if v != 0:
                if v < 0:
                    t = (-t[0], -t[1], -t[2])
                break
        orbit.add(t)
    return min(orbit) if orbit else (int(h_arr[0]), int(h_arr[1]), int(h_arr[2]))


def _enumerate_sg_asu_hkls(cell_tuple, sg, d_min):
    """Enumerate SG-unique (h,k,l) at resolution ≥ d_min, excluding
    systematic absences and (0,0,0).

    Uses `generate_hkls` (Friedel-unique) + `_canonical_hkl` to deduplicate
    by SG orbit, then filters out HKLs forbidden by centering / glide /
    screw operators via `gemmi.GroupOps.is_systematically_absent`.
    Avoids gemmi's `ReciprocalAsu`/`prepare_asu_data` (both known to
    misclassify HKLs for certain SGs).
    """
    fu_hkls = generate_hkls(cell_tuple, d_min)
    if len(fu_hkls) and tuple(fu_hkls[0]) == (0, 0, 0):
        fu_hkls = fu_hkls[1:]
    proper_ops = get_proper_symops(sg.xhm())
    R_ints = [np.round(R).astype(int) for R, _ in proper_ops] or [np.eye(3, dtype=int)]
    ops = sg.operations()
    canon = set()
    for h in fu_hkls:
        ch = _canonical_hkl(tuple(int(x) for x in h), R_ints)
        if ops.is_systematically_absent(ch):
            continue
        canon.add(ch)
    return np.array(sorted(canon), dtype=np.int32)


def _mtz_completeness(mtz_path):
    """Return (n_obs, n_expected, fraction) for an MTZ.

    n_expected = unique SG-ASU reflections at the input's resolution limit,
    excluding systematic absences and (0,0,0).  Counts only rows with a
    finite observation in the input's amplitude or intensity column
    (F preferred: FP / F / FOBS; falls back to I: IMEAN / I / IOBS).
    Rows with NaN are treated as missing even when the HKL is present.

    If neither F nor I is present (e.g. a post-refmac map-coef MTZ with
    only FWT/PHWT), return (0, 0, 1.0) — refmac will reject it anyway,
    and the completeness gate shouldn't preempt that with a confusing
    error.
    """
    mtz = gemmi.read_mtz_file(mtz_path)
    data = np.array(mtz.array)
    lbls = mtz.column_labels()
    obs_lbl = next((l for l in ('FP', 'F', 'FOBS', 'IMEAN', 'I', 'IOBS')
                    if l in lbls), None)
    if obs_lbl is None:
        return 0, 0, 1.0
    hkls    = data[:, :3].astype(int)
    obs_col = data[:, lbls.index(obs_lbl)]
    proper_ops = get_proper_symops(mtz.spacegroup.xhm())
    R_ints = [np.round(R).astype(int) for R, _ in proper_ops] or [np.eye(3, dtype=int)]
    in_canon_with_obs = {_canonical_hkl(tuple(int(x) for x in h), R_ints)
                         for h, v in zip(hkls, obs_col) if np.isfinite(v)}
    cell = mtz.cell
    d_min = float(mtz.resolution_high())
    asu = _enumerate_sg_asu_hkls(
        (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma),
        mtz.spacegroup, d_min)
    n_expected = len(asu)
    n_obs = sum(1 for h in asu
                if tuple(int(x) for x in h) in in_canon_with_obs)
    frac = n_obs / n_expected if n_expected else 0.0
    return n_obs, n_expected, frac


def _pick_free_column(mtz):
    """Pick the first sensible free-R-flag column in `mtz`.

    A sensible column has at least one row with value 0 (the free set) AND
    at least one with a non-zero positive integer (the working set).  This
    rejects deposit MTZs that carry a degenerate FREE column (e.g. dhfr
    1rx2.mtz has FREE = all-1s in its observed rows, with the real
    cross-validation set in a separate FreeR_flag column) — picking that
    column would make refmac see 'no free reflections' and collapse SigmaA.

    Falls back to the first existing column by name preference if none of
    the candidates is sensible — caller can still proceed; refmac will
    warn but not crash.
    """
    candidates = [c.label for c in mtz.columns
                  if c.label in ('FREE', 'RFREE', 'FreeR_flag')]
    if not candidates:
        return None
    data = mtz.array
    lbls = mtz.column_labels()
    for lbl in candidates:
        col = data[:, lbls.index(lbl)]
        finite = col[np.isfinite(col)]
        if finite.size == 0:
            continue
        ints = finite.astype(int)
        if (ints == 0).any() and (ints > 0).any():
            return lbl
    return candidates[0]


def _pad_mtz_to_full_asu(in_mtz_path, out_mtz_path, d_min=None):
    """Add HKL+work-set-FreeR rows for HKLs missing from the SG ASU.

    Used by `run_refinement` to give refmac a complete-ASU input so its
    FWT/PHWT/FC_ALL output covers every SG-ASU HKL — refmac's output
    rowcount equals its input rowcount, so any HKL missing from the
    input becomes a hole in the downstream bent map.

    Added rows carry HKL + a work-set FreeR_flag and leave F/SIGF/I/SIGI
    NaN.  Refmac and phenix.refine skip NaN-Fobs rows in their LSQ
    target and R sum but still compute FC/FWT/2FOFCWT from the model for
    those HKLs — on their own internal scale, so no Fobs scale mismatch
    is possible.  Same behaviour they already apply to natively-empty
    rows in some deposits (e.g. lyso 3aw7 ships with 683 such rows).

    If the input has no F or I column, returns the input path unchanged.
    """
    mtz = gemmi.read_mtz_file(in_mtz_path)
    cell = mtz.cell
    sg   = mtz.spacegroup
    data = np.array(mtz.array, copy=True)
    lbls = mtz.column_labels()
    types = [c.type for c in mtz.columns]
    in_hkls = data[:, :3].astype(int)
    if d_min is None:
        d_min = float(mtz.resolution_high())

    # Identify amplitude / intensity / sigma / free columns.  An MTZ may
    # carry F (post-truncate / refined output) or I (raw deposit straight
    # from cif2mtz) or both — we fill whichever obs column refmac will use.
    F_lbl    = next((l for l in ('FP', 'F', 'FOBS') if l in lbls), None)
    SIG_lbl  = next((l for l in ('SIGFP', 'SIGF', 'SIGFOBS') if l in lbls), None)
    I_lbl    = next((l for l in ('IMEAN', 'I', 'IOBS') if l in lbls), None)
    SIGI_lbl = next((l for l in ('SIGIMEAN', 'SIGI', 'SIGIOBS') if l in lbls), None)
    FREE_lbl = _pick_free_column(mtz)
    if F_lbl is None and I_lbl is None:
        # Nothing to fill — no obs column refmac would consume.
        import shutil as _sh
        _sh.copyfile(in_mtz_path, out_mtz_path)
        return out_mtz_path

    # Strip systematic absences from input rows.  They're zero by symmetry;
    # presence is either a bug in the source MTZ or stale rows from a tool
    # that wrote them in.  Refmac and other CCP4 tools either silently
    # ignore them or flag them as soft errors — cleaner to drop them here.
    ops = sg.operations()
    keep_mask = np.array([not ops.is_systematically_absent(tuple(int(x) for x in h))
                          for h in in_hkls])
    n_stripped = int((~keep_mask).sum())
    if n_stripped:
        data = data[keep_mask]
        in_hkls = data[:, :3].astype(int)

    proper_ops = get_proper_symops(sg.xhm())
    R_ints = [np.round(R).astype(int) for R, _ in proper_ops] or [np.eye(3, dtype=int)]

    # Map existing input HKLs to canonical form so we don't double-count
    in_canon = {_canonical_hkl(tuple(int(x) for x in h), R_ints)
                for h in in_hkls}
    cell_tuple = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    asu_hkls = _enumerate_sg_asu_hkls(cell_tuple, sg, d_min)
    missing = [tuple(int(x) for x in h) for h in asu_hkls
               if tuple(int(x) for x in h) not in in_canon]
    if not missing:
        # All SG-ASU HKLs already present.  If we stripped any absences above,
        # write out the cleaned MTZ — otherwise just copy the input.
        if n_stripped:
            mtz.set_data(data.astype(np.float32))
            mtz.write_to_file(out_mtz_path)
        else:
            import shutil as _sh
            _sh.copyfile(in_mtz_path, out_mtz_path)
        return out_mtz_path

    missing_hkl = np.array(missing, dtype=np.int32)

    # Filled rows go in the working set (free-flag > 0), not the free set.
    # Missing HKLs cluster at high resolution and assigning them as free
    # reflections would dump a non-random shell into the cross-validation
    # set.  Use the most common positive value in the existing column,
    # falling back to 1.
    if FREE_lbl:
        existing = data[:, lbls.index(FREE_lbl)]
        existing_int = existing[np.isfinite(existing)].astype(int)
        work_vals = existing_int[existing_int > 0]
        if work_vals.size:
            vals, counts = np.unique(work_vals, return_counts=True)
            free_default = np.float32(int(vals[counts.argmax()]))
        else:
            free_default = np.float32(1)

    new_rows = np.full((len(missing_hkl), data.shape[1]),
                       np.nan, dtype=np.float32)
    new_rows[:, :3] = missing_hkl
    if FREE_lbl:
        new_rows[:, lbls.index(FREE_lbl)] = free_default

    out_data = np.vstack([data, new_rows]).astype(np.float32)
    mtz.set_data(out_data)
    mtz.write_to_file(out_mtz_path)
    return out_mtz_path


# ─── sigmaA-weighted Fc fill (SIGMAA_FC_FILL.md) ──────────────────────────────
# Principled version of `_pad_mtz_to_full_asu`.  The NaN-fill path leaves FP
# empty and relies on refmac to regenerate FWT for missing rows from its
# model.  This path instead:
#
#   1. Recomputes Fc via `gemmi sfcalc --scale-to` — bulk solvent + aniso
#      scaling puts Fc on the Fo amplitude scale (fixes the DHFR ~1/18 raw
#      DensityCalculatorX bug documented in CLAUDE.md).
#   2. Estimates per-shell σA / D / SigN / SigP on FREE==0 (held-out).
#   3. Fits ln σA (Luzzati) and ln SigN (Wilson) linearly in s²; extrapolates
#      past d_min so extrapolated σA → 0 and D·Fc → 0 (self-limiting fill).
#   4. Adds rows for every missing SG-ASU HKL with FP = D·|Fc|,
#      SIGFP = √((1-σA²)·SigN).  Present rows leave Fo/SIGFP untouched.
#   5. Emits diagnostic FC / PHIC columns on every SG-ASU HKL (scaled Fc + phi
#      from sfcalc) so downstream can build maps directly.
#
# The blend F_best = (Fo/σo² + D·Fc/σD²) / (1/σo² + 1/σD²) is deferred — refmac
# ingests the filled MTZ and produces its own FWT/PHWT downstream.

def _relabel_mtz_fp_to_f(in_mtz_path, out_mtz_path):
    """Copy MTZ with FP/SIGFP (or FOBS/SIGFOBS) relabeled to F/SIGF.

    `gemmi sfcalc --scale-to` defaults to reading F/SIGF; deposited MTZs
    typically carry FP/SIGFP.  This helper is a no-op copy when F/SIGF
    already exist.
    """
    import shutil as _sh
    _sh.copyfile(in_mtz_path, out_mtz_path)
    mtz = gemmi.read_mtz_file(out_mtz_path)
    lbls = {c.label for c in mtz.columns}
    if 'F' in lbls and 'SIGF' in lbls:
        return out_mtz_path
    for c in mtz.columns:
        if c.label in ('FP', 'FOBS') and 'F' not in lbls:
            c.label = 'F'
        elif c.label in ('SIGFP', 'SIGFOBS') and 'SIGF' not in lbls:
            c.label = 'SIGF'
    mtz.write_to_file(out_mtz_path)
    return out_mtz_path


def _recompute_fc_scaled(pdb_path, ref_mtz_path, out_fc_mtz, d_min=None,
                        verbose=True):
    """Recompute Fc with bulk solvent + anisotropic scaling to ref_mtz's F.

    Wraps `gemmi sfcalc --dmin=<d> --scale-to=<mtz>:F --to-mtz=<out>`.
    Result MTZ has columns FC / PHIC on the ref-Fo amplitude scale, so
    D·Fc is directly usable to fill missing rows without introducing a
    second Fobs population (the ~1/18-scale bug that killed the raw
    DensityCalculatorX approach on DHFR 1rx1/1rx2).

    ref_mtz_path is relabeled to F/SIGF into a scratch file first if it
    carries FP/SIGFP or FOBS/SIGFOBS — `--scale-to` doesn't accept
    arbitrary labels through its default parser.
    """
    import shutil as _sh
    gemmi_exe = (_sh.which('gemmi')
                 or os.path.join(os.environ.get('CCP4', ''), 'bin', 'gemmi'))
    if not (gemmi_exe and os.path.exists(gemmi_exe)):
        raise RuntimeError('gemmi executable not found on PATH or $CCP4/bin')

    mtz = gemmi.read_mtz_file(ref_mtz_path)
    cols = {c.label for c in mtz.columns}
    if not any(l in cols for l in ('F', 'FP', 'FOBS')):
        raise RuntimeError(f'{ref_mtz_path}: no F/FP/FOBS amplitude column')
    if d_min is None:
        d_min = float(mtz.resolution_high())

    # Relabel to F/SIGF (gemmi sfcalc --scale-to default labels) if needed.
    scratch = out_fc_mtz + '.scale_input.mtz'
    _relabel_mtz_fp_to_f(ref_mtz_path, scratch)

    cmd = [gemmi_exe, 'sfcalc', f'--dmin={d_min:.4f}',
           f'--scale-to={scratch}:F',
           f'--to-mtz={out_fc_mtz}', pdb_path]
    if verbose:
        print(f'  gemmi sfcalc --scale-to → {os.path.basename(out_fc_mtz)}',
              flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    log_path = out_fc_mtz + '.sfcalc.log'
    with open(log_path, 'w') as _lf:
        _lf.write(f'# cmd: {" ".join(cmd)}\n')
        _lf.write(f'# returncode: {r.returncode}\n\n--- stdout ---\n')
        _lf.write(r.stdout or '')
        _lf.write('\n--- stderr ---\n')
        _lf.write(r.stderr or '')
    try:
        os.remove(scratch)
    except OSError:
        pass
    if r.returncode != 0 or not os.path.exists(out_fc_mtz):
        raise RuntimeError(
            f'gemmi sfcalc failed (rc={r.returncode}); see {log_path}')
    return out_fc_mtz


def _sigmaa_shell_stats(fo, fc, s2, free_flag, n_shells=20, use_free=True,
                        min_per_shell=10):
    """Per-shell σA, D, SigN, SigP estimators on FREE==0 (held-out).

    Estimators (Read 1986, per-shell, uncentered):
        σA    = Σ(Fo·Fc) / √(ΣFo² · ΣFc²)          (amplitude correlation)
        D     = Σ(Fo·Fc) / ΣFc²                    (LS slope; the fill factor)
        SigN  = <Fo²>                              (Wilson variance of Fo)
        SigP  = <Fc²>                              (variance of Fc)

    Only finite-Fo, finite-Fc, Fc>0 rows enter.  When `use_free` and a FreeR
    column is provided, restricts to `free_flag == 0` (the held-out set) —
    working-set σA is overfit-optimistic (8sf1: 0.958 free vs 0.971 working).
    Falls back to all working rows if the free subset is too small.

    Returns dict with s²_mid (Å⁻²), per-shell σA/D/SigN/SigP arrays (NaN in
    under-populated shells), and the number of reflections used.
    """
    finite = np.isfinite(fo) & np.isfinite(fc) & (fc > 0)
    mask = finite.copy()
    n_used_free = 0
    if use_free and free_flag is not None:
        free_mask = mask & (free_flag == 0)
        if free_mask.sum() >= 4 * n_shells * min_per_shell:
            mask = free_mask
            n_used_free = int(free_mask.sum())
    fo_ = fo[mask]; fc_ = fc[mask]; s2_ = s2[mask]

    edges = np.linspace(s2_.min(), s2_.max(), n_shells + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    sigmaA = np.full(n_shells, np.nan)
    D      = np.full(n_shells, np.nan)
    SigN   = np.full(n_shells, np.nan)
    SigP   = np.full(n_shells, np.nan)
    n_per  = np.zeros(n_shells, dtype=int)
    for i in range(n_shells):
        m = (s2_ >= edges[i]) & (s2_ < edges[i+1])
        if i == n_shells - 1:
            m = (s2_ >= edges[i]) & (s2_ <= edges[i+1])
        n_per[i] = int(m.sum())
        if n_per[i] < min_per_shell:
            continue
        Fo_i = fo_[m]; Fc_i = fc_[m]
        S_ffc  = float(np.sum(Fo_i * Fc_i))
        S_fofo = float(np.sum(Fo_i ** 2))
        S_fcfc = float(np.sum(Fc_i ** 2))
        if S_fofo > 0 and S_fcfc > 0:
            sigmaA[i] = S_ffc / np.sqrt(S_fofo * S_fcfc)
            D[i]      = S_ffc / S_fcfc
            SigN[i]   = S_fofo / n_per[i]
            SigP[i]   = S_fcfc / n_per[i]
    return dict(edges=edges, s2_mid=mid, sigmaA=sigmaA, D=D,
                SigN=SigN, SigP=SigP, n_per=n_per,
                n_used_free=n_used_free, n_used=int(mask.sum()))


def _fit_sigmaa_wilson(stats, luzzati_d_max=5.0):
    """Fit ln σA (Luzzati) and ln SigN (Wilson) linearly in s².

    Both falloffs are ~Gaussian in real space → linear in s².  Fit on shells
    with valid stats; skip d > `luzzati_d_max` Å for Luzzati (bulk-solvent
    dominated at low res per Read).  Returns (sigmaA_fit, wilson_fit) where
    each is (intercept, slope) in {ln σA = intercept + slope·s²} form; slope
    is negative for Gaussian falloff.  None entries when the fit can't be
    made (too few valid shells).
    """
    ok = np.isfinite(stats['sigmaA']) & (stats['sigmaA'] > 0)
    s2 = stats['s2_mid']

    sigmaA_fit = None
    ok_luz = ok & (s2 > 1.0 / (luzzati_d_max ** 2))
    if ok_luz.sum() >= 3:
        x = s2[ok_luz]
        y = np.log(stats['sigmaA'][ok_luz])
        A = np.vstack([np.ones_like(x), x]).T
        coefs, *_ = np.linalg.lstsq(A, y, rcond=None)
        sigmaA_fit = (float(coefs[0]), float(coefs[1]))

    wilson_fit = None
    ok_w = ok & np.isfinite(stats['SigN']) & (stats['SigN'] > 0)
    if ok_w.sum() >= 3:
        x = s2[ok_w]
        y = np.log(stats['SigN'][ok_w])
        A = np.vstack([np.ones_like(x), x]).T
        coefs, *_ = np.linalg.lstsq(A, y, rcond=None)
        wilson_fit = (float(coefs[0]), float(coefs[1]))

    return sigmaA_fit, wilson_fit


def _eval_sigmaa_wilson(sigmaA_fit, wilson_fit, s2, sigmaA_stats=None,
                       clip_lo=1e-3, clip_hi=1.0):
    """Evaluate σA(s²) and SigN(s²) at arbitrary s² values.

    Uses the fitted Luzzati/Wilson lines; σA clipped to [clip_lo, clip_hi].
    D = σA · √(SigN / SigP) — SigP is known everywhere via Fc, so we return
    (σA, SigN) here and let the caller multiply by √(SigN/SigP) row-wise.

    If `sigmaA_stats` is provided and s² lies inside its measured range,
    linearly interpolate the empirical per-shell σA/SigN instead of the fit
    — cheaper than extrapolating.  Extrapolation always uses the fit.
    """
    if sigmaA_fit is None or wilson_fit is None:
        raise RuntimeError('sigmaA / Wilson fit missing — cannot evaluate')
    a_sa, b_sa = sigmaA_fit
    a_w,  b_w  = wilson_fit
    sigmaA = np.exp(a_sa + b_sa * s2)
    SigN   = np.exp(a_w  + b_w  * s2)
    if sigmaA_stats is not None:
        smid = sigmaA_stats['s2_mid']
        ok = (np.isfinite(sigmaA_stats['sigmaA'])
              & (sigmaA_stats['sigmaA'] > 0))
        if ok.sum() >= 2:
            x = smid[ok]
            in_range = (s2 >= x.min()) & (s2 <= x.max())
            sigmaA_int = np.interp(s2, x, sigmaA_stats['sigmaA'][ok])
            SigN_int   = np.interp(s2, x, sigmaA_stats['SigN'][ok])
            sigmaA = np.where(in_range, sigmaA_int, sigmaA)
            SigN   = np.where(in_range, SigN_int,   SigN)
    sigmaA = np.clip(sigmaA, clip_lo, clip_hi)
    SigN = np.maximum(SigN, 0.0)
    return sigmaA, SigN


def _sigmaa_fill_mtz(in_mtz_path, pdb_path, out_mtz_path, d_min=None,
                    n_shells=20, luzzati_d_max=5.0, use_free=True,
                    verbose=True, workdir=None):
    """Fill missing SG-ASU rows with sigmaA-weighted D·|Fc| (Fo-scaled).

    See SIGMAA_FC_FILL.md for method + rationale.  Pipeline:
      1. `gemmi sfcalc --scale-to` → Fo-scaled Fc/PHIC for all SG-ASU HKLs.
      2. Match input Fo ↔ Fc per canonical HKL; per-shell σA/D/SigN/SigP
         on FREE==0 (held-out; falls back to working set if free is empty
         or degenerate).
      3. Fit ln σA (Luzzati, skip d > `luzzati_d_max` Å) + ln SigN (Wilson)
         linearly in s².  Extrapolate past input d_min → σA → 0 → self-
         limiting fill.
      4. Missing HKLs get FP = D·|Fc|, SIGFP = √((1-σA²)·SigN),
         FreeR = work-set default; other columns NaN.  Present rows leave
         Fo/SIGFP untouched.
      5. Add diagnostic columns FC (type F) + PHIC (type P) on every row,
         holding the scaled Fc / phi from sfcalc.

    Systematic absences are stripped from the input rows (they're zero by
    symmetry).  Same convention as `_pad_mtz_to_full_asu`.

    Requires:
      - `gemmi` executable (for sfcalc --scale-to).
      - Input MTZ with F/FP/FOBS amplitude column.
      - Cell-consistent pdb_path (same cell as input MTZ; caller stretches
        beforehand if needed).
    """
    if workdir is None:
        workdir = os.path.dirname(os.path.abspath(out_mtz_path)) or '.'
    os.makedirs(workdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(out_mtz_path))[0]

    mtz_in = gemmi.read_mtz_file(in_mtz_path)
    cell = mtz_in.cell
    sg   = mtz_in.spacegroup
    data = np.array(mtz_in.array, copy=True)
    lbls = mtz_in.column_labels()
    types = [c.type for c in mtz_in.columns]
    if d_min is None:
        d_min = float(mtz_in.resolution_high())

    F_lbl    = next((l for l in ('FP', 'F', 'FOBS')          if l in lbls), None)
    SIG_lbl  = next((l for l in ('SIGFP', 'SIGF', 'SIGFOBS') if l in lbls), None)
    FREE_lbl = _pick_free_column(mtz_in)
    if F_lbl is None:
        raise RuntimeError(f'{in_mtz_path}: no F/FP/FOBS column — cannot fill')

    # ── Recompute Fo-scaled Fc via gemmi sfcalc --scale-to.
    fc_mtz_path = os.path.join(workdir, f'{stem}_fc_scaled.mtz')
    _recompute_fc_scaled(pdb_path, in_mtz_path, fc_mtz_path,
                        d_min=d_min, verbose=verbose)
    mtz_fc = gemmi.read_mtz_file(fc_mtz_path)
    fc_data = np.array(mtz_fc.array)
    fc_lbls = mtz_fc.column_labels()
    if 'FC' not in fc_lbls or 'PHIC' not in fc_lbls:
        raise RuntimeError(f'{fc_mtz_path}: expected FC/PHIC columns from '
                           f'sfcalc, got {sorted(fc_lbls)}')
    fc_hkls  = fc_data[:, :3].astype(int)
    fc_F     = fc_data[:, fc_lbls.index('FC')]
    fc_PHI   = fc_data[:, fc_lbls.index('PHIC')]

    # ── Canonicalize both HKL sets for matching.
    ops = sg.operations()
    R_ints = [np.round(R).astype(int)
              for R, _ in get_proper_symops(sg.xhm())] or [np.eye(3, dtype=int)]

    # Strip systematic absences from input.
    in_hkls_raw = data[:, :3].astype(int)
    keep_mask = np.array([not ops.is_systematically_absent(tuple(int(x) for x in h))
                          for h in in_hkls_raw])
    n_stripped = int((~keep_mask).sum())
    if n_stripped:
        data = data[keep_mask]
    in_hkls = data[:, :3].astype(int)
    in_canon = np.array([_canonical_hkl(tuple(int(x) for x in h), R_ints)
                         for h in in_hkls])
    in_canon_set = {tuple(int(x) for x in h) for h in in_canon}

    fc_canon = np.array([_canonical_hkl(tuple(int(x) for x in h), R_ints)
                         for h in fc_hkls])
    fc_canon_to_idx = {tuple(int(x) for x in h): i
                       for i, h in enumerate(fc_canon)}

    # ── Match Fc onto input row order.  Any input row whose HKL doesn't
    # appear in the sfcalc output (should not happen at same d_min, but
    # guard anyway) gets Fc/PHIC = NaN and drops out of shell stats.
    n_in = data.shape[0]
    Fc_at_in  = np.full(n_in, np.nan)
    PHI_at_in = np.full(n_in, np.nan)
    for i, h in enumerate(in_canon):
        idx = fc_canon_to_idx.get(tuple(int(x) for x in h))
        if idx is not None:
            Fc_at_in[i]  = fc_F[idx]
            PHI_at_in[i] = fc_PHI[idx]

    # ── Per-shell σA/D/SigN/SigP on FREE==0.  s² = 1/d² for each row.
    Fo   = data[:, lbls.index(F_lbl)]
    s2_in = np.array([1.0 / cell.calculate_d(tuple(int(x) for x in h)) ** 2
                      for h in in_hkls])
    free_in = (data[:, lbls.index(FREE_lbl)].astype(int)
               if FREE_lbl is not None else None)
    stats = _sigmaa_shell_stats(Fo, Fc_at_in, s2_in, free_in,
                                n_shells=n_shells, use_free=use_free)
    sigmaA_fit, wilson_fit = _fit_sigmaa_wilson(stats,
                                                luzzati_d_max=luzzati_d_max)
    if sigmaA_fit is None or wilson_fit is None:
        raise RuntimeError('sigmaA/Wilson fit failed — too few valid shells.  '
                           'Ensure input MTZ has enough finite Fo across '
                           f'the resolution range (n_shells={n_shells}, '
                           f'n_finite={int(np.isfinite(Fo).sum())}).')

    if verbose:
        n_valid = int((np.isfinite(stats['sigmaA'])
                       & (stats['sigmaA'] > 0)).sum())
        src = ('FREE==0' if stats['n_used_free']
               else 'all working rows (free set empty/degenerate)')
        a_sa, b_sa = sigmaA_fit
        a_w,  b_w  = wilson_fit
        rms_coord = float(np.sqrt(max(0.0, -b_sa) / (2.0 * np.pi ** 2)))
        print(f'  sigmaA fit ({n_valid}/{n_shells} shells on {src}): '
              f'ln σA = {a_sa:+.3f} + {b_sa:+.3f}·s²  '
              f'(rms coord err ≈ {rms_coord:.2f} Å)', flush=True)
        print(f'  Wilson fit: ln SigN = {a_w:+.3f} + {b_w:+.3f}·s²',
              flush=True)

    # ── Enumerate SG-ASU HKLs at d_min; identify missing.
    cell_tuple = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    asu_hkls = _enumerate_sg_asu_hkls(cell_tuple, sg, d_min)
    missing = [tuple(int(x) for x in h) for h in asu_hkls
               if tuple(int(x) for x in h) not in in_canon_set]

    # Work-set default free flag.
    free_default = None
    if FREE_lbl:
        col = data[:, lbls.index(FREE_lbl)]
        finite = col[np.isfinite(col)].astype(int)
        work_vals = finite[finite > 0]
        if work_vals.size:
            vals, counts = np.unique(work_vals, return_counts=True)
            free_default = np.float32(int(vals[counts.argmax()]))
        else:
            free_default = np.float32(1)

    # ── Build filled rows.  Evaluate σA/SigN per HKL, look up Fc/PHIC from
    # the sfcalc output, compute FP = D·|Fc|, SIGFP = √((1-σA²)·SigN).
    if missing:
        miss_hkl = np.array(missing, dtype=np.int32)
        miss_Fc   = np.full(len(missing), np.nan, dtype=np.float32)
        miss_PHIC = np.full(len(missing), np.nan, dtype=np.float32)
        for i, h in enumerate(missing):
            idx = fc_canon_to_idx.get(h)
            if idx is not None:
                miss_Fc[i]   = fc_F[idx]
                miss_PHIC[i] = fc_PHI[idx]

        s2_miss = np.array([1.0 / cell.calculate_d(h) ** 2 for h in missing])
        sigmaA_miss, SigN_miss = _eval_sigmaa_wilson(
            sigmaA_fit, wilson_fit, s2_miss, sigmaA_stats=stats)
        # D = σA · √(SigN/SigP) with SigP = |Fc|² locally (per-HKL, not shell
        # average) — SigP is known everywhere from Fc so we use it directly
        # rather than the shell-average SigP fit.  This makes D·Fc exactly
        # σA·√(SigN)·sign(Fc), independent of the SigP shell binning.
        with np.errstate(invalid='ignore', divide='ignore'):
            D_miss = np.where(miss_Fc > 0,
                              sigmaA_miss * np.sqrt(SigN_miss) / miss_Fc,
                              0.0)
        Fp_miss  = D_miss * miss_Fc      # = σA · √(SigN) · sign(Fc)
        sigFp_miss = np.sqrt(np.maximum(1.0 - sigmaA_miss ** 2, 0.0)
                             * np.maximum(SigN_miss, 0.0))
    else:
        miss_hkl = np.zeros((0, 3), dtype=np.int32)
        Fp_miss = np.zeros(0, dtype=np.float32)
        sigFp_miss = np.zeros(0, dtype=np.float32)
        miss_Fc = np.zeros(0, dtype=np.float32)
        miss_PHIC = np.zeros(0, dtype=np.float32)

    # ── Assemble new rows in the existing column order + append 2 diagnostic
    # columns (FC as F-type, PHIC as P-type) if not already present.
    new_col_labels  = ['FC', 'PHIC']
    new_col_types   = ['F', 'P']
    add_labels = [(l, t) for l, t in zip(new_col_labels, new_col_types)
                  if l not in lbls]

    n_extra = len(add_labels)
    ncol_out = data.shape[1] + n_extra
    n_present = data.shape[0]
    n_miss    = miss_hkl.shape[0]
    out_data = np.full((n_present + n_miss, ncol_out), np.nan,
                       dtype=np.float32)

    # Copy present rows (input columns), then add present-row FC/PHIC.
    out_data[:n_present, :data.shape[1]] = data
    col_out_of = {lbl: i for i, lbl in enumerate(lbls)}
    for j, (lbl, _t) in enumerate(add_labels):
        col_out_of[lbl] = data.shape[1] + j
    out_data[:n_present, col_out_of['FC']]   = Fc_at_in.astype(np.float32)
    out_data[:n_present, col_out_of['PHIC']] = PHI_at_in.astype(np.float32)

    # Missing rows: HKL, FP, SIGFP, FreeR, FC, PHIC; everything else NaN.
    if n_miss:
        out_data[n_present:, 0:3] = miss_hkl.astype(np.float32)
        out_data[n_present:, col_out_of[F_lbl]]   = Fp_miss.astype(np.float32)
        if SIG_lbl is not None:
            out_data[n_present:, col_out_of[SIG_lbl]] = sigFp_miss.astype(np.float32)
        if FREE_lbl is not None and free_default is not None:
            out_data[n_present:, col_out_of[FREE_lbl]] = free_default
        out_data[n_present:, col_out_of['FC']]   = miss_Fc.astype(np.float32)
        out_data[n_present:, col_out_of['PHIC']] = miss_PHIC.astype(np.float32)

    # ── Build the output MTZ.  Copy the input structure exactly (all
    # column metadata, dataset, history) via a scratch clone, add the two
    # new columns, then set the extended data array.
    mtz_out = gemmi.read_mtz_file(in_mtz_path)
    # Existing FC/PHIC in the input would collide — the sfcalc-scaled values
    # take precedence.  If they already exist, we still write our fresh
    # values into them (same column indices in col_out_of).
    for lbl, t in add_labels:
        # Attach to the first dataset (index 1 for CRYST-like MTZs; gemmi's
        # `add_column` puts it at the end of the current column list).
        dataset_id = mtz_out.datasets[-1].id if mtz_out.datasets else 0
        mtz_out.add_column(lbl, t, dataset_id=dataset_id, expand_data=False)
    mtz_out.set_data(out_data)
    mtz_out.history = list(getattr(mtz_out, 'history', []) or [])
    mtz_out.history.append(
        f'BENDFINDER sigmaA-fill: added {n_miss} rows (of {len(asu_hkls)} '
        f'SG-ASU total), stripped {n_stripped} systematic absences; '
        f'σA fit ln σA = {sigmaA_fit[0]:+.3f} + {sigmaA_fit[1]:+.3f}·s², '
        f'Wilson ln SigN = {wilson_fit[0]:+.3f} + {wilson_fit[1]:+.3f}·s²')
    mtz_out.write_to_file(out_mtz_path)
    if verbose:
        print(f'  sigmaA-fill: added {n_miss} rows, stripped {n_stripped} '
              f'absences → {os.path.basename(out_mtz_path)} '
              f'({n_present + n_miss} total rows)', flush=True)
    return out_mtz_path


def _fft_box_d_min(cell_obj, NX, NY, NZ):
    """Resolution limit (Å) supported by an FFT grid (NX,NY,NZ) on `cell_obj`.

    Conservative: takes the per-axis Nyquist min, so HKLs enumerated to this
    d_min always fit inside the FFT box.
    """
    return float(min(2.0 * cell_obj.a / NX,
                     2.0 * cell_obj.b / NY,
                     2.0 * cell_obj.c / NZ))


def _dipole_score_from_diff_norm(diff_norm, ref_d, cell, sg_name,
                                  sigma_min=4.0, n_top=30, min_sep_A=2.0):
    """Continuous embossing / dipole content of a diff map.

    For each of the top |σ|-ranked SG-unique peaks in ``diff_norm``:
      1. Find the nearest local maximum in the reference density.
      2. Compute the mirror position through that ref-density peak.
      3. Sample the diff map at the mirror; a clean dipole has
         mirror value opposite in sign to the peak.

    Per-peak dipole-ness = (-v_mirror / v_peak) * min(1, d_dens/0.5).
    The gate suppresses peaks that sit on the ref-density max (d_dens
    ≈ 0, where mirror trivially = self).  Aggregate is a σ²-weighted
    mean across the top-n peaks.

    Returns a scalar in ~[-1, +1]:
      0    = no dipole content (peaks sit on density or have no partner)
      +1   = perfect dipoles for every peak
      < 0  = same-sign mirrors (rare; usually noise-only maps)

    ``diff_norm`` and ``ref_d`` are numpy (NX, NY, NZ) grids on the same
    lattice (cell + dims).  ``sg_name`` is the H-M symbol used for
    SG-orbit dedup of peaks.
    """
    from scipy.ndimage import maximum_filter, minimum_filter, map_coordinates

    if diff_norm.shape != ref_d.shape:
        return 0.0
    dims = diff_norm.shape
    # σ-normalise ref density for the local-maximum search (relative scale
    # only; the value itself doesn't enter the dipole formula, just the
    # position of the local max).
    ref_rms = float(np.sqrt(np.mean(ref_d.astype(np.float64) ** 2)))
    if ref_rms <= 0:
        return 0.0
    ref_n = ref_d / ref_rms

    # SG symops for peak dedup (fractional).  gemmi.SpaceGroup accepts
    # either H-M or gemmi's xhm() string; fall through to trivial group
    # if lookup fails.
    try:
        sg = gemmi.SpaceGroup(sg_name)
        sg_ops = list(sg.operations())
    except Exception:
        sg_ops = None

    spacings = [cell.a / dims[0], cell.b / dims[1], cell.c / dims[2]]
    fp = tuple(max(3, int(np.ceil(min_sep_A / s)) | 1) for s in spacings)

    def _find(arr, positive):
        if positive:
            filt = maximum_filter(arr, size=fp, mode='wrap')
            mask = (arr == filt) & (arr >= sigma_min)
        else:
            filt = minimum_filter(arr, size=fp, mode='wrap')
            mask = (arr == filt) & (arr <= -sigma_min)
        idxs = np.argwhere(mask)
        vals = arr[mask]
        order = np.argsort(-np.abs(vals))[:200]      # cap
        out = []
        seen = []
        for k in order:
            i, j, m = idxs[k]
            frac = np.array([i/dims[0], j/dims[1], m/dims[2]])
            is_dup = False
            if sg_ops is not None:
                for prev_frac in seen:
                    for op in sg_ops:
                        tf = op.apply_to_xyz(list(prev_frac))
                        tf = np.array([tf[0] % 1.0, tf[1] % 1.0, tf[2] % 1.0])
                        df = frac - tf
                        df -= np.round(df)
                        p = cell.orthogonalize(gemmi.Fractional(*df))
                        if float(np.hypot(np.hypot(p.x, p.y), p.z)) < 1.5:
                            is_dup = True
                            break
                    if is_dup:
                        break
            if is_dup:
                continue
            seen.append(frac)
            out.append((frac, float(vals[k])))
            if len(out) >= 60:
                break
        return out

    peaks = sorted(_find(diff_norm, True) + _find(diff_norm, False),
                   key=lambda x: -abs(x[1]))[:n_top]
    if not peaks:
        return 0.0

    def _sample(arr, frac):
        idx = np.array([frac[i] * dims[i] for i in range(3)])
        return float(map_coordinates(arr, idx.reshape(3, 1),
                                      order=3, mode='wrap')[0])

    def _nearest_local_max(seed_frac, search_radius_A=3.0):
        r_vox = [max(1, int(np.ceil(search_radius_A / s))) for s in spacings]
        seed_idx = np.array([seed_frac[i] * dims[i] for i in range(3)])
        di = np.arange(-r_vox[0], r_vox[0]+1)
        dj = np.arange(-r_vox[1], r_vox[1]+1)
        dk = np.arange(-r_vox[2], r_vox[2]+1)
        ii = (int(seed_idx[0]) + di) % dims[0]
        jj = (int(seed_idx[1]) + dj) % dims[1]
        kk = (int(seed_idx[2]) + dk) % dims[2]
        sub = ref_n[np.ix_(ii, jj, kk)]
        lm = (sub == maximum_filter(sub, size=3, mode='reflect'))
        lm_pts = np.argwhere(lm)
        if len(lm_pts) == 0:
            return seed_frac, 0.0
        best_dist = None
        best_frac = seed_frac
        for lp in lm_pts:
            i2 = int(ii[lp[0]]); j2 = int(jj[lp[1]]); k2 = int(kk[lp[2]])
            lp_frac = np.array([i2/dims[0], j2/dims[1], k2/dims[2]])
            df = lp_frac - seed_frac
            df -= np.round(df)
            p = cell.orthogonalize(gemmi.Fractional(*df))
            d = float(np.hypot(np.hypot(p.x, p.y), p.z))
            if best_dist is None or d < best_dist:
                best_dist = d
                best_frac = lp_frac
        return best_frac, best_dist

    d_num = 0.0
    d_den = 0.0
    for frac, v in peaks:
        q_frac, d_dens = _nearest_local_max(frac)
        p_prime = np.array([(2*q_frac[k] - frac[k]) % 1.0 for k in range(3)])
        v_mirror = _sample(diff_norm, p_prime)
        gate = min(1.0, d_dens / 0.5)
        dip_ness = (-v_mirror / v) * gate if v != 0 else 0.0
        w = v * v
        d_num += w * dip_ness
        d_den += w
    return float(d_num / d_den) if d_den > 0 else 0.0


def _map2mtz(mapfile, mtzfile):
    """Convert a CCP4 map to an SG-ASU MTZ via direct numpy FFT.

    Replaces the gemmi `transform_map_to_f_phi + prepare_asu_data` pipeline
    which silently drops ~50% of unique reflections for trigonal/hexagonal
    spacegroups (v0.6 wedge bug; v0.7 missing-k=0-plane bug; see CLAUDE.md).

    Writes one row per SG-unique HKL — NOT a P1 / Friedel-unique full-box
    MTZ.  Encoding the same density in N point-group copies would multiply
    the inverse-FFT'd map by N.
    """
    data, hdr = read_ccp4(mapfile)
    cell_tuple = hdr['cell']
    cell_obj   = gemmi.UnitCell(*cell_tuple)
    sg_num     = hdr.get('spacegroup')
    if sg_num is None:
        ispg = struct.unpack_from('<i', hdr['raw_header'], 88)[0]
        sg_num = ispg if ispg > 0 else 1
    sg = gemmi.find_spacegroup_by_number(int(sg_num)) or gemmi.SpaceGroup('P 1')

    grid_xyz = _grid_xyz_fullcell(data, dict(hdr, spacegroup=sg_num))
    NX, NY, NZ = grid_xyz.shape
    d_min = _fft_box_d_min(cell_obj, NX, NY, NZ)
    hkls  = _enumerate_sg_asu_hkls(cell_tuple, sg, d_min)
    F     = _fft_lookup_at_hkls(grid_xyz, hkls, cell_obj.volume)
    mtz  = gemmi.Mtz(with_base=True)
    mtz.spacegroup = sg
    mtz.cell       = cell_obj
    mtz.add_dataset('dataset')
    mtz.add_column('F',   'F')
    mtz.add_column('PHI', 'P')
    arr = np.zeros((len(hkls), 5), dtype=np.float32)
    arr[:, :3] = hkls
    arr[:,  3] = np.abs(F).astype(np.float32)
    arr[:,  4] = np.degrees(np.angle(F)).astype(np.float32)
    mtz.set_data(arr)
    mtz.write_to_file(mtzfile)


def _fspace_scale_and_diff(bent_map, ref_h, ref_mtz_path, ref_f_col, ref_phi_col,
                            mtzfile, anisotropic=False, n_cycles=4,
                            subtract='ref'):
    """FFT bent_map; F-space fit (k, B) of bent → ref MTZ; compute scaled bent
    F + F-space DELFWT; write bent.mtz; return (diff_real, riso, k, B).

    subtract : 'ref'  → diff = bent_scaled − ref   (default)
                        positive peak = density present in bent but absent (or
                        weaker) in ref
               'bent' → diff = ref − bent_scaled
                        positive peak = density present in ref but absent (or
                        weaker) in bent

    diff_real is the inverse-FFT of the F-space diff on ref's grid — i.e. the
    residual density after absorbing the resolution-dependent scale.  See
    compute_riso docstring for the LS target details.
    """
    if subtract not in ('ref', 'bent'):
        raise ValueError(f"subtract must be 'ref' or 'bent', got {subtract!r}")
    from scipy.optimize import least_squares

    cell = gemmi.UnitCell(*ref_h['cell'])
    sg   = gemmi.find_spacegroup_by_number(ref_h['spacegroup'])

    # ── Load ref MTZ first — its HKL list is the lookup target ────────────
    ref_mtz_obj = gemmi.read_mtz_file(ref_mtz_path)
    ref_data    = np.array(ref_mtz_obj.array)
    ref_lbls    = ref_mtz_obj.column_labels()
    ref_hkl     = ref_data[:, :3].astype(np.int32)
    ref_F_all   = ref_data[:, ref_lbls.index(ref_f_col)]
    ref_PH_all  = ref_data[:, ref_lbls.index(ref_phi_col)]

    # ── FFT bent_map and look up F directly at the ref HKLs.  Replaces the
    # broken gemmi `transform_map_to_f_phi + prepare_asu_data` pipeline,
    # which dropped ~50% of unique reflections for trigonal/hexagonal SGs
    # (v0.6 wedge bug; v0.7 missing-k=0-plane bug; see CLAUDE.md).  Direct
    # lookup gives 100% HKL match by construction.
    grid_xyz = _grid_xyz_fullcell(bent_map, ref_h)
    bent_c   = _fft_lookup_at_hkls(grid_xyz, ref_hkl, cell.volume)

    hkl = ref_hkl
    bF  = np.abs(bent_c)
    bPH = np.degrees(np.angle(bent_c))
    rF  = ref_F_all
    rPH = ref_PH_all

    ok = (rF > 0) & (bF > 0) & np.isfinite(rF) & np.isfinite(bF)
    hkl, rF, rPH, bF, bPH = hkl[ok], rF[ok], rPH[ok], bF[ok], bPH[ok]
    if len(rF) < 10:
        return None, None, None, None

    # ── Reciprocal-cart h, s² = 1/d²
    M_frac = np.array(cell.frac.mat.tolist())
    h_orth = hkl.astype(float) @ M_frac
    s_sq   = np.sum(h_orth ** 2, axis=1)

    # ── Fit (k, B) [or (k, V) anisotropic] via F-space LS, iterated
    use = np.ones(len(rF), dtype=bool)
    if not anisotropic:
        params = np.array([np.log(max(np.dot(rF, bF) / np.dot(bF, bF), 1e-10)), 0.0])
        def _scale_fn(p, m=None):
            ss = s_sq if m is None else s_sq[m]
            return np.exp(p[0]) * np.exp(-p[1] / 4.0 * ss)
    else:
        params = np.array([np.log(max(np.dot(rF, bF) / np.dot(bF, bF), 1e-10)),
                           0., 0., 0., 0., 0., 0.])
        def _Vmat(p):
            V11, V22, V33, V12, V13, V23 = p[1:]
            return np.array([[V11, V12, V13],
                             [V12, V22, V23],
                             [V13, V23, V33]])
        def _scale_fn(p, m=None):
            hh = h_orth if m is None else h_orth[m]
            return np.exp(p[0]) * np.exp(-np.einsum('ij,ni,nj->n',
                                                     _Vmat(p), hh, hh))

    for _ in range(n_cycles):
        try:
            res = least_squares(lambda p: rF[use] - _scale_fn(p, use) * bF[use],
                                 params, method='lm', max_nfev=200)
            params = res.x
        except Exception:
            break

    scale = _scale_fn(params)
    if not anisotropic:
        k_fit = float(np.exp(params[0])); B_fit = float(params[1])
    else:
        k_fit = float(np.exp(params[0]))
        B_fit = float(4.0 * np.trace(_Vmat(params)) / 3.0)

    bF_scaled = scale * bF
    riso = float(np.sum(np.abs(rF - bF_scaled)) / np.sum(np.abs(rF)))

    # ── F-space diff (sign controlled by `subtract`)
    ref_c       = rF        * np.exp(1j * np.radians(rPH))
    bent_c_scl  = bF_scaled * np.exp(1j * np.radians(bPH))
    if subtract == 'ref':
        diff_c = bent_c_scl - ref_c           # positive peak = bent has more
    else:
        diff_c = ref_c - bent_c_scl           # positive peak = ref has more
    diff_F      = np.abs(diff_c).astype(np.float32)
    diff_PH     = np.degrees(np.angle(diff_c)).astype(np.float32)

    # ── Write bent.mtz with scaled FDM and F-space DELFWT
    out = gemmi.Mtz(with_base=True)
    out.cell, out.spacegroup = cell, sg
    out.add_dataset('dataset')
    out.add_column('FDM',     'F')
    out.add_column('PHIDM',   'P')
    out.add_column('DELFWT',  'F')
    out.add_column('PHDELWT', 'P')
    arr = np.zeros((len(hkl), 7), dtype=np.float32)
    arr[:, :3] = hkl
    arr[:,  3] = bF_scaled.astype(np.float32)
    arr[:,  4] = bPH.astype(np.float32)
    arr[:,  5] = diff_F
    arr[:,  6] = diff_PH
    out.set_data(arr)
    out.write_to_file(mtzfile)

    # ── Inverse-FFT diff → real-space at ref grid
    diff_mtz = gemmi.Mtz(with_base=True)
    diff_mtz.cell, diff_mtz.spacegroup = cell, sg
    diff_mtz.add_dataset('d')
    diff_mtz.add_column('F',   'F')
    diff_mtz.add_column('PHI', 'P')
    darr = np.zeros((len(hkl), 5), dtype=np.float32)
    darr[:, :3] = hkl
    darr[:,  3] = diff_F
    darr[:,  4] = diff_PH
    diff_mtz.set_data(darr)
    # Sample at the same grid as bent_map (ref_h)
    nu, nv, nw = ref_h['nx'], ref_h['ny'], ref_h['nz']
    grid = diff_mtz.transform_f_phi_to_map('F', 'PHI', exact_size=(nu, nv, nw))
    full_diff = np.array(grid, copy=True)
    # full_diff is (nu, nv, nw) = (NX, NY, NZ); rearrange to ref_h's
    # (ns, nr, nc) = (maps-axis, mapr-axis, mapc-axis) and slice to ASU
    mapc, mapr, maps = ref_h['mapc'], ref_h['mapr'], ref_h['maps']
    diff_full = np.transpose(full_diff, (maps - 1, mapr - 1, mapc - 1))
    Ng = {1: nu, 2: nv, 3: nw}
    ncs, nrs, nss = ref_h['ncstart'], ref_h['nrstart'], ref_h['nsstart']
    nc,  nr,  ns  = ref_h['nc'],      ref_h['nr'],      ref_h['ns']
    # Slice ASU window (ncstart..ncstart+nc-1 etc, wrapping)
    def _slc(start, n, N):
        idx = (np.arange(n) + start) % N
        return idx
    iss = _slc(nss, ns, Ng[maps])
    irr = _slc(nrs, nr, Ng[mapr])
    icc = _slc(ncs, nc, Ng[mapc])
    diff_real = diff_full[np.ix_(iss, irr, icc)].astype(np.float32)
    return diff_real, riso, k_fit, B_fit


def _write_scan_mtz(bent_arr, diff_arr, ref_h, mtzfile):
    """Write bent.mtz with FDM/PHIDM (bent map) + DELFWT/PHDELWT (diff map).

    SG-ASU MTZ via direct numpy FFT (not P1 / Friedel-box — see
    `_map2mtz` docstring for why).
    """
    cell_obj = gemmi.UnitCell(*ref_h['cell'])
    sg_num   = ref_h.get('spacegroup', 1)
    sg = gemmi.find_spacegroup_by_number(int(sg_num)) or gemmi.SpaceGroup('P 1')

    grid_b = _grid_xyz_fullcell(bent_arr, ref_h)
    grid_d = _grid_xyz_fullcell(diff_arr, ref_h)
    NX, NY, NZ = grid_b.shape
    d_min = _fft_box_d_min(cell_obj, NX, NY, NZ)
    hkls = _enumerate_sg_asu_hkls(
        (cell_obj.a, cell_obj.b, cell_obj.c,
         cell_obj.alpha, cell_obj.beta, cell_obj.gamma),
        sg, d_min)
    Fb = _fft_lookup_at_hkls(grid_b, hkls, cell_obj.volume)
    Fd = _fft_lookup_at_hkls(grid_d, hkls, cell_obj.volume)

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = sg
    mtz.cell       = cell_obj
    mtz.add_dataset('dataset')
    mtz.add_column('FDM',     'F')
    mtz.add_column('PHIDM',   'P')
    mtz.add_column('DELFWT',  'F')
    mtz.add_column('PHDELWT', 'P')
    arr = np.zeros((len(hkls), 7), dtype=np.float32)
    arr[:, :3] = hkls
    arr[:,  3] = np.abs(Fb).astype(np.float32)
    arr[:,  4] = np.degrees(np.angle(Fb)).astype(np.float32)
    arr[:,  5] = np.abs(Fd).astype(np.float32)
    arr[:,  6] = np.degrees(np.angle(Fd)).astype(np.float32)
    mtz.set_data(arr)
    mtz.write_to_file(mtzfile)


def _eval_chunked(ref_pts, hkls, AB_xyz, chunk=50000):
    delta = np.zeros_like(ref_pts)
    for i in range(0, len(ref_pts), chunk):
        delta[i:i+chunk] = eval_shift_field(ref_pts[i:i+chunk], hkls, AB_xyz)
    return delta


def fitreso_scan(
    mov_pdb, ref_pdb, mov_mtz, ref_mtz, scan_dir,
    f_col=None, phi_col=None,
    run_refinement_flag=False, refine_cycles=5,
    fill_asu=False, fill_method='nan',
    sample_rate=0.0,
    mov_fullcell=None,
    fitreso_list=(20, 15, 12, 10, 8, 7, 6, 5),
    max_hkl_scan=10,
    outlier_sigma=2.5, b_sigma=3.0, drop_snr=0.0, od_margin=1.5,
    batch_hkls=100, chunk_size=50000,
    riso_n_cycles=4, riso_sigma_cut=float('inf'),
    subtract='ref', atom_sel='all',
    bound_by_obs=True, pnn_mode=None,
    preserve_mov_numbering=True,
    scan_all_fr=False,
    early_stop_tol=0.01, early_stop_n=2,
    verbose=True,
):
    """Run the standard fitreso scan and write outputs to scan_dir.

    Section 1  — hkl00: zero shift (unregistered baseline)
    Section 2  — hkl01..hklNN: single progressive call, one canonical HKL at a time
    Section 3  — fr<N>: separate progressive calls at each resolution in fitreso_list

    Parameters
    ----------
    mov_pdb, ref_pdb : str   Moving and reference PDB paths
    mov_mtz, ref_mtz : str   Moving and reference MTZ paths (or CCP4 .map for compat)
    scan_dir         : str   Output directory (created if needed)
    f_col, phi_col   : str   MTZ column overrides (default: auto-detect FWT/2FOFCWT)
    run_refinement_flag : bool  Run refmac/phenix to generate FWT/PHWT first
    refine_cycles    : int   Refinement cycles when run_refinement_flag=True
    sample_rate      : float FFT oversampling for MTZ → map (default 0.0 = gemmi auto)
    mov_fullcell     : bool  CCP4 only: expand mov to full unit cell so wrap-mode
                              interpolation avoids ASU-boundary noise (default
                              None = auto-enable for CCP4 input)
    fitreso_list     : seq   Resolution endpoints for section 3 (Å)
    max_hkl_scan     : int   Number of non-DC canonical HKLs in section 2 (default 10)
    chunk_size       : int   Voxel batch size for eval_shift_field (default 50000)
    subtract         : str   Diff sign convention.  'ref' (default) → diff =
                              bent − ref so positive peaks = density present
                              in bent but absent (or weaker) in ref.  'bent' →
                              diff = ref − bent (opposite sign).
    """
    if subtract not in ('ref', 'bent'):
        raise ValueError(f"subtract must be 'ref' or 'bent', got {subtract!r}")
    import time

    _CANONICAL_F = {'FWT', '2FOFCWT'}

    os.makedirs(scan_dir, exist_ok=True)

    # ── normalize mov PDB chain/resnum vs ref (apply BEFORE anything else) ───
    # Some deposits use a different chain ID (e.g. mov='B' vs ref='A') or a
    # constant resnum offset (e.g. old-PDB +1000 convention).  Without this
    # both resolve_altindex's _collect_matched_ca and downstream match_atoms
    # fail to find any pairs.  The replacement PDB has identical Cartesian
    # coordinates — only the (chain, resnum) labels change — so all downstream
    # consumers (pre/, refmac, altindex, bend) see a consistent view.
    _norm_pdb = os.path.join(scan_dir, '_mov_normalized.pdb')
    orig_mov_pdb = mov_pdb                           # untouched user-supplied path
    mov_pdb, norm_plan = _normalize_mov_pdb(
        mov_pdb, ref_pdb, _norm_pdb, verbose=verbose)

    # ── snapshot the truly-raw inputs (used by the pre/ section below) ───────
    # pre/ shows the "old-fashioned" diff map: mov sampled on ref's grid with
    # NO resolve_altindex / origin transformation.  ref is fine to take post-
    # refmac (refinement settles atoms by <0.1 Å vs deposited); mov MUST stay
    # pre-resolve so the pre row honestly reports what the unaligned mov would
    # look like.  If the raw mov MTZ lacks FWT/PHWT we refmac it once below so
    # pre/ has phases to FFT, but the pdb+frame stay unchanged.
    raw_mov_pdb_in = mov_pdb
    raw_mov_mtz_in = mov_mtz

    # ── optional refinement (parallel: ref + pre-mov-raw) ────────────────────
    # Refmac can drop disordered atoms during the run, so subsequent steps
    # (altindex resolution + re-refinement, write_bent_pdb, peak hunting)
    # all use the *refined* PDB to keep the atom set consistent.
    #
    # ref refmac and pre-mov-raw refmac are fully independent — no shared
    # inputs, no shared outputs — so we dispatch both to a
    # ThreadPoolExecutor and wait for the pair.  Halves the front-end
    # wall time on the common case where BOTH need refining.  If only
    # one is needed, the executor still runs it (single-slot) and the
    # code path is identical.
    pre_mov_pdb = raw_mov_pdb_in
    pre_mov_mtz = raw_mov_mtz_in
    _need_ref_refmac = (run_refinement_flag
                        and ref_mtz.lower().endswith('.mtz'))
    _need_premov_refmac = (run_refinement_flag
                           and pre_mov_mtz.lower().endswith('.mtz')
                           and 'FWT' not in {c.label for c in
                               gemmi.read_mtz_file(pre_mov_mtz).columns})
    if _need_ref_refmac or _need_premov_refmac:
        from concurrent.futures import ThreadPoolExecutor as _TPE
        with _TPE(max_workers=2) as _refmac_pool:
            _ref_fut = _premov_fut = None
            if _need_ref_refmac:
                _ref_fut = _refmac_pool.submit(
                    run_refinement, ref_pdb, ref_mtz,
                    outdir=os.path.join(scan_dir, 'refine_ref'),
                    n_cycles=refine_cycles, fill_asu=fill_asu,
                    fill_method=fill_method)
            if _need_premov_refmac:
                # pre/ FFTs raw_mov_mtz to a map.  If FWT/PHWT aren't
                # present, refmac one pass on raw mov against its own
                # data — same idea as the post-resolve mov refmac below,
                # but starting from the un-transformed mov so pre/
                # reflects the un-aligned state.
                _premov_fut = _refmac_pool.submit(
                    run_refinement, pre_mov_pdb, pre_mov_mtz,
                    outdir=os.path.join(scan_dir, 'refine_mov_raw'),
                    n_cycles=refine_cycles, fill_asu=fill_asu,
                    fill_method=fill_method)
            if _ref_fut is not None:
                ref_pdb, ref_mtz = _ref_fut.result()
            if _premov_fut is not None:
                pre_mov_pdb, pre_mov_mtz = _premov_fut.result()

    # ── _load helper (used by both the ref/pre_mov phase below and the
    # post-altindex mov phase further down) ─────────────────────────────────
    def _load(path, tag, fullcell=False):
        ext = os.path.splitext(path)[1].lower()
        if ext == '.mtz':
            fc, ph = detect_2fofc_cols(path, f_col, phi_col)
            if verbose:
                print(f'  {tag}: {os.path.basename(path)}  columns {fc}/{ph}',
                      flush=True)
            data, hdr = mtz_to_map_data(path, fc, ph, sample_rate=sample_rate)
            return data, hdr, path, fc
        if fullcell:
            data, hdr = read_ccp4_fullcell(path)
        else:
            data, hdr = read_ccp4(path)
        return data, hdr, None, None

    # ── Altindex-independent setup: ref density + pre_mov density + grids ──
    # Everything below (up to the altindex launch) uses ONLY the pre-altindex
    # inputs (raw mov_pdb, raw mov_mtz, pre_mov_pdb, pre_mov_mtz, ref_pdb,
    # ref_mtz).  This lets us dispatch resolve_altindex to a background
    # thread and run Section 0 (pre) concurrently — the pre section needs
    # none of altindex's outputs.
    ref_d, ref_h, ref_mtz_resolved, ref_f_col = _load(ref_mtz, 'ref')

    pre_mov_fullcell = pre_mov_mtz.lower().endswith(('.map', '.ccp4', '.mrc'))
    pre_mov_d, pre_mov_h, _, pre_mov_f_col = _load(
        pre_mov_mtz, 'pre_mov', fullcell=pre_mov_fullcell)
    pre_mov_pad_mode = ('wrap' if pre_mov_mtz.lower().endswith('.mtz')
                                   or pre_mov_fullcell else 'reflect')
    # col_suffix — filename suffix for bent{col_suffix}.map / .mtz / diff etc.
    # Almost always '' (canonical FWT input); non-canonical F cols produce
    # e.g. _F suffix.  Set from pre_mov here so Section 0 (pre) has a
    # value; re-derived from mov_f_col after altindex for Section 1+
    # (usually identical because mov and pre_mov share columns pre-altindex).
    col_suffix = ('' if (pre_mov_f_col is None or pre_mov_f_col in _CANONICAL_F)
                  else f'_{pre_mov_f_col}')

    ns2, nr2, nc2 = ref_h['ns'], ref_h['nr'], ref_h['nc']
    nx2, ny2, nz2 = ref_h['nx'], ref_h['ny'], ref_h['nz']
    mc, mr, ms    = ref_h['mapc'], ref_h['mapr'], ref_h['maps']

    g_s, g_r, g_c = np.meshgrid(np.arange(ns2), np.arange(nr2), np.arange(nc2),
                                 indexing='ij')
    amap  = {mc: g_c, mr: g_r, ms: g_s}
    start = {mc: ref_h['ncstart'], mr: ref_h['nrstart'], ms: ref_h['nsstart']}
    Nxyz  = {1: nx2, 2: ny2, 3: nz2}
    frac_x = (amap[1].ravel() + start[1]) / Nxyz[1]
    frac_y = (amap[2].ravel() + start[2]) / Nxyz[2]
    frac_z = (amap[3].ravel() + start[3]) / Nxyz[3]
    ref_pts = np.stack([frac_x, frac_y, frac_z], axis=1)

    M         = _scan_cell_matrix(ref_h)
    _ref_cell = gemmi.UnitCell(*ref_h['cell'])
    # SG name used by _dipole_score_from_diff_norm and the pre-row dipole
    # calculation.  Read from ref_pdb (authoritative CRYST1 record) —
    # ref_h.get('spacegroup') gives the SG *number* which gemmi can look up.
    try:
        _sg_name = gemmi.read_pdb(ref_pdb).spacegroup_hm
    except Exception:
        _sg_name = None

    # pre_mov_display_pdb / pre_atoms — used by Section 0 (pre), don't depend
    # on altindex.  mov_display_pdb / atoms are computed AFTER altindex.
    if preserve_mov_numbering and norm_plan is not None:
        pre_mov_display_pdb = _reverse_normalize_pdb(
            pre_mov_pdb, norm_plan,
            os.path.join(scan_dir, '_pre_mov_display.pdb'))
    else:
        pre_mov_display_pdb = pre_mov_pdb

    # Pre/ uses the un-transformed mov pdb (resolve_altindex hasn't run on it).
    # For cross-cell pairs (lipox, porin) the raw mov atoms are NOT in ref's
    # frame — peak labelling for pre will pick whatever ref or raw-mov atom
    # happens to land closest in fractional space, which is the honest "what
    # would the un-aligned diff map look like?" answer.
    pre_atoms = (_scan_read_pdb_atoms_p1(ref_pdb, src='ref',
                                          cell_override=_ref_cell)
                 + _scan_read_pdb_atoms_p1(pre_mov_display_pdb, src='mov',
                                            cell_override=_ref_cell))

    # ── ref MTZ for Riso (FT of reference density map) ────────────────────────
    _PHI_LOOKUP = {'FWT': 'PHWT', '2FOFCWT': 'PH2FOFCWT', 'F': 'PHI'}
    if ref_mtz_resolved is not None:
        # ref input is already an MTZ — use it directly for Riso
        riso_ref_mtz = ref_mtz_resolved
        riso_ref_col = ref_f_col
        riso_ref_phi = _PHI_LOOKUP.get(ref_f_col, 'PHWT')
    else:
        riso_ref_mtz = os.path.join(scan_dir, 'ref_riso.mtz')
        riso_ref_col = 'F'
        riso_ref_phi = 'PHI'
        if not os.path.exists(riso_ref_mtz):
            if verbose:
                print('Converting ref map to MTZ for Riso...', flush=True)
            _map2mtz(ref_mtz, riso_ref_mtz)

    def _relsymlink_dir(src, dst):
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)
    _relsymlink_dir(os.path.abspath(ref_pdb), os.path.join(scan_dir, 'ref.pdb'))
    # ref.mtz: always a symlink — riso_ref_mtz is either the user-supplied ref MTZ
    # or the FFT'd-from-CCP4-map fallback we constructed above.
    _relsymlink_dir(os.path.abspath(riso_ref_mtz),
                     os.path.join(scan_dir, 'ref.mtz'))

    _h0  = np.zeros((0, 3), dtype=int)
    _AB0 = np.zeros((3, 0, 2))

    # ── compute pre-origin (raw) RMSD ────────────────────────────────────────
    # Use the *truly raw* mov pdb (pre-resolve_altindex) so this number really
    # is "what you'd see if you just downloaded the two structures and
    # subtracted."  For systems where resolve_altindex took action='none' this
    # equals hkl00_rmsd; for systems where it applied an origin/SG/altindex
    # op, raw_rmsd will be substantially larger.
    with open(os.devnull, 'w') as _dev, _contextlib.redirect_stdout(_dev):
        _a1, _cell1, _ = expand_to_p1(pre_mov_pdb)
        _a2, _cell2, _ = expand_to_p1(ref_pdb)
        _fm, _ca, _uid, _bf = match_atoms(_a1, _a2)
        _fm, _ca, _, _, _ = reject_outliers(_fm, _ca, _uid, _cell1,
                                            mad_sigma=outlier_sigma,
                                            b_sigma=b_sigma, bfacs=_bf)
    raw_rmsd = rmsd_ca(_fm, _ca, _h0, _AB0, _cell2)

    # ── Launch resolve_altindex in a background thread ──────────────────────
    # Runs concurrently with Section 0 (pre) below.  Section 0 depends only
    # on ref + raw/pre_mov (all set above), so we don't need altindex's
    # result until Section 1 (hkl00) starts.  Post-altindex mov refmac
    # (if needed), post-altindex mov density load, mov_display_pdb, atoms,
    # and hkl00_rmsd all happen after we wait, before Section 1.
    from concurrent.futures import ThreadPoolExecutor as _TPE2
    _altindex_pool = _TPE2(max_workers=1)
    _altindex_fut = None
    if mov_mtz.lower().endswith('.mtz'):
        _altindex_fut = _altindex_pool.submit(
            resolve_altindex, mov_pdb, ref_pdb, mov_mtz,
            outdir=os.path.join(scan_dir, 'altindex_resolve'),
            refine_cycles=refine_cycles, verbose=verbose,
            fill_asu=fill_asu, fill_method=fill_method)

    # unbent.pdb: written further down (after hkl00 is set up) as a symlink
    # to hkl00/bent.pdb, keeping it in the same cell+frame as unbent.mtz
    # and overlaying cleanly with ref structures in Coot.

    # ── print header + "pre" row ─────────────────────────────────────────────
    if verbose:
        _sign_note = ('diff = bent − ref  (+peak ⇒ density in bent absent from ref)'
                      if subtract == 'ref' else
                      'diff = ref − bent  (+peak ⇒ density in ref absent from bent)')
        print(f"\nsign convention: {_sign_note}", flush=True)
        print(f"\n{'label':>7}  {'RMSD':>6}  {'active':>6}  "
              f"{'Rbent':>6}  {'Rbend':>6}  {'bondZ/angZ':>13}  "
              f"{'peak':>7}  {'atom':>20}", flush=True)
        print('-' * 96, flush=True)

    # The PSDVF maps every fractional point in the moving cell to its
    # corresponding fractional point in the reference cell.  For PSDVF=0
    # (hkl00) the transform preserves fractional coordinates exactly —
    # bent_map(ref_frac) = mov_d(SAME ref_frac).  Cell differences
    # between mov and ref are reflected only in the cell header of the
    # output map, not in any cartesian rescaling.  Hence: sample mov_d
    # at the same fractional positions as ref_pts, no cell-aware
    # cartesian round-trip.
    def _map_query(ref_pts, delta=None):
        return ref_pts if delta is None else ref_pts - delta

    def _atoms1_pts(ref_pts):
        return ref_pts

    # Path to the hkl00 (unbent-resampled) bent.mtz — set after first _save_point
    _unbent_bent = [None]
    # Accumulate one dict per scan point for the final scan_fitreso.log table
    _log_rows = []

    def _save_point(label, outdir, bent_map, rmsd_before, rmsd_after,
                    n_active, t_elapsed, hkls=None, AB_xyz=None):
        bent_map_name    = f'bent{col_suffix}.map'
        diffnorm_name    = f'diff_norm{col_suffix}.map'
        bent_mtz_name    = f'bent{col_suffix}.mtz'
        write_ccp4(f'{outdir}/{bent_map_name}', bent_map, ref_h)
        geom = None
        if hkls is not None and AB_xyz is not None:
            write_bent_pdb(mov_display_pdb, ref_pdb, hkls, AB_xyz,
                            f'{outdir}/bent.pdb')
            geom = check_geometry(f'{outdir}/bent.pdb')
        def _relsymlink(src, dst):
            if os.path.lexists(dst):
                os.unlink(dst)
            os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)
        _relsymlink(os.path.abspath(ref_pdb), f'{outdir}/ref.pdb')
        _relsymlink(os.path.abspath(os.path.join(scan_dir, 'unbent.pdb')),
                    f'{outdir}/unbent.pdb')
        _relsymlink(os.path.abspath(riso_ref_mtz), f'{outdir}/ref.mtz')
        _relsymlink(os.path.abspath(os.path.join(scan_dir, 'unbent.mtz')),
                    f'{outdir}/unbent.mtz')
        # F-space scale (k, B) bent → ref, then compute diff in F-space and
        # inverse-FFT to diff.map.  Absorbs resolution-dependent amplitude
        # mismatch so the diff highlights structural changes, not scaling.
        # Rbent reported here IS post-scaling (riso returned from the scaled fit).
        bent_mtz_path = f'{outdir}/{bent_mtz_name}'
        diff, riso, kF, B = _fspace_scale_and_diff(
            bent_map, ref_h, riso_ref_mtz, riso_ref_col, riso_ref_phi,
            bent_mtz_path, n_cycles=riso_n_cycles, subtract=subtract)
        if diff is None:
            # F-space path failed (too few matched reflections) — fall back to
            # z-scored real-space diff (matches the requested sign convention)
            ref_n     = (ref_d    - ref_d.mean())    / ref_d.std()
            bent_n    = (bent_map - bent_map.mean()) / bent_map.std()
            diff      = (bent_n - ref_n) if subtract == 'ref' else (ref_n - bent_n)
            _write_scan_mtz(bent_map, diff, ref_h, bent_mtz_path)
            riso, kF, B = compute_riso(riso_ref_mtz, riso_ref_col,
                                       bent_mtz_path, 'FDM',
                                       n_cycles=riso_n_cycles,
                                       sigma_cut=riso_sigma_cut)
        diff_norm = diff / diff.std()
        write_ccp4(f'{outdir}/{diffnorm_name}', diff_norm, ref_h)
        # Rbend = R-factor of bent map vs hkl00-resampled map (how much bending
        # changed the resampled mov density).
        if _unbent_bent[0] is None:
            rbend = 0.0  # hkl00 row: bent IS unbent
            _unbent_bent[0] = f'{outdir}/{bent_mtz_name}'
        else:
            rbend, _, _ = compute_riso(_unbent_bent[0], 'FDM',
                                        f'{outdir}/{bent_mtz_name}', 'FDM',
                                        n_cycles=riso_n_cycles,
                                        sigma_cut=riso_sigma_cut)
        peaks = _scan_find_peaks(diff_norm, ref_h, atoms, M)
        # largest peak by absolute σ (sign preserved)
        p_top = peaks['pos'] if abs(peaks['pos'][0]) >= abs(peaks['neg'][0]) else peaks['neg']
        rbent_str  = f'{riso*100:.1f}%'  if riso  is not None else '  N/A'
        rbend_str  = f'{rbend*100:.1f}%' if rbend is not None else '  N/A'
        after_str  = f'{rmsd_after:.3f}'  if rmsd_after  is not None else '  ---'
        atom_str   = _atom_label(p_top[1])
        # Bond/angle RMSZ vs CCP4 monomer-library ideals.  Format:
        # "<bondRMSZ>/<angleRMSZ>".  Refined-quality < 1; loose deposit ~2;
        # geometry breaking down > 5; catastrophic > 100.
        if geom is not None:
            geom_str = f'{geom["bond_rmsz"]:7.2f}/{geom["angle_rmsz"]:5.2f}'
        else:
            geom_str = '        --'

        # ── Dipole content + combined d_opt score ─────────────────────────
        # dipole_score: σ²-weighted mean of per-peak dipole-ness on the
        # top-30 SG-unique peaks in diff_norm (see
        # _dipole_score_from_diff_norm).  Combined score:
        #     S = Rbent + 0.1·RMSD + 0.5·max(0,dip) + 0.05·max(0,bondZ−1)
        # Higher = worse.  argmin(S) across fr-rows is an alternative
        # d_opt pick that combines fit quality, geometry, and
        # misregistration content on a single axis.
        try:
            dipole = _dipole_score_from_diff_norm(
                diff_norm, ref_d,
                gemmi.UnitCell(*ref_h['cell']),
                _sg_name)
        except Exception:
            dipole = None
        S = None
        if riso is not None and rmsd_after is not None:
            _bond_pen = (max(0.0, geom['bond_rmsz'] - 1.0)
                          if geom is not None else 0.0)
            _dip_pen  = max(0.0, dipole) if dipole is not None else 0.0
            S = float(riso + 0.1 * rmsd_after + 0.5 * _dip_pen
                       + 0.05 * _bond_pen)

        row_str = (f"{label:>7}  {after_str:>6}  {n_active:>6d}  "
                   f"{rbent_str:>6}  {rbend_str:>6}  {geom_str:>13}  "
                   f"{p_top[0]:>+7.2f}σ  {atom_str:>20} {p_top[2]:.2f}Å  "
                   f"[{t_elapsed:.0f}s]")
        _log_rows.append(dict(label=label, rmsd=rmsd_after, n_active=n_active,
                              rbent=riso, rbend=rbend, k=kF, B=B,
                              peak_sigma=p_top[0], peak_atom=atom_str,
                              peak_dist=p_top[2], peak_idx=p_top[3],
                              geom=geom,
                              dipole=dipole, S=S,
                              t=t_elapsed, row_str=row_str))
        if verbose:
            print(row_str, flush=True)

    # ── Section 0: pre — raw mov vs ref, NO resolve_altindex, no shift field ─
    # The "old-fashioned" difference map.  hkl00 below is the same idea after
    # resolve_altindex has aligned mov to ref (origin shift / SG op / cell
    # stretch / altindex).  For systems where resolve_altindex took
    # action='none' (lyso, dhfr, myo, raddam — already aligned), pre ≈ hkl00.
    # For systems with action != 'none' (insulin sg_op_origin, porin altindex,
    # lipox altindex+stretch), pre vs hkl00 quantifies the value of the
    # alignment step alone (before any bending).
    outdir_pre = os.path.join(scan_dir, 'pre')
    os.makedirs(outdir_pre, exist_ok=True)
    t_pre = time.time()
    pre_bent_vals = interpolate_map(pre_mov_d, pre_mov_h, ref_pts,
                                     pad_mode=pre_mov_pad_mode)
    pre_bent_map  = pre_bent_vals.reshape(ns2, nr2, nc2)
    pre_bent_map_name = f'bent{col_suffix}.map'
    pre_diff_name     = f'diff_norm{col_suffix}.map'
    pre_bent_mtz_name = f'bent{col_suffix}.mtz'
    write_ccp4(f'{outdir_pre}/{pre_bent_map_name}', pre_bent_map, ref_h)
    # bent.pdb / unbent.pdb / unbent.mtz in pre/ all point to the RAW mov
    # (display copy under the default `preserve_mov_numbering=True`).
    # ref.pdb / ref.mtz still come from the (possibly refmac'd) ref.
    _relsymlink_dir(os.path.abspath(pre_mov_display_pdb),
                    f'{outdir_pre}/bent.pdb')
    _relsymlink_dir(os.path.abspath(ref_pdb),     f'{outdir_pre}/ref.pdb')
    _relsymlink_dir(os.path.abspath(pre_mov_display_pdb),
                    f'{outdir_pre}/unbent.pdb')
    _relsymlink_dir(os.path.abspath(riso_ref_mtz),f'{outdir_pre}/ref.mtz')
    _relsymlink_dir(os.path.abspath(pre_mov_mtz), f'{outdir_pre}/unbent.mtz')
    pre_bent_mtz_path = f'{outdir_pre}/{pre_bent_mtz_name}'
    pre_diff, pre_riso, pre_kF, pre_B = _fspace_scale_and_diff(
        pre_bent_map, ref_h, riso_ref_mtz, riso_ref_col, riso_ref_phi,
        pre_bent_mtz_path, n_cycles=riso_n_cycles, subtract=subtract)
    if pre_diff is None:
        ref_n  = (ref_d        - ref_d.mean())        / ref_d.std()
        bent_n = (pre_bent_map - pre_bent_map.mean()) / pre_bent_map.std()
        pre_diff = (bent_n - ref_n) if subtract == 'ref' else (ref_n - bent_n)
        _write_scan_mtz(pre_bent_map, pre_diff, ref_h, pre_bent_mtz_path)
        pre_riso, pre_kF, pre_B = compute_riso(
            riso_ref_mtz, riso_ref_col, pre_bent_mtz_path, 'FDM',
            n_cycles=riso_n_cycles, sigma_cut=riso_sigma_cut)
    pre_diff_norm = pre_diff / pre_diff.std()
    write_ccp4(f'{outdir_pre}/{pre_diff_name}', pre_diff_norm, ref_h)
    pre_peaks  = _scan_find_peaks(pre_diff_norm, ref_h, pre_atoms, M)
    pre_p_top  = (pre_peaks['pos'] if abs(pre_peaks['pos'][0]) >= abs(pre_peaks['neg'][0])
                  else pre_peaks['neg'])
    pre_atom_s = _atom_label(pre_p_top[1])

    # bondZ/angZ + dipole + S for the pre row — same helpers as _save_point.
    # bondZ from check_geometry on the raw mov PDB (the pre row's model
    # source); dipole from the same _dipole_score_from_diff_norm helper
    # applied to pre_diff_norm + ref_d.
    try:
        pre_geom = check_geometry(pre_mov_pdb)
    except Exception:
        pre_geom = None
    try:
        pre_dipole = _dipole_score_from_diff_norm(
            pre_diff_norm, ref_d, gemmi.UnitCell(*ref_h['cell']), _sg_name)
    except Exception:
        pre_dipole = None
    pre_S = None
    if pre_riso is not None:
        _bond_pen = (max(0.0, pre_geom['bond_rmsz'] - 1.0)
                      if pre_geom is not None else 0.0)
        _dip_pen  = max(0.0, pre_dipole) if pre_dipole is not None else 0.0
        pre_S = float(pre_riso + 0.1 * raw_rmsd + 0.5 * _dip_pen
                       + 0.05 * _bond_pen)

    pre_t_el   = time.time() - t_pre
    pre_rbent_str = f'{pre_riso*100:.1f}%' if pre_riso is not None else '  N/A'
    pre_row_str = (f"{'pre':>7}  {raw_rmsd:>6.3f}  {0:>6d}  "
                   f"{pre_rbent_str:>6}  {'  ---':>6}  {'       --':>13}  "
                   f"{pre_p_top[0]:>+7.2f}σ  {pre_atom_s:>20} {pre_p_top[2]:.2f}Å  "
                   f"[{pre_t_el:.0f}s]")
    _log_rows.append(dict(label='pre', rmsd=raw_rmsd, n_active=0,
                          rbent=pre_riso, rbend=None, k=pre_kF, B=pre_B,
                          peak_sigma=pre_p_top[0], peak_atom=pre_atom_s,
                          peak_dist=pre_p_top[2], peak_idx=pre_p_top[3],
                          geom=pre_geom,
                          dipole=pre_dipole, S=pre_S,
                          t=pre_t_el, row_str=pre_row_str))
    if verbose:
        print(pre_row_str, flush=True)

    # ── Wait for background altindex; post-altindex mov refmac + load ──────
    # Everything below depends on the post-altindex mov (density, atoms,
    # display PDB, hkl00_rmsd).  Section 0 (pre) ran in parallel with
    # altindex — now block and rebind mov_pdb/mov_mtz to whatever
    # resolve_altindex chose.
    if _altindex_fut is not None:
        _res = _altindex_fut.result()
        _altindex_pool.shutdown(wait=True)
        mov_pdb = _res['mov_pdb_out']
        mov_mtz = _res['mov_mtz_out']

    # When altindex returned action='none' (already aligned, no altindex
    # needed), we still need refmac on mov to generate FWT/PHWT if it
    # doesn't already carry them.
    if (run_refinement_flag and mov_mtz.lower().endswith('.mtz')
            and 'FWT' not in {c.label
                              for c in gemmi.read_mtz_file(mov_mtz).columns}):
        mov_pdb, mov_mtz = run_refinement(
            mov_pdb, mov_mtz,
            outdir=os.path.join(scan_dir, 'refine_mov'),
            n_cycles=refine_cycles, fill_asu=fill_asu,
            fill_method=fill_method)

    # Load post-altindex mov density.  Default mov_fullcell to True for
    # CCP4 input (avoids ASU-boundary noise in diff_norm); MTZ input
    # ignores the flag (already full-cell).
    if mov_fullcell is None:
        mov_fullcell = mov_mtz.lower().endswith(('.map', '.ccp4', '.mrc'))
    mov_d, mov_h, mov_mtz_resolved, mov_f_col = _load(mov_mtz, 'mov',
                                                        fullcell=mov_fullcell)
    col_suffix = ('' if (mov_f_col is None or mov_f_col in _CANONICAL_F)
                  else f'_{mov_f_col}')
    pad_mode   = ('wrap' if (mov_mtz_resolved is not None or mov_fullcell)
                  else 'reflect')

    # mov_display_pdb + `atoms` (post-altindex; used by _save_point and by
    # write_bent_pdb via closure).
    if preserve_mov_numbering and norm_plan is not None:
        mov_display_pdb = _reverse_normalize_pdb(
            mov_pdb, norm_plan,
            os.path.join(scan_dir, '_mov_display.pdb'))
    else:
        mov_display_pdb = mov_pdb
    atoms = (_scan_read_pdb_atoms_p1(ref_pdb, src='ref',
                                      cell_override=_ref_cell)
             + _scan_read_pdb_atoms_p1(mov_display_pdb, src='mov',
                                        cell_override=_ref_cell))

    # hkl00 RMSD — uses post-altindex mov_pdb.
    precomputed_origin = None
    with open(os.devnull, 'w') as _dev, _contextlib.redirect_stdout(_dev):
        _a1o, _cell1o, _sg1o = expand_to_p1(mov_pdb)
        _a2o, _cell2o, _sg2o = expand_to_p1(ref_pdb)
        _fmo, _cao, _uido, _bfo = match_atoms(_a1o, _a2o)
        _fmo, _cao, _, _, _ = reject_outliers(_fmo, _cao, _uido, _cell1o,
                                               mad_sigma=outlier_sigma,
                                               b_sigma=b_sigma, bfacs=_bfo)
    hkl00_rmsd = rmsd_ca(_fmo, _cao, _h0, _AB0, _cell2o)

    # ── Section 1: hkl00 — zero shift field (post-resolve baseline) ──────────
    outdir0 = os.path.join(scan_dir, 'hkl00')
    os.makedirs(outdir0, exist_ok=True)
    t0 = time.time()
    bent_vals = interpolate_map(mov_d, mov_h, _map_query(ref_pts),
                                 pad_mode=pad_mode)
    bent_map0 = bent_vals.reshape(ns2, nr2, nc2)
    _save_point('hkl00', outdir0, bent_map0, hkl00_rmsd, hkl00_rmsd, 0,
                time.time() - t0, hkls=_h0, AB_xyz=_AB0)

    # unbent.{pdb,mtz}: symlink to the post-resolve mov PDB and MTZ
    # (whatever resolve_altindex chose — refine_mov/<stem>_refined.* for
    # the simple case, or altindex_resolve/<stem>_sgop.* / *_alt_refined.*
    # for the SG-op / true-altindex cases).  Both files are in mov's cell
    # and are bit-identical to the refmac output (no gemmi re-write).  For
    # mov cells that differ from ref's (e.g. insulin H 3) the unbent pair
    # won't overlay with ref structures in Coot — use the per-subdir
    # bent.{pdb,mtz} files for that overlay (they're written on ref's grid
    # by definition).
    _relsymlink_dir(os.path.abspath(mov_display_pdb),
                     os.path.join(scan_dir, 'unbent.pdb'))
    if mov_mtz_resolved is not None:
        _relsymlink_dir(os.path.abspath(mov_mtz_resolved),
                         os.path.join(scan_dir, 'unbent.mtz'))
    else:
        unbent_mtz_path = os.path.join(scan_dir, 'unbent.mtz')
        if not os.path.exists(unbent_mtz_path):
            if verbose:
                print(f'  Converting mov map → MTZ for unbent.mtz...',
                      flush=True)
            _map2mtz(mov_mtz, unbent_mtz_path)

    # ── Section 2: hkl01..hklNN — progressive, one canon HKL at a time ───────
    _last_hkl_n = [0]   # track last n_non_dc saved; list so closure can mutate

    def _hkl_callback(iter_i, nhkls, n_canon, rmsd, hkls_now, AB_xyz, active, snr):
        n_non_dc = n_canon - 1
        if n_non_dc < 1 or n_non_dc > max_hkl_scan:
            return
        if n_non_dc == _last_hkl_n[0]:
            return   # same canonical count as last save — skip duplicate
        t0_cb = time.time()
        label  = f'hkl{n_non_dc:02d}'
        outdir = os.path.join(scan_dir, label)
        os.makedirs(outdir, exist_ok=True)
        a1pts     = _atoms1_pts(ref_pts)
        delta     = _eval_chunked(a1pts, hkls_now, AB_xyz, chunk_size)
        bent_vals = interpolate_map(mov_d, mov_h,
                                     _map_query(ref_pts, delta),
                                     pad_mode=pad_mode)
        bent_map  = bent_vals.reshape(ns2, nr2, nc2)
        _save_point(label, outdir, bent_map, hkl00_rmsd, rmsd,
                    int(active.sum()), time.time() - t0_cb,
                    hkls=hkls_now, AB_xyz=AB_xyz)
        save_fitparams(os.path.join(outdir, 'PSDVF.mtz'),
                       hkls_now, AB_xyz, active, snr,
                       ref_h['cell'], ref_h['cell'], 'xyz', rmsd)
        _last_hkl_n[0] = n_non_dc

    bend_fit_progressive(
        mov_pdb, ref_pdb,
        fitreso_start=100.0, fitreso_end=2.0,
        max_canon=max_hkl_scan + 1, batch_hkls=1, od_margin=od_margin,
        drop_snr=drop_snr, outlier_sigma=outlier_sigma, b_sigma=b_sigma, atom_sel=atom_sel, bound_by_obs=bound_by_obs, pnn_mode=pnn_mode,
        verbose=False,
        iter_callback=_hkl_callback,
        precomputed_origin=precomputed_origin,
    )

    # ── Section 3: fr<N> — resolution scans (single progressive call) ─────────
    # All fr-rows share ONE bend_fit_progressive call at fitreso_end =
    # min(fitreso_list).  We pre-compute the exact n_used (= 1 DC + count of
    # non-DC HKLs with d ≥ thr) for each requested fr-row from the actual HKL
    # pool, pass it as iteration_schedule, and snapshot in iter_callback —
    # each iteration lands exactly on one fr-row.  Result is bit-equivalent
    # to N independent bend_fit_progressive calls (same HKL set per fr-row,
    # same SVD per iter), but eliminates the redundant coarse-iter work each
    # separate call repeats.  For lipox (~4× cell expansion, ~6700 P1 atoms)
    # this cuts fr5 wall time from days to hours.
    _early_stop_reason = [None]     # written by _fr_callback, read by log footer
    if verbose:
        print('-' * 122, flush=True)
    if fitreso_list:
        fr_thresholds_desc = sorted([float(x) for x in fitreso_list],
                                     reverse=True)   # coarsest first

        # Build the actual HKL pool to fr_min and count HKLs at d ≥ thr for
        # each requested fr-row.  Must match generate_hkls inside
        # bend_fit_progressive: same cell (ref/mov are the same protein +
        # post-resolve so cell1 from mov_pdb is the relevant one), same
        # pool_reso = fr_min.  Computed once, mapped to n_used = count + 1
        # (+1 for the DC term at all_hkls[0]).
        fr_min = fr_thresholds_desc[-1]
        from gemmi import read_structure as _read_structure
        st_mov = _read_structure(mov_pdb)
        st_mov.setup_entities()
        cell_for_pool = (st_mov.cell.a, st_mov.cell.b, st_mov.cell.c,
                         st_mov.cell.alpha, st_mov.cell.beta, st_mov.cell.gamma)
        pool_hkls = generate_hkls(cell_for_pool, fr_min)
        Gs_pool   = _reciprocal_metric(cell_for_pool)
        hv_pool   = pool_hkls[1:].astype(float)
        inv_d2    = np.sum((hv_pool @ Gs_pool) * hv_pool, axis=1)
        d_pool    = np.where(inv_d2 > 0, 1.0 / np.sqrt(inv_d2), np.inf)

        # thr → n_used (sorted ascending).  Multiple fr-rows may collapse to
        # the same n_used (e.g. fr6 and fr5 both exhaust the pool) — fr_thrs
        # maps each n_used to ALL thresholds that resolve to it.
        n_to_thrs = {}
        for thr in fr_thresholds_desc:
            n_used_thr = int(np.sum(d_pool >= thr - 1e-9)) + 1
            n_used_thr = max(n_used_thr, 2)              # at least DC + 1
            n_used_thr = min(n_used_thr, len(pool_hkls))
            n_to_thrs.setdefault(n_used_thr, []).append(thr)
        schedule = sorted(n_to_thrs.keys())

        t_fit0 = time.time()
        fr_done = set()
        final_state = [None]

        def _fr_label(thr):
            return f'fr{int(thr) if thr == int(thr) else thr}'

        def _save_fr(thr, hkls, AB_xyz, active, snr, rmsd, t_cum):
            label  = _fr_label(thr)
            outdir = os.path.join(scan_dir, label)
            os.makedirs(outdir, exist_ok=True)
            a1pts     = _atoms1_pts(ref_pts)
            delta     = _eval_chunked(a1pts, hkls, AB_xyz, chunk_size)
            bent_vals = interpolate_map(mov_d, mov_h,
                                         _map_query(ref_pts, delta),
                                         pad_mode=pad_mode)
            bent_map  = bent_vals.reshape(ns2, nr2, nc2)
            _save_point(label, outdir, bent_map, hkl00_rmsd, rmsd,
                        int(active.sum()), t_cum,
                        hkls=hkls, AB_xyz=AB_xyz)
            save_fitparams(os.path.join(outdir, 'PSDVF.mtz'),
                           hkls, AB_xyz, active, snr,
                           ref_h['cell'], ref_h['cell'], 'xyz', rmsd)

        # Early-stop machinery (default: scan_all_fr=False).  When the
        # combined `score` column has risen for 2 consecutive fr-rows
        # relative to the running argmin, stop the progressive fit — we've
        # already passed the score minimum and further fine-fitreso rows
        # take much longer while adding no value.  `scan_all_fr=True`
        # disables this and runs the full 8-row fr schedule (used by the
        # test gamut for reference).
        class _EarlyStop(Exception):
            pass
        _stop_reason = _early_stop_reason

        def _fr_callback(iter_i, nhkls, n_canon, rmsd,
                         hkls_now, AB_xyz, active, snr):
            t_cum = time.time() - t_fit0
            final_state[0] = (hkls_now, AB_xyz, active, snr, rmsd, t_cum)
            # This iter's n_used = nhkls + 1 (nhkls is the non-DC count); save
            # every fr-row whose schedule entry matches.  od_margin may stop
            # the loop early — remaining fr-rows are handled by the final-
            # state fallback below.
            n_used = nhkls + 1
            for thr in n_to_thrs.get(n_used, []):
                if thr in fr_done:
                    continue
                _save_fr(thr, hkls_now, AB_xyz, active, snr, rmsd, t_cum)
                fr_done.add(thr)

            # Early-stop check.  Runs after each fr-row is saved; needs at
            # least 3 fr-rows before it can trigger (need a baseline + 2
            # worse-than-baseline in a row).
            if scan_all_fr:
                return
            fr_scored = [r for r in _log_rows
                         if r['label'].startswith('fr')
                         and r.get('S') is not None]
            if len(fr_scored) < 3:
                return
            best_S = min(r['S'] for r in fr_scored)
            # Relative tolerance around the running argmin; rows within
            # this fraction of the min don't count as "worse".  Default
            # 1% (early_stop_tol=0.01) — tune coarser for looser
            # stopping (0.02 = 2%), tighter for aggressive stopping.
            tol    = max(1e-4, early_stop_tol * abs(best_S))
            # Count consecutive tail rows worse than best (beyond tol).
            n_worse = 0
            for r in reversed(fr_scored):
                if r['S'] > best_S + tol:
                    n_worse += 1
                else:
                    break
            if n_worse >= early_stop_n:
                _stop_reason[0] = (
                    f"score rose above argmin {best_S:.3f} for {n_worse} "
                    f"consecutive rows — skipping fr-rows finer than "
                    f"{fr_scored[-1]['label']}")
                raise _EarlyStop()

        try:
            bend_fit_progressive(
                mov_pdb, ref_pdb,
                fitreso_start=20.0, fitreso_end=fr_min,
                drop_snr=drop_snr, batch_hkls=batch_hkls, od_margin=od_margin,
                outlier_sigma=outlier_sigma, b_sigma=b_sigma, atom_sel=atom_sel, bound_by_obs=bound_by_obs, pnn_mode=pnn_mode,
                verbose=False,
                iter_callback=_fr_callback,
                precomputed_origin=precomputed_origin,
                iteration_schedule=schedule,
            )
        except _EarlyStop:
            if verbose:
                print(f'  early-stop: {_stop_reason[0]}', flush=True)

        # od_margin (or early-stop) stopped us before reaching some
        # thresholds → for the OD-stop case, snapshot those with the final
        # iter's state (the deepest fit we got, same as what a separate
        # fitreso_end=thr call would produce under the same OD stop).
        # For early-stop, DO NOT snapshot — the whole point is those rows
        # aren't worth computing, and filling them with the deepest state
        # would just clone the argmin row's fit under a lie of a label.
        if _stop_reason[0] is None:
            for thr in fr_thresholds_desc:
                if thr in fr_done:
                    continue
                if final_state[0] is not None:
                    _save_fr(thr, *final_state[0])
                    fr_done.add(thr)

    # ── Section 4: best — parabola fit of Rbent vs 1/d², re-run at vertex ────
    # Across all systems, Rbent has a clear U-shape vs fitreso: dropping with
    # finer resolution as the shift field absorbs more structural detail,
    # then climbing again as high-frequency HKLs add noise outside the smooth
    # PSDVF's natural bandwidth.  Fit a parabola in x=1/d² to the 3-5 fr-rows
    # nearest the empirical minimum, find the vertex, and re-fit at d_opt.
    fr_rows_all = [r for r in _log_rows
                    if r['label'].startswith('fr') and r['rbent'] is not None]
    fr_resos_all = []
    for r in fr_rows_all:
        try:
            fr_resos_all.append(float(r['label'][2:]))
        except ValueError:
            fr_resos_all.append(None)

    # ── RMSD-baseline filter ─────────────────────────────────────────────
    # The fr-row list is in admission order (coarse → fine).  Walking
    # forward, accept rows as long as their CA RMSD remains ≤ the hkl00
    # (unbent) baseline.  A row with RMSD > baseline means the field is
    # actively pulling CAs further from their targets than no-bending would
    # — clearly the wrong direction; stop the parabola fit there.  We do
    # NOT require strict monotone decrease, because a parabola minimum
    # legitimately sits between two rows where one is slightly higher than
    # the other; both can still be useful.
    hkl00_rmsd = next((r['rmsd'] for r in _log_rows
                       if r['label'] == 'hkl00' and r['rmsd'] is not None), None)
    fr_rows  = []
    fr_resos = []
    cliff_reason = None
    for r, d in zip(fr_rows_all, fr_resos_all):
        if d is None or r['rmsd'] is None:
            continue
        if hkl00_rmsd is not None and r['rmsd'] > hkl00_rmsd:
            cliff_reason = (f'{r["label"]} RMSD {r["rmsd"]:.3f} > baseline '
                            f'{hkl00_rmsd:.3f}')
            break
        fr_rows.append(r)
        fr_resos.append(d)
    if verbose and cliff_reason:
        print(f'\nbest filter: stopped at RMSD cliff ({cliff_reason}); '
              f'using {len(fr_rows)} of {len(fr_rows_all)} fr-rows',
              flush=True)

    fr_pts = [(d, r['rbent']) for d, r in zip(fr_resos, fr_rows) if d is not None]
    d_opt   = None
    r_pred  = None
    clamped = ''
    ds_used, lo, hi = None, 0, 0

    if len(fr_pts) >= 3:
        ds  = np.array([p[0] for p in fr_pts], dtype=float)
        rbs = np.array([p[1] for p in fr_pts], dtype=float)
        i_min = int(np.argmin(rbs))
        # 3-5 points centered on argmin, clipped to bracket
        lo = max(0, i_min - 2)
        hi = min(len(ds), i_min + 3)
        if hi - lo < 3:                      # extend if argmin is at an end
            if lo == 0:
                hi = min(len(ds), lo + 3)
            else:
                lo = max(0, hi - 3)
        x = 1.0 / (ds[lo:hi] ** 2)
        y = rbs[lo:hi]
        a, b, c = np.polyfit(x, y, 2)        # y = a x² + b x + c (np convention)
        if a > 0:
            x_vert_raw = -b / (2.0 * a)
            x_vert     = float(np.clip(x_vert_raw, x.min(), x.max()))
            d_opt      = float(np.sqrt(1.0 / x_vert))
            r_pred     = a * x_vert**2 + b * x_vert + c
            clamped    = ' (clamped to bracket — true vertex outside scan)' \
                         if x_vert != x_vert_raw else ''
        else:
            # Concave-down — no interior min; fall back to argmin scan point
            d_opt  = float(ds[i_min])
            r_pred = float(rbs[i_min])
            clamped = ' (parabola concave-down; using argmin)'
        ds_used = ds[lo:hi]
        if verbose:
            print(f'\nbest-fit parabola over {hi - lo} fr-rows '
                  f'({", ".join(f"fr{int(d) if d == int(d) else d}" for d in ds[lo:hi])}): '
                  f'd_opt = {d_opt:.2f} Å, predicted Rbent = {r_pred*100:.1f}%{clamped}',
                  flush=True)
    elif len(fr_rows_all) > 0:
        # Not enough RMSD-monotone rows for a parabola — use the row with
        # the smallest RMSD across ALL fr-rows as the best fitreso end.
        rmsd_arr = [r['rmsd'] for r in fr_rows_all if r['rmsd'] is not None]
        i_min_rmsd = int(np.argmin(rmsd_arr))
        d_opt  = float(fr_resos_all[i_min_rmsd])
        r_pred = float(fr_rows_all[i_min_rmsd]['rbent'])
        clamped = (f' (parabola skipped — only {len(fr_rows)} RMSD-monotone '
                    f'rows; using argmin-RMSD fr-row)')
        ds_used = np.array([d_opt])
        lo, hi = 0, 1
        if verbose:
            print(f'\nbest: argmin(RMSD) fallback at d={d_opt:.2f} Å, '
                  f'Rbent={r_pred*100:.1f}%{clamped}', flush=True)

    # If d_opt was determined by either branch, re-fit at it and save the
    # 'best' point with full per-point outputs.  Otherwise (no fr-rows at
    # all — shouldn't happen) skip the best section.
    if d_opt is not None:
        outdir = os.path.join(scan_dir, 'best')
        os.makedirs(outdir, exist_ok=True)
        t0 = time.time()
        result = bend_fit_progressive(
            mov_pdb, ref_pdb,
            fitreso_start=20.0, fitreso_end=d_opt,
            drop_snr=drop_snr, batch_hkls=batch_hkls, od_margin=od_margin,
            outlier_sigma=outlier_sigma, b_sigma=b_sigma, atom_sel=atom_sel, bound_by_obs=bound_by_obs, pnn_mode=pnn_mode,
            verbose=False,
            precomputed_origin=precomputed_origin,
        )
        t_fit = time.time() - t0
        dim_list = list(result.dimensions)
        xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
        AB_xyz   = np.zeros((3, len(result.hkls), 2))
        for new_i, old_i in enumerate(xyz_dims):
            AB_xyz[new_i] = result.AB[old_i]
        a1pts     = _atoms1_pts(ref_pts)
        delta     = _eval_chunked(a1pts, result.hkls, AB_xyz, chunk_size)
        bent_vals = interpolate_map(mov_d, mov_h,
                                     _map_query(ref_pts, delta),
                                     pad_mode=pad_mode)
        bent_map  = bent_vals.reshape(ns2, nr2, nc2)
        _save_point('best', outdir, bent_map, hkl00_rmsd, result.rmsd,
                    int(result.active.sum()), t_fit,
                    hkls=result.hkls, AB_xyz=AB_xyz)
        save_fitparams(f'{outdir}/PSDVF.mtz',
                       result.hkls, result.AB, result.active, result.snr,
                       result.cell1, result.cell2, result.dimensions, result.rmsd)
        _best_meta = dict(d_opt=d_opt, r_pred=r_pred, n_fit=int(hi - lo),
                          fr_used=[float(d) for d in (ds_used if ds_used is not None
                                                       else [])],
                          clamped=clamped)
        # (Best-row pre-σ is now filled in the per-row loop below.)
    else:
        _best_meta = None

    # ── Per-row pre@peak: σ in pre/diff_norm.map at each row's own top-peak
    # voxel.  Small |pre@peak| relative to |peak_sigma| ⇒ the row's peak
    # was hidden pre-bending and revealed by alignment+bending.  Similar
    # magnitude ⇒ feature was already visible pre-bending and neither
    # bending revealed nor damped it (chemistry or persistent residual).
    # Trivially equals peak_sigma for the pre row and ~equal for hkl00.
    for _row in _log_rows:
        pidx = _row.get('peak_idx')
        if pidx is None:
            continue
        i, j, k = pidx
        _row['pre_at_peak_sigma'] = float(pre_diff_norm[i, j, k])
    if verbose and _best_meta is not None:
        _br = _log_rows[-1]
        if _br.get('pre_at_peak_sigma') is not None:
            print(f"  pre diff-map σ at best-peak voxel: "
                  f"{_br['pre_at_peak_sigma']:+.2f}σ "
                  f"(best peak: {_br['peak_sigma']:+.2f}σ at "
                  f"{_br['peak_atom']})", flush=True)

    # ── Final summary log ───────────────────────────────────────────────────
    _sign_note = ('diff = bent − ref  (+peak ⇒ density in bent absent from ref)'
                  if subtract == 'ref' else
                  'diff = ref − bent  (+peak ⇒ density in ref absent from bent)')
    log_path = os.path.join(scan_dir, 'scan_fitreso.log')
    with open(log_path, 'w') as fh:
        fh.write(f"# fitreso_scan summary\n")
        fh.write(f"# mov_pdb : {mov_pdb}\n")
        fh.write(f"# ref_pdb : {ref_pdb}\n")
        fh.write(f"# mov_mtz : {mov_mtz}\n")
        fh.write(f"# ref_mtz : {ref_mtz}\n")
        fh.write(f"# scan_dir: {scan_dir}\n")
        fh.write(f"# sign    : {_sign_note}\n")
        fh.write(f"# Rbent   : Rfac after F-space (k+B) scaling of bent → ref\n")
        fh.write(f"# Rbend   : Rfac of bent vs hkl00 bent (how much the PSDVF moved density)\n")
        fh.write(f"# k, B    : F-space scale + isotropic B (Å²) used for Rbent\n")
        fh.write(f"# pre@peak: σ in pre/diff_norm.map at THIS row's top-peak voxel\n")
        fh.write(f"#           (|pre@peak| ≪ |peak| ⇒ feature hidden pre-bending, revealed by fit;\n")
        fh.write(f"#            |pre@peak| ≈ |peak|  ⇒ feature was already visible pre-bending)\n\n")
        fh.write(f"# bondZ/  : RMSZ of bond lengths / bond angles vs CCP4 monomer-lib\n")
        fh.write(f"# angZ      ideals (gemmi.prepare_topology).  Refined ~<1; loose ~2;\n")
        fh.write(f"#           breaking >5; catastrophic >100 — bent.pdb is unphysical.\n")
        fh.write(f"# dipole  : σ²-weighted mean of per-peak dipole-ness on the top-30 diff peaks.\n")
        fh.write(f"#           ~0 = no misregistration/embossing content; +1 = perfect dipoles.\n")
        fh.write(f"# score   : combined d_opt score, score = Rbent + 0.1·RMSD + 0.5·max(0,dipole)\n")
        fh.write(f"#           + 0.05·max(0,bondZ−1).  Lower = better; argmin(score) across fr-rows\n")
        fh.write(f"#           is an alternative d_opt pick that combines fit / geometry / align.\n\n")
        hdr = (f"{'label':>7}  {'RMSD':>6}  {'active':>6}  "
               f"{'Rbent':>6}  {'Rbend':>6}  "
               f"{'k':>6}  {'B(Å²)':>7}  "
               f"{'bondZ':>7}/{'angZ':>5}  "
               f"{'dipole':>7}  {'score':>6}  "
               f"{'peak':>7}  {'atom':>20}  {'dist':>6}  "
               f"{'pre@peak':>9}  {'t(s)':>5}\n")
        fh.write(hdr)
        fh.write('-' * (len(hdr) - 1) + '\n')
        for r in _log_rows:
            rmsd_s = f"{r['rmsd']:.3f}" if r['rmsd']  is not None else '  ---'
            rbe_s  = f"{r['rbent']*100:.1f}%" if r['rbent'] is not None else '  N/A'
            rbd_s  = f"{r['rbend']*100:.1f}%" if r['rbend'] is not None else '  N/A'
            k_s    = f"{r['k']:.3f}"   if r['k'] is not None else '  N/A'
            B_s    = f"{r['B']:+6.2f}" if r['B'] is not None else '   N/A'
            g      = r.get('geom')
            if g is not None:
                g_s = f"{g['bond_rmsz']:7.2f}/{g['angle_rmsz']:5.2f}"
            else:
                g_s = f"{'--':>7}/{'--':>5}"
            dip    = r.get('dipole')
            dip_s  = f"{dip:>+7.3f}" if dip is not None else f"{'—':>7}"
            S_val  = r.get('S')
            S_s    = f"{S_val:>6.3f}" if S_val is not None else f"{'—':>6}"
            pab    = r.get('pre_at_peak_sigma')
            pab_s  = f"{pab:>+7.2f}σ" if pab is not None else f"{'—':>9}"
            fh.write(f"{r['label']:>7}  {rmsd_s:>6}  {r['n_active']:>6d}  "
                     f"{rbe_s:>6}  {rbd_s:>6}  "
                     f"{k_s:>6}  {B_s:>7}  "
                     f"{g_s}  "
                     f"{dip_s}  {S_s}  "
                     f"{r['peak_sigma']:>+7.2f}σ  {r['peak_atom']:>20}  "
                     f"{r['peak_dist']:>5.2f}Å  "
                     f"{pab_s:>9}  {r['t']:>4.0f}s\n")

        if _early_stop_reason[0] is not None:
            fh.write(f"\n# early-stop: {_early_stop_reason[0]}\n")

        # ── Alternative d_opt pick from argmin(score) across fr-rows ──────
        fr_score_rows = [r for r in _log_rows
                         if r['label'].startswith('fr') and r.get('S') is not None]
        if fr_score_rows:
            best_score = min(fr_score_rows, key=lambda r: r['S'])
            fh.write(f"\n# combined-score d_opt: {best_score['label']}  "
                     f"score={best_score['S']:.3f}  "
                     f"(Rbent={best_score['rbent']*100:.1f}%, "
                     f"RMSD={best_score['rmsd']:.3f}Å, "
                     f"dipole={best_score.get('dipole') or 0.0:+.3f}, "
                     f"bondZ={(best_score['geom']['bond_rmsz'] if best_score.get('geom') else 0.0):.2f})\n")

            # Replay early-stop rule over the fr rows so debug/full runs
            # still show where the (default) early-stop would have fired
            # and which row it would have picked as d_opt.  Same rule as
            # the live callback: after each row, if the *last 2* rows are
            # worse than the current running min by >1% relative, stop.
            _run_best_S   = float('inf')
            _run_best_row = None
            for _i, _r in enumerate(fr_score_rows):
                _tol = max(1e-4, early_stop_tol * abs(_run_best_S)
                                if _run_best_S != float('inf') else 1e-4)
                if _r['S'] < _run_best_S - _tol:
                    _run_best_S = _r['S']
                    _run_best_row = _r
                _seen = fr_score_rows[:_i+1]
                if len(_seen) >= early_stop_n + 1:
                    _n_worse = 0
                    for _rt in reversed(_seen):
                        if _rt['S'] > _run_best_S + _tol:
                            _n_worse += 1
                        else:
                            break
                    if _n_worse >= early_stop_n:
                        _same = (_run_best_row['label'] == best_score['label'])
                        _agree = ' (same as full-scan argmin)' if _same else \
                                 f' (full-scan argmin was {best_score["label"]})'
                        fh.write(f"# early-stop would have fired at "
                                 f"{_r['label']} — picked "
                                 f"{_run_best_row['label']} score="
                                 f"{_run_best_S:.3f}{_agree} "
                                 f"[tol={early_stop_tol}, n={early_stop_n}]\n")
                        break
            else:
                fh.write(f"# early-stop would NOT have fired — every row "
                         f"stayed within {100*early_stop_tol:g}% of the "
                         f"running argmin ({best_score['label']} score="
                         f"{best_score['S']:.3f}) "
                         f"[tol={early_stop_tol}, n={early_stop_n}]\n")

        if _best_meta is not None:
            fr_used_str = ', '.join(f'fr{int(d) if d == int(d) else d}'
                                     for d in _best_meta['fr_used'])
            fh.write(f"\n# best row: parabola in 1/d² over {_best_meta['n_fit']} "
                     f"rows ({fr_used_str}); vertex at d_opt = "
                     f"{_best_meta['d_opt']:.2f} Å"
                     f"{_best_meta['clamped']}"
                     f", predicted Rbent = "
                     f"{_best_meta['r_pred']*100:.1f}%\n")
    if verbose:
        print(f'\nDone.  Summary table: {log_path}', flush=True)


def compute_riso(ref_mtz_path, ref_col, test_mtz_path, test_col,
                 n_cycles=4, sigma_cut=float('inf'), anisotropic=False):
    """Scale F_test to F_ref with k + isotropic (default) or anisotropic B;
    return (riso, k, B_iso_eq) or (None, None, None).

    Model (same as CCP4 scaleit):

        F_test_scaled = k · exp(-Q(h)) · F_test
        Q(h) = (B/4) · s²                       (isotropic, s² = 1/d²)
        Q(h) = h_orth^T · V · h_orth             (anisotropic; V symmetric 3×3,
                                                  B_iso_eq = (4/3)·trace(V))

    NB: the fit objective here is **not** the same as scaleit's.  Scaleit does
    Fox-Holmes intensity-space LS while refining both SC(1) and SC(2) — and
    the (SC(1), SC(2)) → (α·SC(1), α·SC(2)) degeneracy is broken only by a
    rank-revealing solver (EIGSOL with ELIM=3.5e-4), so its reported (k, B)
    sit somewhere arbitrary along the degenerate ray.

    Here we instead fix SC(1)=1 and refine k (=SC(2)) and B (or V) by F-space
    nonlinear least squares (scipy least_squares, LM) — well-posed and
    directly minimizes the F-space (= real-space L2, by Parseval) residual.
    Iterated n_cycles times with sigma_cut · rms-residual outlier rejection
    between cycles.  Default sigma_cut=inf keeps all reflections (no reliable
    σ for map-derived F — diff.com's `EXCLUDE SIG 3` is intentionally skipped).

    The reported B will not match scaleit's exactly; expect F-space-LS to give
    a *lower* F-Rfac for the same model.

    Riso = Σ|F_ref − F_test_scaled| / Σ|F_ref| over remaining reflections.
    """
    try:
        from scipy.optimize import least_squares

        ref_mtz = gemmi.read_mtz_file(ref_mtz_path)
        tst_mtz = gemmi.read_mtz_file(test_mtz_path)
        cell    = ref_mtz.cell

        ref_data = np.array(ref_mtz.array)
        tst_data = np.array(tst_mtz.array)
        ref_lbls = ref_mtz.column_labels()
        tst_lbls = tst_mtz.column_labels()
        ref_hkl  = ref_data[:, :3].astype(int)
        tst_hkl  = tst_data[:, :3].astype(int)
        ref_F    = ref_data[:, ref_lbls.index(ref_col)]
        tst_F    = tst_data[:, tst_lbls.index(test_col)]

        # Match by HKL
        tdict = {(int(h[0]), int(h[1]), int(h[2])): i
                 for i, h in enumerate(tst_hkl)}
        keep_idx = [tdict.get((int(h[0]), int(h[1]), int(h[2])))
                    for h in ref_hkl]
        mask = np.array([k is not None for k in keep_idx])
        hkl  = ref_hkl[mask].astype(float)
        r    = ref_F[mask]
        t    = np.array([tst_F[i] for i in keep_idx if i is not None])

        with np.errstate(invalid='ignore'):
            ok = (r > 0) & (t > 0) & np.isfinite(r) & np.isfinite(t)
        hkl, r, t = hkl[ok], r[ok], t[ok]
        if len(r) < 10:
            return None, None, None

        # Reciprocal-cartesian h vectors: h_orth_row = hkl_row · M_frac
        # (cell.fractionalization_matrix takes direct cart→frac; its rows are
        # the reciprocal cartesian basis vectors a*, b*, c*.)
        M_frac = np.array(cell.frac.mat.tolist())
        h_orth = hkl @ M_frac           # (N, 3)
        s_sq   = np.sum(h_orth ** 2, axis=1)   # 1/d²

        use = np.ones(len(r), dtype=bool)

        # initial k from Karle-Hauptman; B (or V) = 0
        k0 = float(np.dot(r, t) / np.dot(t, t))

        if not anisotropic:
            params = np.array([np.log(max(k0, 1e-10)), 0.0])

            def _scale(p, mask_sub=None):
                lk, B = p
                ss = s_sq if mask_sub is None else s_sq[mask_sub]
                return np.exp(lk) * np.exp(-B / 4.0 * ss)
        else:
            params = np.array([np.log(max(k0, 1e-10)),
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

            def _Vmat(p):
                V11, V22, V33, V12, V13, V23 = p[1:]
                return np.array([[V11, V12, V13],
                                 [V12, V22, V23],
                                 [V13, V23, V33]])

            def _scale(p, mask_sub=None):
                lk = p[0]
                V  = _Vmat(p)
                hh = h_orth if mask_sub is None else h_orth[mask_sub]
                Q  = np.einsum('ij,ni,nj->n', V, hh, hh)
                return np.exp(lk) * np.exp(-Q)

        for _ in range(n_cycles):
            def _resid(p):
                return r[use] - _scale(p, mask_sub=use) * t[use]
            try:
                res = least_squares(_resid, params, method='lm', max_nfev=200)
                params = res.x
            except Exception:
                break
            if np.isfinite(sigma_cut):
                resid_all = np.abs(r - _scale(params) * t)
                rms = float(np.std(resid_all[use])) if use.sum() > 1 else np.inf
                use = resid_all < sigma_cut * rms
                if use.sum() < 10:
                    break

        scaled = _scale(params) * t
        riso   = float(np.sum(np.abs(r - scaled)) / np.sum(np.abs(r)))
        k_fit  = float(np.exp(params[0]))
        if anisotropic:
            B_iso_eq = float(4.0 * np.trace(_Vmat(params)) / 3.0)
        else:
            B_iso_eq = float(params[1])
        return riso, k_fit, B_iso_eq
    except Exception:
        return None, None, None


# ══════════════════════════════════════════════════════════════════════════════
# MTZ helpers — column detection, inverse FFT, and optional refinement
# ══════════════════════════════════════════════════════════════════════════════

_CANONICAL_2FOFC_PAIRS = [('FWT', 'PHWT'), ('2FOFCWT', 'PH2FOFCWT')]


def detect_2fofc_cols(mtz_path, f_col=None, phi_col=None):
    """Return (f_col, phi_col) suitable for a likelihood-weighted 2Fo-Fc map.

    If both f_col and phi_col are supplied, validate they exist and return them.
    Otherwise auto-detect from canonical pairs (FWT/PHWT, then 2FOFCWT/PH2FOFCWT).
    Raises ValueError with a helpful message if no acceptable columns are found.
    """
    mtz    = gemmi.read_mtz_file(mtz_path)
    labels = {c.label for c in mtz.columns}

    if f_col is not None and phi_col is not None:
        for c in (f_col, phi_col):
            if c not in labels:
                raise ValueError(
                    f"Column '{c}' not found in {mtz_path}; "
                    f"available: {sorted(labels)}")
        return f_col, phi_col

    for f, ph in _CANONICAL_2FOFC_PAIRS:
        if f in labels and ph in labels:
            return f, ph

    raise ValueError(
        f"No likelihood-weighted 2Fo-Fc columns in {mtz_path}.\n"
        f"Tried: {_CANONICAL_2FOFC_PAIRS}\n"
        f"Available: {sorted(labels)}\n"
        f"Use labels=F,PHI (CLI) to specify columns, "
        f"or use run_refinement to generate FWT/PHWT first.")


def mtz_to_map_data(mtz_path, f_col, phi_col, sample_rate=0.0):
    """Inverse-FFT an MTZ F/PHI column pair; return (ndarray, header_dict).

    The returned header matches the read_ccp4 format: full unit cell
    (starts=0, mapc=1/mapr=2/maps=3).  No mov_fullcell workaround needed.
    """
    mtz  = gemmi.read_mtz_file(mtz_path)
    grid = mtz.transform_f_phi_to_map(f_col, phi_col, sample_rate=sample_rate)
    # gemmi gives (nu, nv, nw) = (na, nb, nc); CCP4 with MAPC=1/MAPR=2/MAPS=3
    # stores data as (sections=c, rows=b, cols=a) = (nw, nv, nu)
    arr  = np.array(grid, copy=False)
    nu, nv, nw = arr.shape
    data = arr.transpose(2, 1, 0).copy()   # → (nw, nv, nu) = (ns, nr, nc)
    cell = grid.unit_cell
    hdr  = dict(
        nc=nu, nr=nv, ns=nw,
        nx=nu, ny=nv, nz=nw,
        ncstart=0, nrstart=0, nsstart=0,
        mapc=1, mapr=2, maps=3,
        cell=(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma),
        spacegroup=grid.spacegroup.number if grid.spacegroup else 1,
    )
    return data, hdr


def run_refinement(pdb_path, mtz_path, outdir='.', n_cycles=5,
                   fill_asu=False, fill_method='nan',
                   completeness_threshold=0.99):
    """Run refmac5 or phenix.refine on pdb_path+mtz_path; return
    (refined_pdb_path, refined_mtz_path).

    Tries refmac5 first ($CCP4/bin or PATH), then phenix.refine.
    Calls sys.exit(1) if both fail.

    Returning both paths lets the caller chain subsequent operations
    (altindex resolution, peak hunting, etc.) on the refined atom set —
    which can differ from the input atom set when refmac drops
    disordered residues during the run.

    fill_asu : if True, pad input MTZ to full SG-ASU coverage by adding
                 rows for HKLs missing from the deposited dataset before
                 refinement.  Refmac/phenix only emit FWT for HKLs present
                 in input, so an incomplete input means an incomplete
                 bent.mtz downstream.  When False and input is below
                 `completeness_threshold` (default 99%), sys.exit(1)
                 with instructions.
    fill_method : how to populate the filled rows.
        'nan' (default) — HKL + work-set FreeR only; FP/SIGFP left NaN.
                   Refmac skips NaN-Fobs rows in the LSQ target but still
                   emits FC/FWT/2FOFCWT for them from the model on the
                   refined scale.  Scale-safe by construction (no Fobs
                   population added), matches the historical CLAUDE.md
                   behaviour.
        'sigmaa' — see SIGMAA_FC_FILL.md.  Recomputes Fc via
                   `gemmi sfcalc --scale-to` (bulk solvent + aniso
                   scaling onto the Fo scale), estimates per-shell σA
                   and Wilson SigN on FREE==0, fits ln σA + ln SigN
                   linearly in s², then fills each missing HKL with
                   FP = D·|Fc|, SIGFP = √((1-σA²)·SigN).  Also emits
                   diagnostic FC/PHIC columns on every row.  Fixes the
                   raw-DensityCalculatorX ~1/18-scale bug flagged in
                   CLAUDE.md.
    completeness_threshold : minimum SG-ASU completeness required when
                 fill_asu=False.
    """
    import shutil as _sh2
    os.makedirs(outdir, exist_ok=True)
    stem    = os.path.splitext(os.path.basename(pdb_path))[0]
    out_pdb = os.path.join(outdir, f'{stem}_refined.pdb')
    out_mtz = os.path.join(outdir, f'{stem}_refined.mtz')

    # ── Completeness check.  Refmac writes FWT/PHWT only for HKLs present
    # in the input MTZ; an incomplete input means an incomplete FWT and
    # downstream bent.mtz inherits the gaps (Coot shows missing chunks).
    # Below `completeness_threshold` (default 99%) we refuse to proceed
    # unless the caller opts into `fill_asu=True`, which pads the input
    # MTZ to full SG-ASU coverage with HKL+work-set-FreeR rows (FP=NaN).
    # Refmac/phenix then compute FC/FWT/2FOFCWT from the model on those
    # rows during refinement.  Opt-in because silently padding changes
    # the rowcount refmac reports and the FWT coverage downstream tools
    # see — not the refinement target itself, but worth being explicit.
    n_obs, n_exp, frac = _mtz_completeness(mtz_path)
    print(f'  input completeness: {n_obs}/{n_exp} = {100*frac:.1f}%',
          flush=True)
    mtz_path_for_refmac = mtz_path
    if frac < completeness_threshold:
        if not fill_asu:
            print(f'ERROR: {os.path.basename(mtz_path)} is {100*frac:.1f}% '
                  f'complete (< {100*completeness_threshold:.0f}% threshold).\n'
                  f'  Re-run with fill_asu=True (or pass --fill-asu on CLI) '
                  f'to pad the input MTZ to full SG-ASU coverage so refmac '
                  f'writes FWT/FC_ALL for every HKL (added rows carry only HKL '
                  f'+ FreeR; the refinement target sees only deposited Fobs).',
                  file=sys.stderr)
            sys.exit(1)
        mtz_filled = os.path.join(outdir, f'{stem}_filled.mtz')
        try:
            if fill_method == 'sigmaa':
                _sigmaa_fill_mtz(mtz_path, pdb_path, mtz_filled,
                                workdir=outdir, verbose=True)
                # sigmaa fill logs its own per-shell summary; still report
                # the row-count change for parity with the NaN path.
                n_in_rows  = len(gemmi.read_mtz_file(mtz_path).array)
                n_out_rows = len(gemmi.read_mtz_file(mtz_filled).array)
                print(f'  sigmaA-fill added {n_out_rows - n_in_rows} SG-ASU '
                      f'rows with FP=D·|Fc|, SIGFP=√((1-σA²)·SigN) on the '
                      f'Fo scale (deposited Fobs untouched); FC/PHIC '
                      f'diagnostic columns emitted → {n_out_rows} total rows',
                      flush=True)
            elif fill_method == 'nan':
                _pad_mtz_to_full_asu(mtz_path, mtz_filled)
                # Report row-count change, not Fobs completeness: the filled
                # rows carry only HKL + FreeR (FP/SIGFP stay NaN by design —
                # refmac / phenix generate Fcalc/FWT/2FOFCWT for them).
                # _mtz_completeness measures Fobs completeness so it would
                # still print the input's ~90% post-fill and read like the
                # fill silently no-op'd.
                n_in_rows  = len(gemmi.read_mtz_file(mtz_path).array)
                n_out_rows = len(gemmi.read_mtz_file(mtz_filled).array)
                print(f'  padded {n_out_rows - n_in_rows} missing SG-ASU rows '
                      f'with HKL + work-set FreeR only (FP left NaN — refmac '
                      f'generates FC/FWT from the model, Fobs unchanged) → '
                      f'{n_out_rows} total rows in MTZ', flush=True)
            else:
                raise ValueError(
                    f"fill_method must be 'nan' or 'sigmaa', got "
                    f"{fill_method!r}")
            mtz_path_for_refmac = mtz_filled
        except Exception as e:
            print(f'  fill ({fill_method}) failed ({e!r}); aborting',
                  file=sys.stderr)
            sys.exit(1)

    # ── refmac5 ──────────────────────────────────────────────────────────────
    refmac = (_sh2.which('refmac5')
              or os.path.join(os.environ.get('CCP4', ''), 'bin', 'refmac5'))
    if refmac and os.path.exists(refmac):
        mtz_in = gemmi.read_mtz_file(mtz_path_for_refmac)
        cols   = {c.label for c in mtz_in.columns}
        fp     = next((l for l in ('FP', 'F', 'FOBS')              if l in cols), None)
        sigfp  = next((l for l in ('SIGFP', 'SIGF', 'SIGFOBS')     if l in cols), None)
        i_lbl  = next((l for l in ('IMEAN', 'I', 'IOBS')           if l in cols), None)
        sigi   = next((l for l in ('SIGIMEAN', 'SIGI', 'SIGIOBS')  if l in cols), None)
        # I-only MTZs (e.g. RCSB raw deposits straight from cif2mtz) need
        # an F column for REFI TYPE RIGID — refmac's intensity input path
        # doesn't run French-Wilson in rigid mode (reports "All observed
        # amplitudes are 0").  Auto-ctruncate to add F/SIGF, then point
        # refmac at the truncated MTZ.
        if not (fp and sigfp) and (i_lbl and sigi):
            ctruncate = _sh2.which('ctruncate')
            if not ctruncate:
                print('ERROR: input MTZ has I but not F, and `ctruncate` is not '
                      'on PATH.  Source the CCP4 setup, or pre-convert the MTZ.',
                      file=sys.stderr)
                sys.exit(1)
            mtz_F = os.path.join(outdir, f'{stem}_truncated.mtz')
            print(f'  ctruncate {i_lbl}/{sigi} → F/SIGF ({mtz_F})', flush=True)
            r0 = subprocess.run(
                [ctruncate, '-hklin', mtz_path_for_refmac, '-hklout', mtz_F,
                 '-colin', f'/*/*/[{i_lbl},{sigi}]'],
                capture_output=True, text=True)
            log_path = os.path.join(outdir, f'{stem}_ctruncate.log')
            with open(log_path, 'w') as _lf:
                _lf.write(f'# ctruncate returncode: {r0.returncode}\n\n--- stdout ---\n')
                _lf.write(r0.stdout or '')
                _lf.write('\n--- stderr ---\n')
                _lf.write(r0.stderr or '')
            if r0.returncode != 0 or not os.path.exists(mtz_F):
                print(f'  ctruncate failed (rc={r0.returncode}); see {log_path}',
                      file=sys.stderr)
                sys.exit(1)
            mtz_path_for_refmac = mtz_F
            # Re-detect F columns from the truncated MTZ.
            mtz_in = gemmi.read_mtz_file(mtz_path_for_refmac)
            cols   = {c.label for c in mtz_in.columns}
            fp     = next((l for l in ('FP', 'F', 'FOBS')          if l in cols), None)
            sigfp  = next((l for l in ('SIGFP', 'SIGF', 'SIGFOBS') if l in cols), None)
        if fp and sigfp:
            # No FREE= in LABIN — refmac does not require a free set for
            # rigid-body refinement, and propagating the deposit's free
            # flags to LABIN sometimes triggered divergence (e.g. dhfr
            # 1rx2.mtz carries a degenerate FREE column where all observed
            # rows are flagged as work, leaving the "free" set populated
            # only by the asu-pad rows after run_refinement → Rfree=0.91
            # → SigmaA collapse → rigid-body shifts to R=0.53).  We don't
            # use Rfree downstream and refmac handles its absence cleanly.
            # (The phenix.refine fallback below picks free flags directly
            # from the MTZ column types and is unaffected.)
            labin = f'FP={fp} SIGFP={sigfp}'
        else:
            print(f'ERROR: {os.path.basename(mtz_path_for_refmac)} has no usable '
                  f'amplitude columns (looked for FP/F/FOBS+SIGFP/SIGF/SIGFOBS).  '
                  f'Columns present: {sorted(cols)}.', file=sys.stderr)
            sys.exit(1)
        # REFI TYPE RIGID: fast, no atomic refinement, insensitive to
        # unknown ligands (no restraints needed).  We only need refmac to
        # produce clean FWT/PHWT/SIGFP-scaled output; the input PDB is
        # already at the correct geometry.
        script = (f'REFI TYPE RIGID\n'
                  f'NCYC {n_cycles}\n'
                  f'LABIN {labin}\n'
                  f'END\n')
        cmd    = [refmac, 'xyzin', pdb_path, 'xyzout', out_pdb,
                  'hklin', mtz_path_for_refmac, 'hklout', out_mtz]
        print(f'  run_refinement: {os.path.basename(refmac)} rigid '
              f'({n_cycles} cycles)...', flush=True)
        r = subprocess.run(cmd, input=script, capture_output=True, text=True)
        log_path = os.path.join(outdir, f'{stem}_refmac.log')
        with open(log_path, 'w') as _lf:
            _lf.write(f'# refmac5 command: {" ".join(cmd)}\n')
            _lf.write(f'# script:\n{script}\n')
            _lf.write(f'# returncode: {r.returncode}\n\n--- stdout ---\n')
            _lf.write(r.stdout or '')
            _lf.write('\n--- stderr ---\n')
            _lf.write(r.stderr or '')
        if r.returncode == 0 and os.path.exists(out_mtz):
            print(f'  refmac5 done → {out_mtz}  (log: {log_path})', flush=True)
            return out_pdb, out_mtz
        print(f'  refmac5 failed (rc={r.returncode}); see {log_path}',
              file=sys.stderr)

    # ── phenix.refine ─────────────────────────────────────────────────────────
    phenix = _sh2.which('phenix.refine')
    if phenix:
        # phenix runs in cwd=outdir so its scratch files stay there.  All
        # file paths in the cmd must therefore be absolute (or resolved
        # relative to the new cwd), otherwise phenix prints "No files
        # found" and exits with "arguments are not recognized: <basename>".
        abs_pdb = os.path.abspath(pdb_path)
        abs_mtz = os.path.abspath(mtz_path_for_refmac)
        # output.prefix=<stem> makes phenix write <stem>_001.{pdb,mtz}; default
        # prefix would produce <stem>_refine_001.* which clutters downstream
        # filename matching.  We then look for the _001 suffix below.
        cmd = [phenix, abs_pdb, abs_mtz,
               f'refinement.main.number_of_macro_cycles={n_cycles}',
               f'output.prefix={stem}']
        print(f'  run_refinement: phenix.refine ({n_cycles} cycles)...', flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=outdir)
        log_path = os.path.join(outdir, f'{stem}_phenix.log')
        with open(log_path, 'w') as _lf:
            _lf.write(f'# phenix.refine command: {" ".join(cmd)}\n')
            _lf.write(f'# cwd: {outdir}\n')
            _lf.write(f'# returncode: {r.returncode}\n\n--- stdout ---\n')
            _lf.write(r.stdout or '')
            _lf.write('\n--- stderr ---\n')
            _lf.write(r.stderr or '')
        candidate_mtz = os.path.join(outdir, f'{stem}_001.mtz')
        candidate_pdb = os.path.join(outdir, f'{stem}_001.pdb')
        if r.returncode == 0 and os.path.exists(candidate_mtz):
            print(f'  phenix.refine done → {candidate_mtz}  (log: {log_path})',
                  flush=True)
            return candidate_pdb, candidate_mtz
        print(f'  phenix.refine failed (rc={r.returncode}); see {log_path}',
              file=sys.stderr)

    print('ERROR: neither refmac5 nor phenix.refine succeeded.\n'
          'Provide a refined MTZ with FWT/PHWT or 2FOFCWT/PH2FOFCWT.',
          file=sys.stderr)
    sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# Altindex / origin resolution — discrete enumeration ranked by post-transform
# CA RMSD; on a non-trivial hit, transform mov PDB + reindex mov MTZ Fobs and
# re-refine (because reindexed phases are incompatible with map generation —
# lesson from /home/jamesh/projects/fft_symmetry/claude_test/sfcalc_gpu_collapse.py).
# ══════════════════════════════════════════════════════════════════════════════

def _ortho_pair(cell):
    """Return (O, F) where O is fractional→cartesian and F is cartesian→fractional
    column-vector matrices for the given gemmi.UnitCell."""
    O = np.zeros((3, 3))
    for i, b in enumerate([(1,0,0), (0,1,0), (0,0,1)]):
        p = cell.orthogonalize(gemmi.Fractional(*b))
        O[:, i] = [p.x, p.y, p.z]
    return O, np.linalg.inv(O)


def _collect_matched_ca(mov_st, ref_st):
    """Return (A, B) cartesian arrays of CA atoms matched by (chain, resi, name).
    Drops residues with multiple altlocs."""
    def _ca_dict(st):
        d = {}
        for model in st:
            for chain in model:
                for res in chain:
                    for atom in res:
                        if atom.name.strip() != 'CA':
                            continue
                        key = (chain.name, res.seqid.num, atom.name.strip())
                        if key in d:
                            continue   # skip altlocs / duplicates
                        d[key] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
            break   # first model only
        return d
    a = _ca_dict(mov_st)
    b = _ca_dict(ref_st)
    common = sorted(set(a) & set(b))
    if not common:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return (np.array([a[k] for k in common]),
            np.array([b[k] for k in common]))


def _kabsch(A, B):
    """Optimal R, t such that R @ A.T + t best fits B.T.  Returns (R, t, rmsd)."""
    cA = A.mean(axis=0); cB = B.mean(axis=0)
    H  = (A - cA).T @ (B - cB)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = cB - R @ cA
    rmsd = float(np.sqrt(np.mean(np.sum(((A @ R.T) + t - B) ** 2, axis=1))))
    return R, t, rmsd


def _is_metric_preserving(R_frac, G, atol=1e-3):
    """True iff R_frac preserves the lattice metric tensor G."""
    return np.allclose(R_frac.T @ G @ R_frac, G,
                       atol=atol * np.max(np.abs(G)))


def _rot_deviation_deg(R_a, R_b):
    """Rotation angle (deg) of R_a · R_bᵀ.  Both must be cartesian-orthogonal."""
    cos_t = max(-1.0, min(1.0, (np.trace(R_a @ R_b.T) - 1.0) / 2.0))
    return float(np.degrees(np.arccos(cos_t)))


def _shift_mtz_origin(mtz_in, t_frac, mtz_out):
    """Apply F'(H) = exp(+2πi H·t) F(H) to every (F, PHI) column pair in
    mtz_in — the structure-factor transformation for a real-space shift
    of atoms by +t_cart (matches VALIDATE_ALTINDEX.md and pairs correctly
    with _apply_cart_transform_to_pdb(R=I, t_cart), which moves atoms to
    p' = p + t_cart).  Previous versions used -2πi h·t, which moved the
    map by -t_cart while atoms moved by +t_cart, mis-aligning display
    PDBs from their unbent maps by 2·t_cart (typically ~1-3 Å in the
    CA production).  All other columns (intensities, sigF, FREE flag)
    pass through unchanged.  Pure phase rotation — no re-refinement."""
    mtz = gemmi.read_mtz_file(mtz_in)
    miller = mtz.make_miller_array().astype(np.float64)        # (N, 3)
    delta_phi_deg = (360.0 * (miller @ np.asarray(t_frac, dtype=np.float64))
                     ).astype(np.float32)

    # Find (F, PHI) pairs by column type
    cols = list(mtz.columns)
    f_cols   = [(i, c) for i, c in enumerate(cols) if c.type in ('F', 'G')]
    phi_cols = [(i, c) for i, c in enumerate(cols) if c.type in ('P',)]

    data = np.array(mtz.array, copy=True)
    for fi, fc in f_cols:
        # Match phase column by name suffix convention: F<X>/PH<X>, FWT/PHWT, etc.
        # Try common pairings; otherwise pair by next phase column index > fi.
        target = None
        candidates = ['PH' + fc.label, 'PHI' + fc.label,
                      fc.label.replace('F', 'PH', 1),
                      fc.label.replace('FWT', 'PHWT', 1),
                      fc.label.replace('2FOFCWT', 'PH2FOFCWT', 1),
                      fc.label.replace('DELFWT', 'PHDELWT', 1)]
        for pi, pc in phi_cols:
            if pc.label in candidates:
                target = pi; break
        if target is None and len(phi_cols) == 1:
            target = phi_cols[0][0]
        if target is None:
            continue   # no matching phase column — leave alone
        # NaN-safe phase update: leave NaN rows alone, shift the rest.
        # NaN typically appears in missing-data rows (e.g. SIGFP=NaN in
        # incomplete deposits); applying `%` to NaN raises a benign
        # RuntimeWarning that clutters scan logs.
        shifted = data[:, target] + delta_phi_deg
        finite = np.isfinite(shifted)
        out = data[:, target].copy()
        out[finite] = (shifted[finite] % 360.0).astype(np.float32)
        data[:, target] = out
    mtz.set_data(data)
    mtz.write_to_file(mtz_out)


def _reindex_mtz_fobs(mtz_in, R_frac_int, mtz_out):
    """Reindex HKLs by H_new = R_frac_int @ H_old (col form) and drop all map-
    coefficient columns.  Keeps only Fobs/sigF/intensities/free flag — refmac
    will rebuild map coeffs against the rotated PDB.

    R_frac_int must be an integer 3x3 with det = ±1."""
    R = np.asarray(R_frac_int, dtype=int)
    R_inv = np.round(np.linalg.inv(R.astype(float))).astype(int)
    if not np.allclose(R @ R_inv, np.eye(3), atol=1e-6):
        raise ValueError(f"R_frac_int not invertible to integer: {R}")

    mtz = gemmi.read_mtz_file(mtz_in)
    old_miller = mtz.make_miller_array().astype(int)
    # H_new = R H_old  (col form);  row form: H_new_row = H_old_row @ R^T
    new_miller = old_miller @ R.T

    # Decide which columns to keep: HKL + Fobs/sigF/intensities/free flag.
    # Drop map coeffs (FWT/PHWT/2FOFCWT/PH2FOFCWT/DELFWT/PHDELWT/etc.)
    DROP = {'FWT', 'PHWT', '2FOFCWT', 'PH2FOFCWT', 'DELFWT', 'PHDELWT',
            'FC', 'PHIC', 'FC_ALL', 'PHIC_ALL', 'FC_ALL_LS', 'PHIC_ALL_LS',
            'FOM', 'FOM_LS', 'HLA', 'HLB', 'HLC', 'HLD'}
    keep_idx = [0, 1, 2]
    cols = list(mtz.columns)
    for i in range(3, len(cols)):
        if cols[i].label not in DROP:
            keep_idx.append(i)

    old_data = np.array(mtz.array, copy=True)
    new_data = old_data[:, keep_idx]
    new_data[:, 0:3] = new_miller.astype(np.float32)

    out = gemmi.Mtz(with_base=True)
    out.cell = mtz.cell
    out.spacegroup = mtz.spacegroup
    out.add_dataset('reindexed')
    for i in keep_idx[3:]:
        out.add_column(cols[i].label, cols[i].type)
    out.set_data(new_data)
    out.write_to_file(mtz_out)


def _apply_op_to_mtz(mtz_in, R_frac_int, t_frac, mtz_out):
    """Apply the SF transformation theorem to ALL F/PHI column pairs:

        F'(h) = F(R^T · h) · exp(2πi · h · t)

    i.e. HKL relabel h → (R^T)^-1 · h plus a per-HKL phase shift.  For
    SG-symmetric inputs (FWT/PHWT from refmac), F(R^T h) = F(h) (when R
    is in the SG), so the effect on amplitudes is a permutation of
    identical values; phases additionally shift by exp(2πi h · (t − R^T t_sym))
    where t_sym is the SG translation matching R.  This helper applies
    the *full* (R, t) formula; the caller is responsible for picking the
    right t.

    The formula matches VALIDATE_ALTINDEX.md (validated end-to-end through
    refmac on insulin H 3 to 1.24° mean Δphi).

    R_frac_int : integer 3×3 with det = ±1 (the rotation in fractional coords)
    t_frac     : (3,) fractional translation
    """
    R = np.asarray(R_frac_int, dtype=int)
    R_inv = np.round(np.linalg.inv(R.astype(float))).astype(int)
    if not np.allclose(R @ R_inv, np.eye(3), atol=1e-6):
        raise ValueError(f"R not invertible to integer: {R}")

    mtz = gemmi.read_mtz_file(mtz_in)
    old_miller = mtz.make_miller_array().astype(int)
    # F_new(h_new) = F_orig(R^T h_new); equivalently, an input row at
    # h_orig contributes to h_new = (R^T)^-1 h_orig = R^-T h_orig.
    # Row form: h_new_row = h_orig_row · R^-1.
    new_miller = old_miller @ R_inv

    # Phase shift per row: PHI' = PHI + 360° · h_new · t
    delta_phi_deg = (360.0 * (new_miller.astype(np.float64)
                              @ np.asarray(t_frac, dtype=np.float64))
                     ).astype(np.float32)

    data = np.array(mtz.array, copy=True)
    data[:, 0:3] = new_miller.astype(np.float32)

    cols = list(mtz.columns)
    f_cols = [(i, c) for i, c in enumerate(cols) if c.type in ('F', 'G')]
    phi_cols = [(i, c) for i, c in enumerate(cols) if c.type == 'P']

    for fi, fc in f_cols:
        # Match phase by canonical name suffix (FWT↔PHWT, F↔PHI, etc.)
        target = None
        candidates = ['PH' + fc.label, 'PHI' + fc.label,
                      fc.label.replace('FWT', 'PHWT', 1),
                      fc.label.replace('2FOFCWT', 'PH2FOFCWT', 1),
                      fc.label.replace('DELFWT', 'PHDELWT', 1),
                      fc.label.replace('FC', 'PHIC', 1)]
        for pi, pc in phi_cols:
            if pc.label in candidates:
                target = pi; break
        if target is None and len(phi_cols) == 1:
            target = phi_cols[0][0]
        if target is None:
            continue
        data[:, target] = ((data[:, target] + delta_phi_deg) % 360.0
                            ).astype(np.float32)

    mtz.set_data(data)
    mtz.write_to_file(mtz_out)


def _apply_cart_transform_to_pdb(pdb_in, R_cart, t_cart, pdb_out, ref_cell=None):
    """Write pdb_out with all atoms transformed by p' = R_cart @ p + t_cart."""
    st = gemmi.read_pdb(pdb_in)
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    p = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    q = R_cart @ p + np.asarray(t_cart, dtype=float)
                    atom.pos = gemmi.Position(float(q[0]), float(q[1]), float(q[2]))
    if ref_cell is not None:
        st.cell = ref_cell
    st.write_pdb(pdb_out)


def _enum_alt_rot_origin_candidates(cell, sg, sg_name, metric_tol_rel=1e-6):
    """Yield (R_frac, t_frac) pairs covering (altindex × symop × origin) candidates
    for resolve_altindex.  R_frac is fractional, t_frac is fractional mod 1.

    sg_name is the H-M symbol from the PDB header (e.g. 'H 3'); used to look up
    origin choices in _ORIGINS_TABLE.  gemmi's xhm() returns 'R 3:H' for H3 etc.,
    which doesn't match _ORIGINS_TABLE keys.

    Altindex rotations come from `_get_altindex_ops(sg_name, cell, metric_tol_rel)` —
    metric-tensor enumeration of integer proper rotations preserving the
    cell's lattice that are NOT already in the SG.  metric_tol_rel controls
    how strictly the lattice metric must be preserved; loosening it (e.g.
    0.05) admits ops that hold for a NEARBY higher-symmetry metric, useful
    for pseudo-symmetric cells where the natural alt-cell ops just miss
    strict tolerance (see _get_altindex_ops docstring)."""
    alt_ops = [(np.eye(3), np.zeros(3))] + list(
        _get_altindex_ops(sg_name, cell=cell, metric_tol_rel=metric_tol_rel))
    sym_ops = list(sg.operations())
    key = _sg_origins_key(sg_name)
    entries = _ORIGINS_TABLE.get(key)
    if entries is None:
        # gemmi-name fallback (R 3:H, R 3 2:H, etc.)
        entries = _ORIGINS_TABLE.get(_sg_origins_key(sg.xhm()))
    if entries is None:
        origins = [(dx, dy, dz)
                   for dx in (0, .5) for dy in (0, .5) for dz in (0, .5)]
    else:
        origins = _expand_origin_entries(entries, n_polar=12)
    origins = np.array(origins, dtype=float)

    for alt_R, alt_t in alt_ops:
        for op in sym_ops:
            sym_R = np.array(op.rot, dtype=float) / op.DEN
            sym_t = np.array(op.tran, dtype=float) / op.DEN
            R_comb = sym_R @ alt_R
            t_base = sym_R @ alt_t + sym_t   # carry altindex t through sym
            for o in origins:
                yield R_comb, (t_base + o)


def _no_op_resolve_result(mov_pdb, mov_mtz, drot=0.0, rmsd=0.0):
    return dict(mov_pdb_out=mov_pdb, mov_mtz_out=mov_mtz,
                action='none', R_frac=np.eye(3), t_frac=np.zeros(3),
                drot_deg=drot, rmsd_after=rmsd)


def _cells_match(c1, c2, rtol=1e-4):
    """True iff two gemmi.UnitCell objects agree to relative tolerance rtol
    on all six parameters (a, b, c, α, β, γ)."""
    p1 = np.array(c1.parameters, dtype=float)
    p2 = np.array(c2.parameters, dtype=float)
    return np.allclose(p1, p2, rtol=rtol, atol=0)


def _one_letter_for_resname(name):
    """Return uppercase one-letter code for an amino-acid residue name, or
    None for HETATMs / waters / metals / unknowns.  Modified AAs like MSE
    map to their parent residue's letter (M)."""
    info = gemmi.find_tabulated_residue(name)
    if info is None or not info.is_amino_acid():
        return None
    letter = info.one_letter_code
    if letter is None or letter == ' ':
        return None
    return letter.upper()


def _align_polymer_by_ca(mov_st, ref_st, chain_remap=None):
    """Pairwise-align the mov and ref amino-acid CA sequences (first model,
    single chain each — after applying chain_remap if given) via gemmi's
    Needleman-Wunsch aligner.  Return a mapping and quality metrics.

    Handles constant offsets AND internal INDELs / N-terminal truncations
    that stump the fixed-offset detector (e.g. 6yzt/6yzv/8phm vs 6klz).

    Parameters
    ----------
    mov_st, ref_st : gemmi.Structure
    chain_remap : tuple(str, str) or None
        If given, treat mov chain `chain_remap[0]` as if it had ID
        `chain_remap[1]` (matches stage-1 rename semantics).

    Returns
    -------
    dict or None — dict has keys:
        residue_map  : {mov_resnum: ref_resnum} for matched columns only
                       (mismatches / gaps skipped).
        n_aligned    : alignment columns where both sequences have a residue
                       (whether or not they match).
        n_matches    : aligned columns whose letters match.
        identity     : n_matches / n_aligned in [0, 1].
        mov_coverage : n_aligned / (# CA residues in mov chain), in [0, 1].
        ref_coverage : n_aligned / (# CA residues in ref chain), in [0, 1].
    None if either structure has fewer than 10 mappable CA residues or
    is not single-chain.
    """
    def _ca_seq(st, target_chain=None, rename_from=None):
        # Return parallel lists of (one_letter, resnum) for the first
        # model's target_chain (or first chain if None).  Non-amino-acids
        # (HETATMs, waters) skipped; unknown letters become 'X'.
        out_letters = []
        out_resnums = []
        for model in st:
            for chain in model:
                cname = chain.name
                if rename_from is not None and cname == rename_from[0]:
                    cname = rename_from[1]
                if target_chain is not None and cname != target_chain:
                    continue
                for res in chain:
                    letter = _one_letter_for_resname(res.name)
                    if letter is None:
                        continue
                    out_letters.append(letter)
                    out_resnums.append(res.seqid.num)
                return out_letters, out_resnums
            break
        return out_letters, out_resnums

    mov_chains = sorted({c.name for m in mov_st for c in m
                         for r in c if any(a.name.strip() == 'CA' for a in r)})
    ref_chains = sorted({c.name for m in ref_st for c in m
                         for r in c if any(a.name.strip() == 'CA' for a in r)})
    if len(mov_chains) != 1 or len(ref_chains) != 1:
        return None
    mov_chain0 = mov_chains[0]
    ref_chain0 = ref_chains[0]
    # Under a chain_remap, mov's chain is treated as if renamed; the ref
    # chain we want to align against is chain_remap[1] which should equal
    # ref_chain0.
    if chain_remap is not None:
        if chain_remap[0] != mov_chain0 or chain_remap[1] != ref_chain0:
            return None

    mov_seq, mov_resnums = _ca_seq(mov_st, target_chain=mov_chain0)
    ref_seq, ref_resnums = _ca_seq(ref_st, target_chain=ref_chain0)
    if len(mov_seq) < 10 or len(ref_seq) < 10:
        return None

    scoring = gemmi.AlignmentScoring()
    result  = gemmi.align_string_sequences(list(mov_seq), list(ref_seq),
                                           [], scoring)
    mov_g = result.add_gaps(''.join(mov_seq), 1)
    ref_g = result.add_gaps(''.join(ref_seq), 2)
    if len(mov_g) != len(ref_g):
        return None

    # Map EVERY aligned column (both residues present), regardless of whether
    # the letters match.  Match-atoms downstream still filters by resname via
    # the uid — so mismatched residues won't erroneously pair up — but giving
    # every mov residue in the alignment a definite ref-side resnum avoids
    # collisions between "unmapped mov residues keeping their original resnum"
    # and "mapped mov residues taking that resnum from ref".
    residue_map = {}
    n_aligned = 0
    n_matches = 0
    mov_i = 0
    ref_i = 0
    for ch_m, ch_r in zip(mov_g, ref_g):
        m_present = (ch_m != '-')
        r_present = (ch_r != '-')
        if m_present and r_present:
            n_aligned += 1
            if ch_m == ch_r:
                n_matches += 1
            residue_map[mov_resnums[mov_i]] = ref_resnums[ref_i]
        if m_present:
            mov_i += 1
        if r_present:
            ref_i += 1

    if n_aligned == 0:
        return None
    return dict(
        residue_map=residue_map,
        n_aligned=n_aligned,
        n_matches=n_matches,
        identity=n_matches / n_aligned,
        mov_coverage=n_aligned / len(mov_seq),
        ref_coverage=n_aligned / len(ref_seq),
    )


def _detect_chain_resnum_remap(mov_st, ref_st):
    """Detect a chain rename and/or resnum-offset (constant OR per-residue
    via sequence alignment) that maximizes atom-match uid overlap between
    two structures.

    Three-stage detection with different stringency for each step:

    1. **Chain rename** — if both structures are single-chain with
       different chain IDs, apply the rename if it brings the offset=0
       residue-name match count from <30 to ≥30 with ≥50% identity.
       Lenient threshold (50%) admits same-protein cases where an
       internal indel shifts most resnums by 1 (like 5jdv) — these
       look like 50% identity under constant-offset matching but are
       ~95% identical under proper alignment.

    2. **Constant resnum offset (fast path)** — try non-zero constant
       offsets and apply the best one only if (a) it adds ≥30 matches
       over offset=0 AND (b) the per-offset frac is ≥90%.  Strict
       threshold (90%) avoids false-positive offset shifts on
       related-but-different isoforms where a coincidental offset
       gives 55% identity.

    3. **Per-residue sequence alignment (fallback)** — if stage 2
       rejected, run gemmi Needleman-Wunsch on the CA sequences to
       find a per-residue map that handles internal INDELs (e.g.
       6yzt/6yzv/8phm vs 6klz, all human CA II but with mixed-offset
       numbering due to a missing residue in ref).  Accept if
       identity ≥ 0.90, mov_coverage ≥ 0.50 AND n_matches ≥ 30.

    Returns dict with chain_remap, resnum_offset, residue_map (or None),
    matches, overlap, frac, or None if no justified remap.

    Pure detection — does not mutate inputs.
    """
    def _ca_residues(st):
        out = []
        for model in st:
            for chain in model:
                for res in chain:
                    if any(a.name.strip() == 'CA' for a in res):
                        out.append((chain.name, res.seqid.num, res.name))
            break
        return out

    a = _ca_residues(mov_st)
    b = _ca_residues(ref_st)
    if len(a) < 10 or len(b) < 10:
        return None

    chains_a = sorted(set(c for c, _, _ in a))
    chains_b = sorted(set(c for c, _, _ in b))

    b_by_chain = {}
    for c, r, n in b:
        b_by_chain.setdefault(c, {})[r] = n

    def score(a_renamed, offset):
        overlap = matches = 0
        for c, r, n in a_renamed:
            ref_n = b_by_chain.get(c, {}).get(r + offset)
            if ref_n is not None:
                overlap += 1
                if ref_n == n:
                    matches += 1
        return overlap, matches

    a_resnums = [r for _, r, _ in a]
    b_resnums = [r for _, r, _ in b]
    lo = min(b_resnums) - max(a_resnums) - 5
    hi = max(b_resnums) - min(a_resnums) + 5

    # Stage 1: chain rename (lenient — accepts indel-shifted same proteins)
    chain_remap = None
    base_a = a
    overlap0, matches0 = score(base_a, 0)
    if (len(chains_a) == 1 and len(chains_b) == 1 and chains_a != chains_b
            and matches0 < 30):
        cmap = (chains_a[0], chains_b[0])
        cand_a = [(cmap[1] if c == cmap[0] else c, r, n) for c, r, n in a]
        o_r, m_r = score(cand_a, 0)
        if m_r >= 30 and o_r >= 10 and m_r / o_r >= 0.5:
            chain_remap = cmap
            base_a = cand_a
            overlap0, matches0 = o_r, m_r

    # Stage 2: resnum offset (strict — must be near-perfect)
    best_offset = 0
    best_overlap, best_matches = overlap0, matches0
    for offset in range(lo, hi + 1):
        if offset == 0:
            continue
        o, m = score(base_a, offset)
        if o < 10:
            continue
        # Strict: need ≥30 improvement AND ≥90% frac at the new offset
        if m >= matches0 + 30 and m / o >= 0.9 and m > best_matches:
            best_offset = offset
            best_overlap = o
            best_matches = m

    stage2_frac = (best_matches / best_overlap) if best_overlap else 0.0

    # Stage 3: sequence-alignment fallback — only when stage 2 didn't
    # accept a non-zero offset AND fewer than 90% of mov's residues
    # already line up at offset=0.  Skips the aligner (cheap but not
    # free) for pairs that share identical numbering (lyso, dhfr,
    # magdoff, most of CA production) where the identity map is
    # already correct.
    residue_map = None
    stage3_metrics = None
    frac0_of_mov = (matches0 / len(a)) if len(a) else 0.0
    if best_offset == 0 and frac0_of_mov < 0.90:
        # Stage 2 didn't help.  Try the aligner (may return None for
        # multi-chain, too-short, or degenerate cases).
        al = _align_polymer_by_ca(mov_st, ref_st, chain_remap=chain_remap)
        if (al is not None
                and al['n_matches'] >= 30
                and al['identity'] >= 0.90
                and al['mov_coverage'] >= 0.50):
            residue_map  = al['residue_map']
            best_matches = al['n_matches']
            best_overlap = al['n_aligned']
            stage2_frac  = al['identity']
            stage3_metrics = al

    if chain_remap is None and best_offset == 0 and residue_map is None:
        return None
    return dict(chain_remap=chain_remap, resnum_offset=best_offset,
                residue_map=residue_map,
                matches=best_matches, overlap=best_overlap,
                frac=stage2_frac,
                mov_coverage=(stage3_metrics or {}).get('mov_coverage'))


def _normalize_mov_pdb(mov_pdb, ref_pdb, out_pdb, verbose=True):
    """Apply chain rename, constant resnum offset, or per-residue sequence-
    aligned resnum map to mov_pdb so its (chain, resnum) keys align with
    ref_pdb.  Writes corrected PDB to out_pdb only if a justified remap is
    detected; returns (out_pdb, plan) in that case, else (mov_pdb, None).

    `plan` is the detector's dict (chain_remap, resnum_offset, residue_map,
    ...) — usable by `_reverse_normalize_pdb` later to restore mov's
    original labels on downstream output PDBs.

    Handles three flavours of nomenclature drift:
      - old PDBs with +1000 resnum offsets (1zfk, 3v7x, 3vbd)
      - single-chain entries with B/X chain labels (5j**)
      - same-protein cases where the constant offset drifts across an
        internal INDEL (6yzt/6yzv/8phm — same human CA II as 6klz)"""
    mov_st = gemmi.read_pdb(mov_pdb)
    ref_st = gemmi.read_pdb(ref_pdb)
    plan = _detect_chain_resnum_remap(mov_st, ref_st)
    if plan is None:
        return mov_pdb, None
    cmap  = plan['chain_remap']
    off   = plan['resnum_offset']
    rmap  = plan.get('residue_map')
    if verbose:
        msg = []
        if cmap is not None:
            msg.append(f"chain '{cmap[0]}' → '{cmap[1]}'")
        if off != 0:
            msg.append(f"resnum offset {off:+d}")
        if rmap is not None:
            mov_cov = plan.get('mov_coverage')
            cov_str = f", mov_cov {100*mov_cov:.0f}%" if mov_cov is not None else ""
            msg.append(f"residue-map ({len(rmap)} residues){cov_str}")
        print(f"  _normalize_mov_pdb: {', '.join(msg)} "
              f"({plan['matches']}/{plan['overlap']} CA-residues "
              f"name-match, {100*plan['frac']:.0f}% sequence id)",
              flush=True)
    n_dropped = 0
    for model in mov_st:
        for chain in model:
            if cmap is not None and chain.name == cmap[0]:
                chain.name = cmap[1]
            if rmap is not None:
                # Drop unmapped protein residues first — they'd collide with
                # a mapped residue's target resnum (e.g. mov has an extra
                # N-terminal insertion residue: {3:4, 5:5, ...} leaves mov's
                # original res 4 unmapped, and it would keep num=4 while
                # mov res 3 also lands at num=4).  Unmapped protein
                # residues have no ref counterpart per the alignment, so
                # they can't participate in the fit anyway.  HETATMs
                # (waters, ligands, metals — no CA) are always kept.
                to_drop = []
                for res in chain:
                    if res.seqid.num in rmap:
                        continue
                    has_ca = any(a.name == 'CA' for a in res)
                    if has_ca:
                        to_drop.append(res.seqid)
                for sid in to_drop:
                    # gemmi: remove by seqid; iterate to find + del
                    for i in range(len(chain) - 1, -1, -1):
                        if chain[i].seqid == sid:
                            del chain[i]
                            n_dropped += 1
                            break
                for res in chain:
                    new_num = rmap.get(res.seqid.num)
                    if new_num is not None:
                        res.seqid.num = new_num
            elif off != 0:
                for res in chain:
                    res.seqid.num += off
        break
    if verbose and n_dropped:
        print(f"  _normalize_mov_pdb: dropped {n_dropped} unmapped "
              f"protein residues (insertions relative to ref)", flush=True)
    mov_st.write_pdb(out_pdb)
    return out_pdb, plan


def _reverse_normalize_pdb(pdb_in, plan, pdb_out):
    """Inverse of `_normalize_mov_pdb`: given a PDB whose labels are on
    ref's side (because it descended from a normalized mov file — e.g.
    the resolve_altindex-stretched/refined mov, or any bent PDB written
    from it) and the `plan` dict returned by `_normalize_mov_pdb`, write
    `pdb_out` with mov's ORIGINAL chain / resnum labels restored.
    Coordinates, residue names, atom names, and HETATMs are untouched.

    Used by fitreso_scan under `preserve_mov_numbering=True` (default)
    so that the user-facing `bent.pdb` / `unbent.pdb` retain the naming
    conventions of the moving structure they handed in."""
    if plan is None:
        _shutil.copy(pdb_in, pdb_out)
        return pdb_out
    st = gemmi.read_pdb(pdb_in)
    cmap = plan.get('chain_remap')
    off  = plan.get('resnum_offset', 0)
    rmap = plan.get('residue_map')
    inv_rmap = None
    if rmap is not None:
        inv_rmap = {v: k for k, v in rmap.items()}
    for model in st:
        for chain in model:
            # Reverse per-residue map or constant offset (never both —
            # detector always produces exactly one of the two).
            if inv_rmap is not None:
                for res in chain:
                    orig_num = inv_rmap.get(res.seqid.num)
                    if orig_num is not None:
                        res.seqid.num = orig_num
            elif off != 0:
                for res in chain:
                    res.seqid.num -= off
            # Reverse chain rename LAST (so we don't disturb the
            # chain-name match used by the residue loop above).
            if cmap is not None and chain.name == cmap[1]:
                chain.name = cmap[0]
        break
    st.write_pdb(pdb_out)
    return pdb_out


def _stretch_pdb_to_cell(in_pdb, target_cell, out_pdb):
    """Re-orthogonalize a PDB into target_cell, preserving each atom's
    fractional coordinates.  Port of origins/claude/stretch_to_cell.py.

    Use this when comparing two non-isomorphous crystals: putting the moving
    structure into the reference's cell with fractions preserved absorbs the
    isomorphous cell-distortion shift as a uniform elastic stretch, so the
    remaining real-space difference is just the alt-indexing rotation +
    translation that the discrete enumeration knows how to find."""
    st = gemmi.read_pdb(in_pdb)
    oc = st.cell
    old_cell = gemmi.UnitCell(oc.a, oc.b, oc.c, oc.alpha, oc.beta, oc.gamma)
    st.cell = target_cell
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    frac = old_cell.fractionalize(atom.pos)
                    pos  = target_cell.orthogonalize(frac)
                    atom.pos = gemmi.Position(pos.x, pos.y, pos.z)
        break
    st.write_pdb(out_pdb)


def _relabel_mtz_cell(in_mtz, target_cell, out_mtz):
    """Copy in_mtz to out_mtz with cell (file-level AND per-dataset) replaced
    by target_cell.  HKLs and column values are unchanged.

    The per-dataset cell must also be updated — CCP4 reads the dataset cell,
    not the file-level cell, and a mismatch triggers a refmac "Large
    differences between cells from pdb and mtz" abort even when mtz.cell
    matches pdb.cell."""
    mtz = gemmi.read_mtz_file(in_mtz)
    mtz.cell = target_cell
    for ds in mtz.datasets:
        ds.cell = target_cell
    mtz.write_to_file(out_mtz)


def _finish_after_stretch_only(mov_pdb, mov_mtz, outdir, refine_cycles,
                               fill_asu, fill_method, drot_deg, rmsd,
                               verbose):
    """When cell-stretching alone aligned moving to ref (no discrete altindex
    needed), re-refine the stretched moving model+data under the new cell so
    FWT/PHWT are consistent with ref's cell metric.  The pre-stretch FWT/PHWT
    were computed under moving's old cell; resolution-dependent weights
    (sigma_a, bulk solvent) would otherwise be ~cell-mismatch off."""
    if verbose:
        print(f"  cell-stretch only — re-refining stretched moving model "
              f"under ref cell ...", flush=True)
    out_pdb, out_mtz = run_refinement(mov_pdb, mov_mtz, outdir=outdir,
                                       n_cycles=refine_cycles,
                                       fill_asu=fill_asu,
                                       fill_method=fill_method)
    return dict(mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
                action='cell_stretch_only', R_frac=np.eye(3),
                t_frac=np.zeros(3), drot_deg=drot_deg, rmsd_after=rmsd)


def resolve_altindex(mov_pdb, ref_pdb, mov_mtz, outdir,
                     refine_cycles=5, improve_threshold=0.7,
                     verbose=True, fill_asu=False, fill_method='nan'):
    """Find (altindex × symop × origin) discrete transform that aligns mov onto ref.

    Strategy: Kabsch LSQ fit gives the continuous rigid-body answer; enumerate
    discrete (R_frac, t_frac) candidates and rank by post-transform CA RMSD
    (with fractional wrap-to-nearest-image).  Apply the best discrete op:
      - action='altindex_refine' if R != I: transform PDB cartesian, reindex MTZ
        Fobs (drop map coeffs), re-refine to regenerate FWT/PHWT.
      - action='origin_only' if R = I but ||t|| > 0.01: translate PDB, apply
        phase factor exp(-2πi H·t) directly to MTZ — no re-refinement needed.
      - action='none' otherwise.

    Returns dict with mov_pdb_out, mov_mtz_out, action, R_frac, t_frac,
    drot_deg, rmsd_after.

    Writes ``<outdir>/resolve_altindex.log`` with the full trace + final
    action / R / t / drot / rmsd_after, regardless of the action taken.
    """
    os.makedirs(outdir, exist_ok=True)

    # ── Log capture: tee every user-facing line into a buffer that is
    # dumped to `outdir/resolve_altindex.log` at every return path. ─────
    _log_lines = []
    def _log(msg):
        _log_lines.append(str(msg))
        if verbose:
            print(msg, flush=True)
    def _write_log_and_return(result):
        try:
            with open(os.path.join(outdir, 'resolve_altindex.log'), 'w') as fh:
                fh.write(f"# resolve_altindex log\n")
                fh.write(f"# mov_pdb : {mov_pdb}\n")
                fh.write(f"# ref_pdb : {ref_pdb}\n")
                fh.write(f"# mov_mtz : {mov_mtz}\n")
                fh.write(f"# outdir  : {outdir}\n")
                fh.write("\n".join(_log_lines) + "\n\n")
                fh.write(f"# ── Result ────────────────────────────────\n")
                for k in ('action', 'drot_deg', 'rmsd_after'):
                    v = result.get(k)
                    if v is not None:
                        fh.write(f"# {k:<12} = {v}\n")
                Rf = result.get('R_frac')
                tf = result.get('t_frac')
                if Rf is not None:
                    fh.write(f"# R_frac       =\n"
                             f"#   [[{Rf[0][0]:+.6f} {Rf[0][1]:+.6f} {Rf[0][2]:+.6f}]\n"
                             f"#    [{Rf[1][0]:+.6f} {Rf[1][1]:+.6f} {Rf[1][2]:+.6f}]\n"
                             f"#    [{Rf[2][0]:+.6f} {Rf[2][1]:+.6f} {Rf[2][2]:+.6f}]]\n")
                if tf is not None:
                    fh.write(f"# t_frac       = [{tf[0]:+.6f} {tf[1]:+.6f} {tf[2]:+.6f}]\n")
                fh.write(f"# mov_pdb_out  = {result.get('mov_pdb_out')}\n")
                fh.write(f"# mov_mtz_out  = {result.get('mov_mtz_out')}\n")
        except Exception as _e:
            if verbose:
                print(f"  (warning: failed to write resolve_altindex.log: {_e})",
                      flush=True)
        return result

    mov_st = gemmi.read_pdb(mov_pdb)
    ref_st = gemmi.read_pdb(ref_pdb)
    mov_st.setup_entities()
    ref_st.setup_entities()

    # ── cell-stretch pre-alignment for non-isomorphous pairs ─────────────────
    # When the two crystals have different cells (same SG, isomorphous
    # expansion/contraction, etc.), the discrete (R_frac, t_frac) enumeration
    # below would only ever find lattice translations as the "best" op — no
    # purely-fractional rotation can compensate for a metric change.  Stretch
    # the moving model+data into ref's cell first, preserving fractional
    # coordinates; then the rest of the function runs in matched cells.
    stretched = False
    if not _cells_match(mov_st.cell, ref_st.cell):
        s_pdb_stem = os.path.splitext(os.path.basename(mov_pdb))[0]
        s_mtz_stem = os.path.splitext(os.path.basename(mov_mtz))[0]
        s_pdb = os.path.join(outdir, f'{s_pdb_stem}_stretched.pdb')
        s_mtz = os.path.join(outdir, f'{s_mtz_stem}_stretched.mtz')
        _log(f"resolve_altindex: cells differ — stretching moving into "
             f"ref cell (fractions preserved)")
        _log(f"  moving cell:    {tuple(round(x, 3) for x in mov_st.cell.parameters)}")
        _log(f"  ref    cell:    {tuple(round(x, 3) for x in ref_st.cell.parameters)}")
        _stretch_pdb_to_cell(mov_pdb, ref_st.cell, s_pdb)
        _relabel_mtz_cell(mov_mtz, ref_st.cell, s_mtz)
        mov_pdb = s_pdb
        mov_mtz = s_mtz
        mov_st = gemmi.read_pdb(mov_pdb)
        mov_st.setup_entities()
        stretched = True

    A_cart, B_cart = _collect_matched_ca(mov_st, ref_st)
    if len(A_cart) < 3:
        _log(f"resolve_altindex: <3 CA matches; skipping")
        return _write_log_and_return(_no_op_resolve_result(mov_pdb, mov_mtz))

    rmsd_base = float(np.sqrt(np.mean(np.sum((A_cart - B_cart) ** 2, axis=1))))
    R_lsq, t_lsq, rmsd_lsq = _kabsch(A_cart, B_cart)

    _log(f"resolve_altindex: {len(A_cart)} CA pairs  "
         f"baseline RMSD {rmsd_base:.2f} A  LSQ {rmsd_lsq:.2f} A")

    # ── cross-cell: relax metric tolerance for the altindex enumeration ─────
    # For cross-cell pairs (stretched=True), the natural alt-cell op is
    # often metric-preserving in a HIGHER-symmetry holohedry that the
    # deposited cell only approximates (e.g. nearly-orthorhombic monoclinic
    # admits a 180°-about-z that holds in true orthorhombic but is off by
    # ~few % in the actual monoclinic metric).  Strict 1e-6 tolerance
    # misses these.  Loosening to 5% surfaces them as discrete altindex
    # candidates — preserving experimental Fobs through the reindex+refine
    # pathway rather than falling back to Fcalc-only.  Same-cell pairs are
    # unaffected (their natural ops are strictly metric-preserving).
    metric_tol_rel = 0.05 if stretched else 1e-6

    cell = mov_st.cell
    sg_name = mov_st.spacegroup_hm
    sg = gemmi.find_spacegroup_by_name(sg_name)
    if sg is None:
        _log(f"resolve_altindex: unknown SG {sg_name!r}; skipping")
        return _write_log_and_return(_no_op_resolve_result(mov_pdb, mov_mtz))

    O, Finv = _ortho_pair(cell)
    G = O.T @ O

    A_frac = (Finv @ A_cart.T).T
    B_frac = (Finv @ B_cart.T).T

    # For each discrete R candidate, evaluate the best achievable RMSD using
    # the COM-optimal *continuous* translation (rather than ranking only
    # discrete origin-table entries).  The previous per-atom `diff -=
    # np.round(diff)` wrap was a stand-in for "find the best t"; that works
    # when COM-offset is small (same-cell, aligned origins), but fails for
    # cross-cell pairs after stretch where the COM offset is some continuous
    # value not on the origin-table grid.  COM-optimal continuous t is what
    # the downstream `_apply_cart_transform_to_pdb` already uses to *apply*
    # the transform — using it for ranking too keeps rank and apply
    # consistent, and recovers the cross-cell case.
    #
    # For same-cell, well-aligned pairs this collapses to t_opt ≈ 0 and the
    # is_zero_t guard below falls through to action='none', preserving the
    # old behaviour.  For SG-symop cases the COM-optimal continuous t lands
    # in the same physical equivalence class as the intrinsic discrete t
    # (mod lattice).
    best = (rmsd_base, np.eye(3), np.zeros(3))   # (rmsd, R_frac, t_frac)
    for R_frac, _t_intrinsic in _enum_alt_rot_origin_candidates(
            cell, sg, sg_name, metric_tol_rel=metric_tol_rel):
        if not _is_metric_preserving(R_frac, G, atol=metric_tol_rel):
            continue
        A_rot_frac = (R_frac @ A_frac.T).T
        t_opt_frac = B_frac.mean(axis=0) - A_rot_frac.mean(axis=0)
        diff = (A_rot_frac + t_opt_frac) - B_frac
        cart_diff = (O @ diff.T).T
        rmsd_disc = float(np.sqrt(np.mean(np.sum(cart_diff**2, axis=1))))
        if rmsd_disc < best[0]:
            best = (rmsd_disc, R_frac, t_opt_frac)

    rmsd_top, R_frac, t_frac = best
    R_cart = O @ R_frac @ Finv
    t_cart_opt = (B_cart.mean(axis=0) - (R_cart @ A_cart.T).T.mean(axis=0))
    drot_deg = (_rot_deviation_deg(R_lsq, R_cart)
                if rmsd_top < rmsd_base else 0.0)

    is_identity_R = np.allclose(R_frac, np.eye(3), atol=1e-6)
    # Check the UNWRAPPED t: a t = (1, 0, 0) (whole lattice vector) is an
    # F-space no-op but a real-space ~|a| Å shift that must be applied to
    # the PDB to put moving in the same image of the lattice as ref (this
    # arises post-cell-stretch when COM offset spans a lattice vector).
    is_zero_t = float(np.linalg.norm(t_frac)) < 0.01

    # Is R_frac one of the SG's own proper rotations?  If so, this isn't a
    # genuine altindex — it's just a chain re-assignment via an existing
    # symmetry operator (e.g. insulin H 3, where the two deposits chose
    # different H 3-equivalent ASUs).  Reindexing Fobs by an in-SG op is
    # a no-op (|F| is SG-symmetric) but rotating + translating atoms then
    # re-refining will diverge: refmac sees atoms that have been moved to
    # a non-canonical origin choice while the Fobs still expects them at
    # the canonical origin.  Treat as `action='none'` and let
    # bend_fit_progressive's own symop search handle the chain matching.
    R_frac_int_chk = np.round(R_frac).astype(int)
    R_in_sg = False
    if np.allclose(R_frac, R_frac_int_chk, atol=1e-4):
        for op in sg.operations():
            sg_R = (np.array(op.rot, dtype=int) // op.DEN)
            if np.array_equal(sg_R, R_frac_int_chk):
                R_in_sg = True
                break

    if rmsd_top >= improve_threshold * rmsd_base:
        _log(f"  best discrete op: rmsd={rmsd_top:.2f} A vs baseline "
             f"{rmsd_base:.2f} A — no improvement, skipping")
        if stretched:
            return _write_log_and_return(_finish_after_stretch_only(
                mov_pdb, mov_mtz, outdir, refine_cycles, fill_asu,
                fill_method, drot_deg, rmsd_base, verbose))
        return _write_log_and_return(_no_op_resolve_result(mov_pdb, mov_mtz,
                                                            drot=drot_deg,
                                                            rmsd=rmsd_base))

    if is_identity_R and is_zero_t:
        _log(f"  best discrete op is identity — no transform needed")
        if stretched:
            return _write_log_and_return(_finish_after_stretch_only(
                mov_pdb, mov_mtz, outdir, refine_cycles, fill_asu,
                fill_method, drot_deg, rmsd_top, verbose))
        return _write_log_and_return(_no_op_resolve_result(mov_pdb, mov_mtz,
                                                            drot=drot_deg,
                                                            rmsd=rmsd_top))

    mov_stem = os.path.splitext(os.path.basename(mov_pdb))[0]
    mtz_stem = os.path.splitext(os.path.basename(mov_mtz))[0]

    if R_in_sg and not is_identity_R:
        # R is one of the SG's own rotations + an origin shift — not a
        # genuine altindex but an alternative origin choice + chain
        # re-assignment.  Apply the *discrete* SG op (R + exact discrete
        # t_frac, not the COM-aligned continuous t_cart_opt) to BOTH PDB
        # and MTZ via the SF transformation theorem:
        #     F'(h) = F(R^T h) · exp(2πi h · t)
        # See VALIDATE_ALTINDEX.md (formula validated to 1.24° mean Δphi
        # end-to-end through refmac on insulin H 3).  For in-SG R the
        # F(R^T h) lookup is a no-op on amplitudes (SG-symmetric data)
        # but the phase shift is essential to keep map + atoms consistent.
        action  = 'sg_op_origin'
        out_pdb = os.path.join(outdir, f'{mov_stem}_sgop.pdb')
        out_mtz = os.path.join(outdir, f'{mtz_stem}_sgop.mtz')
        t_cart_discrete = O @ t_frac
        _apply_cart_transform_to_pdb(mov_pdb, R_cart, t_cart_discrete,
                                     out_pdb, ref_cell=cell)
        R_frac_int = np.round(R_frac).astype(int)
        _apply_op_to_mtz(mov_mtz, R_frac_int, t_frac, out_mtz)
        _log(f"  action=sg_op_origin  R is an in-SG symop  "
             f"t_frac={t_frac}  rmsd_after={rmsd_top:.2f} A")
        _log(f"  applied (R, t) to PDB (cartesian) and to MTZ "
             f"(F-space, no re-refinement)")
        return _write_log_and_return(dict(
            mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
            action=action, R_frac=R_frac, t_frac=t_frac,
            drot_deg=drot_deg, rmsd_after=rmsd_top))

    if is_identity_R:
        # Origin-only: translate PDB, phase-shift MTZ.  No re-refinement.
        action = 'origin_only'
        out_pdb = os.path.join(outdir, f'{mov_stem}_origin.pdb')
        out_mtz = os.path.join(outdir, f'{mtz_stem}_origin.mtz')
        _apply_cart_transform_to_pdb(mov_pdb, R_cart, t_cart_opt, out_pdb,
                                     ref_cell=cell)
        _shift_mtz_origin(mov_mtz, t_frac, out_mtz)
        _log(f"  action=origin_only  t_frac={t_frac}  "
             f"rmsd_after={rmsd_top:.2f} A")
        return _write_log_and_return(dict(
            mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
            action=action, R_frac=R_frac, t_frac=t_frac,
            drot_deg=drot_deg, rmsd_after=rmsd_top))

    # Altindex (R != I, with or without origin shift): transform PDB, reindex
    # MTZ Fobs, re-refine.
    action = 'altindex_refine'
    pre_pdb = os.path.join(outdir, f'{mov_stem}_alt.pdb')
    pre_mtz = os.path.join(outdir, f'{mtz_stem}_alt.mtz')
    _apply_cart_transform_to_pdb(mov_pdb, R_cart, t_cart_opt, pre_pdb,
                                 ref_cell=cell)
    R_frac_int = np.round(R_frac).astype(int)
    if not np.allclose(R_frac, R_frac_int, atol=1e-4):
        raise ValueError(f"R_frac not integer: {R_frac}")
    _reindex_mtz_fobs(mov_mtz, R_frac_int, pre_mtz)

    _log(f"  action=altindex_refine  drot={drot_deg:.2f} deg  "
         f"rmsd_after={rmsd_top:.2f} A")
    _log(f"  R_frac =\n{R_frac_int}\n  t_frac = {t_frac}")
    _log(f"  re-refining {os.path.basename(pre_pdb)} against "
         f"{os.path.basename(pre_mtz)} ...")

    # Re-refinement of the altindex-reindexed MTZ: reindex generally drops
    # ~half the reflections (the reindex is a SG-asymmetric reshuffle, and
    # any unmatched HKLs become missing).  Propagate fill_asu so the
    # completeness gate inside run_refinement engages here too.
    out_pdb, out_mtz = run_refinement(pre_pdb, pre_mtz, outdir=outdir,
                                       n_cycles=refine_cycles,
                                       fill_asu=fill_asu,
                                       fill_method=fill_method)

    return _write_log_and_return(dict(
        mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
        action=action, R_frac=R_frac, t_frac=t_frac,
        drot_deg=drot_deg, rmsd_after=rmsd_top))


def main():
    argv = sys.argv[1:]
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)

    pdb1, pdb2, mapfile, mtz1, mtz2, p = _parse_args(argv)

    if pdb1 is None or pdb2 is None:
        print("ERROR: provide two PDB files on the command line.", file=sys.stderr)
        sys.exit(1)

    # ── MTZ handling: column detection (refinement deferred to fitreso_scan
    # when scan_dir is set, or to the single-fit path below)
    mov_mtz = mtz1
    ref_mtz = mtz2

    # ── fitreso_scan mode (refinement + altindex resolution happen inside) ───
    if p['scan_dir']:
        # mapfile may be: None, a single str (one map), or a (mov, ref) tuple
        if isinstance(mapfile, tuple):
            mov_map, ref_map = mapfile
        else:
            mov_map = ref_map = mapfile
        mov_src = mov_mtz or mov_map
        ref_src = ref_mtz or ref_map
        if mov_src is None or ref_src is None:
            print("ERROR: scan_dir requires two map/MTZ inputs.", file=sys.stderr)
            sys.exit(1)
        fitreso_scan(pdb1, pdb2, mov_src, ref_src, p['scan_dir'],
                     f_col=p['f_col'], phi_col=p['phi_col'],
                     run_refinement_flag=p['run_refinement'],
                     refine_cycles=p['refine_cycles'],
                     fill_asu=p['fill_asu'],
                     fill_method=p['fill_method'],
                     sample_rate=p['sample_rate'],
                     subtract=p['subtract'],
                     scan_all_fr=p['scan_all_fr'],
                     early_stop_tol=p['early_stop_tol'],
                     early_stop_n=p['early_stop_n'])
        sys.exit(0)

    # ── Single-fit path: refinement (if requested) then altindex resolution ──
    if mov_mtz and p['run_refinement']:
        pdb1, mov_mtz = run_refinement(pdb1, mov_mtz, outdir='refine_mov',
                                        n_cycles=p['refine_cycles'],
                                        fill_asu=p['fill_asu'],
                                        fill_method=p['fill_method'])
    if ref_mtz and p['run_refinement']:
        pdb2, ref_mtz = run_refinement(pdb2, ref_mtz, outdir='refine_ref',
                                        n_cycles=p['refine_cycles'],
                                        fill_asu=p['fill_asu'],
                                        fill_method=p['fill_method'])

    if mov_mtz:
        _res = resolve_altindex(pdb1, pdb2, mov_mtz,
                                 outdir='altindex_resolve',
                                 refine_cycles=p['refine_cycles'],
                                 fill_asu=p['fill_asu'],
                                 fill_method=p['fill_method'])
        pdb1    = _res['mov_pdb_out']
        mov_mtz = _res['mov_mtz_out']

    col_suffix = ''
    if mov_mtz:
        f_col, phi_col = detect_2fofc_cols(mov_mtz, p['f_col'], p['phi_col'])
        print(f"moving MTZ: {os.path.basename(mov_mtz)}  columns {f_col}/{phi_col}",
              flush=True)
        _CANONICAL_F = {'FWT', '2FOFCWT'}
        col_suffix = '' if f_col in _CANONICAL_F else f'_{f_col}'

    # ── single fit ───────────────────────────────────────────────────────────
    if p['fitparams'] and os.path.exists(p['fitparams']):
        print(f"loading fitparams from {p['fitparams']}")
        hkls, AB, active, snr, cell1, cell2, dims, rmsd0 = load_fitparams(p['fitparams'])
        nhkls = int(active.sum())
    else:
        result = bend_fit(
            pdb1, pdb2,
            nhkls=p['nhkls'], fitreso=p['fitreso'],
            drop_snr=p['drop_snr'],
            dimensions=p['dimensions'], geotest=p['geotest'], frac=p['frac'],
            use_symm=p['use_symm'],
        )
        nhkls = int(result.active.sum())

        fp_path = save_fitparams(
            'psdvf', result.hkls, result.AB, result.active, result.snr,
            result.cell1, result.cell2, result.dimensions, result.rmsd)
        print(f"fitparams saved to {fp_path}")

        hkls, AB, active, snr = result.hkls, result.AB, result.active, result.snr
        cell1, cell2, dims    = result.cell1, result.cell2, result.dimensions

    # Build AB_xyz (3, N, 2)
    dim_list = list(dims)
    xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
    AB_xyz = np.zeros((3, len(hkls), 2))
    for new_i, old_i in enumerate(xyz_dims):
        AB_xyz[new_i] = AB[old_i]

    # Write bent PDB
    bent_pdb = f"bent{nhkls}.pdb"
    write_bent_pdb(pdb1, pdb2, hkls, AB_xyz, bent_pdb, frac=p['frac'])
    print(f"bent PDB: {bent_pdb}")

    # Geometry check
    if p['geotest']:
        print("evaluating chemical geometry after distortion")
        try:
            result_geo = subprocess.run(
                ['refmac5', 'xyzin', bent_pdb, 'xyzout', '/dev/null'],
                input='refi type ideal\nncyc 1\n',
                capture_output=True, text=True
            )
            for line in result_geo.stdout.splitlines():
                if 'rmsBOND' in line or 'rmsANGLE' in line:
                    print(line.strip())
        except FileNotFoundError:
            print("refmac5 not found — skipping geotest")

    # Bent map — from MTZ or CCP4 map
    if mov_mtz:
        mov_d, mov_h = mtz_to_map_data(mov_mtz, f_col, phi_col,
                                        sample_rate=p['sample_rate'])
        bent_map_out = f"bent{nhkls}{col_suffix}.map"
        _make_bent_map_data(mov_d, mov_h, hkls, AB_xyz, bent_map_out, cell2=cell2)
        if p['deltamaps']:
            _make_delta_maps_hdr(mov_h, hkls, AB_xyz, nhkls)
    elif mapfile:
        # Single-fit fallback: use first map only (mov)
        single_map = mapfile[0] if isinstance(mapfile, tuple) else mapfile
        if not os.path.exists(single_map):
            print(f"WARNING: map file not found: {single_map}", file=sys.stderr)
        else:
            bent_map_out = f"bent{nhkls}.map"
            make_bent_map(single_map, hkls, AB_xyz, bent_map_out, cell2=cell2)
            if p['deltamaps']:
                make_delta_maps(single_map, hkls, AB_xyz, nhkls)


if __name__ == '__main__':
    main()
