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

Public API:
    from bendfinder import bend_fit, bend_apply_pdb, bend_apply_map
    from bendfinder import fitreso_scan, detect_2fofc_cols, mtz_to_map_data
    result = bend_fit("moving.pdb", "reference.pdb", nhkls=30)
    print(result.rmsd)
"""

import sys
import os
import time
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import lstsq as sp_lstsq, svd as sp_svd
from scipy.ndimage import map_coordinates
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


def match_atoms(atoms1, atoms2):
    """Match atoms between two P1 lists by uid.

    Returns
    -------
    fitme   : ndarray (N, 8) — xf yf zf  dxf dyf dzf  do  dB
    ca_mask : ndarray (N,) bool
    uids    : list[str]
    bfacs   : ndarray (N, 2) — [B1, B2] for each matched pair
    """
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


def fit_lstsq_symm(X, b, drop_snr=1.0, max_rounds=5):
    """Joint 3D lstsq with per-canonical-HKL SNR pruning.

    X : (3N, 6M) symmetry-adapted design matrix
    b : (3N,)    interleaved shifts [Δx_0,Δy_0,Δz_0, Δx_1,...]

    Returns
    -------
    params       : (6M,) float — fitted (A,B) vectors for each canonical HKL
    active_canon : (M,)  bool  — surviving canonical HKLs
    snr_canon    : (M,)  float — SNR per canonical HKL
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

        params_full[col_mask] = sol
        snr_full[ai] = snr_a

        if drop_snr <= 0:
            break
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


def fit_lstsq(X, shifts, drop_snr=1.0, max_rounds=5):
    """Linear fit with iterative SNR pruning.

    SVD of the design matrix is computed ONCE per round and shared across all
    dimensions — ~3× faster than per-dimension SVD.

    Parameters
    ----------
    X        : (N_atoms, 2*N_hkls) design matrix
    shifts   : (N_atoms, N_dims) fractional shifts
    drop_snr : float — SNR threshold (0 = disable)

    Returns
    -------
    AB     : (N_dims, N_hkls, 2) — A (sin) and B (cos) coefficients
    active : (N_hkls,) bool — surviving HKLs
    snr    : (N_hkls,) float — mean SNR across dimensions
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

    mtz.add_dataset('bendfinder')
    mtz.add_column('H', 'H')
    mtz.add_column('K', 'H')
    mtz.add_column('L', 'H')
    for dim in dimensions:
        mtz.add_column(f'd{dim.upper()}', 'F')   # amplitude of d{x,y,z} shift (Å)
        mtz.add_column(f'PH{dim.upper()}', 'P')  # phase (degrees)
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
    n_cols = 3 + 2 * n_dims + 2
    data = np.zeros((n_hkls, n_cols), dtype=np.float32)
    data[:, 0] = hkls[:, 0]
    data[:, 1] = hkls[:, 1]
    data[:, 2] = hkls[:, 2]
    for d, dim in enumerate(dimensions):
        A, B = AB[d, :, 0], AB[d, :, 1]
        amp_frac = np.sqrt(A**2 + B**2)
        data[:, 3 + 2*d]     = amp_frac * _AXIS_LEN[dim]
        data[:, 3 + 2*d + 1] = np.degrees(np.arctan2(-A, B))
    data[:, 3 + 2*n_dims]     = snr
    data[:, 3 + 2*n_dims + 1] = active.astype(np.float32)

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
def _get_altindex_ops_cached(sg_name, cell_tuple, max_M, det_max):
    """LRU-cached worker for _get_altindex_ops.  Cell is a tuple for
    hashability; cell=None is encoded as None.  ~minute-long enumeration
    runs only once per (sg, cell, max_M, det_max) within a process."""
    return _get_altindex_ops_impl(sg_name, cell_tuple, max_M, det_max)


def _get_altindex_ops(sg_name, cell=None, max_M=1, det_max=1):
    if cell is None:
        cell_tuple = None
    elif hasattr(cell, 'a'):       # gemmi.UnitCell
        cell_tuple = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    else:
        cell_tuple = tuple(cell)
    return _get_altindex_ops_cached(sg_name, cell_tuple, max_M, det_max)


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


def _get_altindex_ops_impl(sg_name, cell=None, max_M=1, det_max=1):
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
    tol_G = 1e-6 * float(np.max(np.abs(G)))
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
    if ca_rmsd > 5.0:
        raise ValueError(f"CA RMSD after origin fix {ca_rmsd:.3f} Å > 5.0 Å. "
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
                         use_symm=True, dimensions='xyz',
                         frac=1.0, verbose=True,
                         altloc_filter=False, altloc_fallback=1.0,
                         iter_callback=None, max_canon=None,
                         precomputed_origin=None):
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

    # Permanent initial outlier rejection (B-factor + shift MAD)
    if verbose:
        print("rejecting outliers...", end=' ', flush=True)
    fitme, ca_mask, uids, bfacs, _ = reject_outliers(
        fitme, ca_mask, uids, cell1, mad_sigma=outlier_sigma,
        bfacs=bfacs, b_sigma=b_sigma, verbose=verbose)
    if verbose:
        print(f"{ca_mask.sum()} CA remaining", flush=True)

    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    ca_orth = frac_to_orth(ca_shifts, cell1)
    ca_rmsd = float(np.sqrt(np.mean(np.sum(ca_orth**2, axis=1))))
    if verbose:
        print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}", flush=True)
        print(f"CA RMSD after origin fix: {ca_rmsd:.3f} Å  ({ca_mask.sum()} CA pairs)", flush=True)
    if ca_rmsd > 5.0:
        raise ValueError(f"CA RMSD after origin fix {ca_rmsd:.3f} Å > 5.0 Å. "
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

    # Initial batch: all HKLs with d >= fitreso_start
    n_initial = int(np.sum(d_all >= fitreso_start - 1e-9))
    n_initial = max(n_initial, batch_hkls)         # at least one batch
    n_initial = min(n_initial, len(hkl_ndc))
    n_used = n_initial + 1                         # +1 for DC at index 0

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
            asu_fit_mask = asu_mask & fit_active
            frac_asu     = fitme[asu_fit_mask, :3]
            shifts_xyz   = fitme[asu_fit_mask][:, [dim_col[d] for d in 'xyz']]
            shifts_flat  = shifts_xyz.ravel()

            X_symm = build_design_matrix_symm(frac_asu, canon_hkls, proper_ops,
                                               hkl_to_canon, hkl_all_ops)
            params, active_c, snr_c = fit_lstsq_symm(
                X_symm, shifts_flat, drop_snr=drop_snr)
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
            frac_fit = fitme[fit_active, :3]
            shifts   = fitme[fit_active][:, dim_indices]
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
            print(f"  [iter {iter_count:2d}] fitreso={eff_reso:.2f}Å  "
                  f"nhkls={n_hkls-1}  canon={n_canon}  OD={od_ratio:.2f}  "
                  f"active={n_active}  RMSD={rmsd:.3f}Å{drop_str}  {dt:.0f}s", flush=True)

        result_AB     = AB
        result_active = active
        result_snr    = snr
        result_hkls   = hkls_now
        result_rmsd   = rmsd
        iter_count   += 1

        if iter_callback is not None:
            iter_callback(iter_count, n_hkls - 1, n_canon, rmsd,
                          hkls_now, AB_xyz_eval, active, snr)

        if max_canon is not None and n_canon >= max_canon:
            break

        # Advance to next batch
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
        'fill_fcalc':     False,
        'sample_rate':    0.0,
        'scan_dir':       None,
        'subtract':       'ref',
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
            elif key in ('fill_fcalc', 'fill-fcalc'):
                params['fill_fcalc'] = val.lower() not in ('false', '0', 'no')
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
        elif arg in ('deltamaps', 'delta', 'nofit', 'run_refinement', 'refine',
                     'fill_fcalc', 'fill-fcalc', '--fill-fcalc'):
            if arg in ('deltamaps', 'delta'):
                params['deltamaps'] = True
            elif arg in ('run_refinement', 'refine'):
                params['run_refinement'] = True
            elif arg in ('fill_fcalc', 'fill-fcalc', '--fill-fcalc'):
                params['fill_fcalc'] = True
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
            out[lbl] = (val, a, dist)
        else:
            out[lbl] = (0., _dummy, 0.)
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
    excluding systematic absences and (0,0,0).  Counts only rows with finite
    F (column auto-detected: FP / F / FOBS); rows with NaN F are treated as
    missing observations even when the HKL is present.
    """
    mtz = gemmi.read_mtz_file(mtz_path)
    data = np.array(mtz.array)
    lbls = mtz.column_labels()
    F_lbl = next((l for l in ('FP', 'F', 'FOBS') if l in lbls), None)
    if F_lbl is None:
        # No Fobs column (already a map-coef MTZ, e.g. post-refmac with just
        # FWT/PHWT).  Refmac will choke on this anyway; return "complete"
        # so the completeness gate doesn't preempt that error with a
        # confusing one.
        return 0, 0, 1.0
    hkls = data[:, :3].astype(int)
    Fcol = data[:, lbls.index(F_lbl)]
    proper_ops = get_proper_symops(mtz.spacegroup.xhm())
    R_ints = [np.round(R).astype(int) for R, _ in proper_ops] or [np.eye(3, dtype=int)]
    in_canon_with_F = {_canonical_hkl(tuple(int(x) for x in h), R_ints)
                       for h, f in zip(hkls, Fcol) if np.isfinite(f)}
    cell = mtz.cell
    d_min = float(mtz.resolution_high())
    asu = _enumerate_sg_asu_hkls(
        (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma),
        mtz.spacegroup, d_min)
    n_expected = len(asu)
    n_obs = sum(1 for h in asu
                if tuple(int(x) for x in h) in in_canon_with_F)
    frac = n_obs / n_expected if n_expected else 0.0
    return n_obs, n_expected, frac


def _fcalc_at_hkls(pdb_path, hkls, d_min, blur=0.0):
    """Fcalc at given HKLs from PDB model via gemmi DensityCalculatorX + numpy FFT.

    `blur` adds Babinet-style B-blur to the density (refmac uses ≈ 60 Å²
    by default to keep the grid compact; set to 0 for unblurred density —
    we deblur via `addends` after FFT).
    Returns (M,) complex128 in gemmi's F-convention.
    """
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    dc = gemmi.DensityCalculatorX()
    dc.d_min = float(d_min)
    if blur > 0:
        dc.blur = float(blur)
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    grid_xyz = np.array(dc.grid.array, dtype=np.float64)
    F = _fft_lookup_at_hkls(grid_xyz, np.asarray(hkls, dtype=np.int32),
                            st.cell.volume)
    if blur > 0:
        # Undo B-blur: F_true = F_blurred * exp(+blur/4 * s²)
        cell = st.cell
        M_frac = np.array(cell.frac.mat.tolist())
        h_orth = np.asarray(hkls, dtype=float) @ M_frac
        s_sq   = np.sum(h_orth ** 2, axis=1)
        F = F * np.exp(blur / 4.0 * s_sq)
    return F


def _fill_missing_with_fcalc(in_mtz_path, pdb_path, out_mtz_path, d_min=None):
    """Augment an MTZ with model-derived Fcalc at HKLs missing from its SG ASU.

    Used by `run_refinement` to give refmac a complete-ASU input so its
    FWT/PHWT output has no missing-data gaps.  Lighter than uniqueify in
    that it preserves the input column layout exactly and assigns FREE=0
    (work set) to filled rows without touching existing free-R flags.

    For filled rows the F column gets Fcalc and SIGFP gets the median of
    observed SIGFP (so refmac sees a finite weight).  If the input has no
    F/SIGFP column, returns the input path unchanged.
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

    # Identify amplitude / sigma / free columns
    F_lbl    = next((l for l in ('FP', 'F', 'FOBS') if l in lbls), None)
    SIG_lbl  = next((l for l in ('SIGFP', 'SIGF', 'SIGFOBS') if l in lbls), None)
    FREE_lbl = next((l for l in ('FREE', 'RFREE', 'FreeR_flag') if l in lbls), None)
    if F_lbl is None:
        # Nothing to fill — input MTZ doesn't have an Fobs column refmac would use.
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
    fc = _fcalc_at_hkls(pdb_path, missing_hkl, d_min)
    fc_amp = np.abs(fc).astype(np.float32)

    # Fill values for ancillary columns
    F_idx = lbls.index(F_lbl)
    sig_default = (np.float32(np.nanmedian(data[:, lbls.index(SIG_lbl)]))
                   if SIG_lbl else np.float32(1.0))
    free_default = np.float32(0.0)   # work set

    new_rows = np.full((len(missing_hkl), data.shape[1]),
                       np.nan, dtype=np.float32)
    new_rows[:, :3] = missing_hkl
    new_rows[:, F_idx] = fc_amp
    if SIG_lbl:
        new_rows[:, lbls.index(SIG_lbl)] = sig_default
    if FREE_lbl:
        new_rows[:, lbls.index(FREE_lbl)] = free_default

    out_data = np.vstack([data, new_rows]).astype(np.float32)
    mtz.set_data(out_data)
    mtz.write_to_file(out_mtz_path)
    return out_mtz_path


def _fft_box_d_min(cell_obj, NX, NY, NZ):
    """Resolution limit (Å) supported by an FFT grid (NX,NY,NZ) on `cell_obj`.

    Conservative: takes the per-axis Nyquist min, so HKLs enumerated to this
    d_min always fit inside the FFT box.
    """
    return float(min(2.0 * cell_obj.a / NX,
                     2.0 * cell_obj.b / NY,
                     2.0 * cell_obj.c / NZ))


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
    fill_fcalc=False,
    sample_rate=0.0,
    mov_fullcell=None,
    fitreso_list=(20, 15, 12, 10, 8, 7, 6, 5),
    max_hkl_scan=10,
    outlier_sigma=2.5, b_sigma=3.0, drop_snr=0.0, od_margin=1.5,
    batch_hkls=100, chunk_size=50000,
    riso_n_cycles=4, riso_sigma_cut=float('inf'),
    subtract='ref',
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

    # ── optional refinement ───────────────────────────────────────────────────
    # Refmac can drop disordered atoms during the run, so subsequent steps
    # (altindex resolution + re-refinement, write_bent_pdb, peak hunting)
    # all use the *refined* PDB to keep the atom set consistent.
    if run_refinement_flag and ref_mtz.lower().endswith('.mtz'):
        ref_pdb, ref_mtz = run_refinement(
            ref_pdb, ref_mtz,
            outdir=os.path.join(scan_dir, 'refine_ref'),
            n_cycles=refine_cycles, fill_fcalc=fill_fcalc)

    # ── altindex / origin resolution ──────────────────────────────────────────
    # Discrete (rotation × symop × origin) enumeration ranked by post-transform
    # CA RMSD.  On a non-trivial hit, mov PDB+MTZ are rewritten and re-refined
    # (or phase-shifted, for origin-only) so the rest of the scan sees an
    # already-aligned mov.
    #
    # Run on raw mov (NOT pre-refmac'd): refmac on raw mov against its own
    # data shifts the model by ~0.5 Å, which can move it out of the basin of
    # attraction of the best discrete altindex op (porin 3poq→3pou: best op
    # gives 2.11 Å on raw mov, 38.8 Å on refined mov — refmac shift crosses
    # a fractional boundary and the wrap-to-nearest-image hops to a wrong
    # equivalent).  resolve_altindex's altindex_refine action re-refines
    # internally on the reindexed MTZ, so mov gets a refinement pass in the
    # correct frame anyway.
    if mov_mtz.lower().endswith('.mtz'):
        _res = resolve_altindex(mov_pdb, ref_pdb, mov_mtz,
                                 outdir=os.path.join(scan_dir, 'altindex_resolve'),
                                 refine_cycles=refine_cycles, verbose=verbose,
                                 fill_fcalc=fill_fcalc)
        mov_pdb = _res['mov_pdb_out']
        mov_mtz = _res['mov_mtz_out']

    # When altindex_resolve returned action='none' (already aligned, no
    # altindex needed), we still need refmac on mov to generate FWT/PHWT.
    if (run_refinement_flag and mov_mtz.lower().endswith('.mtz')
            and 'FWT' not in {c.label for c in gemmi.read_mtz_file(mov_mtz).columns}):
        mov_pdb, mov_mtz = run_refinement(
            mov_pdb, mov_mtz,
            outdir=os.path.join(scan_dir, 'refine_mov'),
            n_cycles=refine_cycles, fill_fcalc=fill_fcalc)

    # ── load maps ─────────────────────────────────────────────────────────────
    # For CCP4 ASU map inputs, expand mov to the full unit cell so interpolation
    # has correct periodic boundary conditions (pad_mode='wrap').  Without this,
    # the bent map's boundary samples (y≈0.5 etc. for half-cell ASU maps) reflect
    # interior density unrelated to ref — producing a noise plane at the ASU
    # boundary in diff_norm.map.  This was historically controlled by the
    # `--fullcell-mov` flag in the standalone fitreso_scan.py; now always on
    # for CCP4 mov input.  MTZ inputs come out full-cell by construction.
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

    # Default mov_fullcell to True for CCP4 input (avoids ASU-boundary noise);
    # MTZ input ignores the flag (already full-cell).
    if mov_fullcell is None:
        mov_fullcell = mov_mtz.lower().endswith(('.map', '.ccp4', '.mrc'))
    mov_d, mov_h, mov_mtz_resolved, mov_f_col = _load(mov_mtz, 'mov',
                                                        fullcell=mov_fullcell)
    ref_d, ref_h, ref_mtz_resolved, ref_f_col = _load(ref_mtz, 'ref')

    col_suffix = ('' if (mov_f_col is None or mov_f_col in _CANONICAL_F)
                  else f'_{mov_f_col}')
    pad_mode   = 'wrap' if (mov_mtz_resolved is not None or mov_fullcell) else 'reflect'

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

    M     = _scan_cell_matrix(ref_h)
    # Atoms used for peak-labelling: union of ref and mov atoms.  Both
    # fractionalized in ref's cell so distances are computed in a single
    # frame (mov has already been moved into ref's frame by resolve_altindex
    # when applicable).  Each atom carries an `src` tag ('r' or 'm') shown
    # in the peak label.
    _ref_cell = gemmi.UnitCell(*ref_h['cell'])
    atoms     = (_scan_read_pdb_atoms_p1(ref_pdb, src='ref',
                                          cell_override=_ref_cell)
                 + _scan_read_pdb_atoms_p1(mov_pdb, src='mov',
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
    with open(os.devnull, 'w') as _dev, _contextlib.redirect_stdout(_dev):
        _a1, _cell1, _ = expand_to_p1(mov_pdb)
        _a2, _cell2, _ = expand_to_p1(ref_pdb)
        _fm, _ca, _uid, _bf = match_atoms(_a1, _a2)
        _fm, _ca, _, _, _ = reject_outliers(_fm, _ca, _uid, _cell1,
                                            mad_sigma=outlier_sigma,
                                            b_sigma=b_sigma, bfacs=_bf)
    raw_rmsd = rmsd_ca(_fm, _ca, _h0, _AB0, _cell2)

    # ── post-resolve hkl00 RMSD ──────────────────────────────────────────────
    # mov has already been aligned by resolve_altindex (above); the hkl00 RMSD
    # is just the matched-atom RMSD with no Fourier shift field applied.
    precomputed_origin = None
    with open(os.devnull, 'w') as _dev, _contextlib.redirect_stdout(_dev):
        _a1o, _cell1o, _sg1o = expand_to_p1(mov_pdb)
        _a2o, _cell2o, _sg2o = expand_to_p1(ref_pdb)
        _fmo, _cao, _uido, _bfo = match_atoms(_a1o, _a2o)
        _fmo, _cao, _, _, _ = reject_outliers(_fmo, _cao, _uido, _cell1o,
                                               mad_sigma=outlier_sigma,
                                               b_sigma=b_sigma, bfacs=_bfo)
    hkl00_rmsd = rmsd_ca(_fmo, _cao, _h0, _AB0, _cell2o)

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
              f"{'Rbent':>6}  {'Rbend':>6}  "
              f"{'peak':>7}  {'atom':>20}", flush=True)
        print('-' * 80, flush=True)
        print(f"{'pre':>7}  {raw_rmsd:>6.3f}  {0:>6d}  "
              f"{'   ---':>6}  {'   ---':>6}  "
              f"{'      ':>7}  (pre-origin alignment)", flush=True)

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
        if hkls is not None and AB_xyz is not None:
            write_bent_pdb(mov_pdb, ref_pdb, hkls, AB_xyz, f'{outdir}/bent.pdb')
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
        row_str = (f"{label:>7}  {after_str:>6}  {n_active:>6d}  "
                   f"{rbent_str:>6}  {rbend_str:>6}  "
                   f"{p_top[0]:>+7.2f}σ  {atom_str:>20} {p_top[2]:.2f}Å  "
                   f"[{t_elapsed:.0f}s]")
        _log_rows.append(dict(label=label, rmsd=rmsd_after, n_active=n_active,
                              rbent=riso, rbend=rbend, k=kF, B=B,
                              peak_sigma=p_top[0], peak_atom=atom_str,
                              peak_dist=p_top[2], t=t_elapsed,
                              row_str=row_str))
        if verbose:
            print(row_str, flush=True)

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
    _relsymlink_dir(os.path.abspath(mov_pdb),
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
        drop_snr=drop_snr, outlier_sigma=outlier_sigma, b_sigma=b_sigma,
        verbose=False,
        iter_callback=_hkl_callback,
        precomputed_origin=precomputed_origin,
    )

    # ── Section 3: fr<N> — resolution scans ──────────────────────────────────
    if verbose:
        print('-' * 122, flush=True)
    for fitreso in fitreso_list:
        outdir = os.path.join(scan_dir, f'fr{fitreso}')
        os.makedirs(outdir, exist_ok=True)
        t0 = time.time()
        result = bend_fit_progressive(
            mov_pdb, ref_pdb,
            fitreso_start=20.0, fitreso_end=float(fitreso),
            drop_snr=drop_snr, batch_hkls=batch_hkls, od_margin=od_margin,
            outlier_sigma=outlier_sigma, b_sigma=b_sigma,
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
        _save_point(f'fr{fitreso}', outdir, bent_map, hkl00_rmsd, result.rmsd,
                    int(result.active.sum()), t_fit,
                    hkls=result.hkls, AB_xyz=AB_xyz)
        save_fitparams(f'{outdir}/PSDVF.mtz',
                       result.hkls, result.AB, result.active, result.snr,
                       result.cell1, result.cell2, result.dimensions, result.rmsd)

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
        fh.write(f"# k, B    : F-space scale + isotropic B (Å²) used for Rbent\n\n")
        hdr = (f"{'label':>7}  {'RMSD':>6}  {'active':>6}  "
               f"{'Rbent':>6}  {'Rbend':>6}  "
               f"{'k':>6}  {'B(Å²)':>7}  "
               f"{'peak':>7}  {'atom':>20}  {'dist':>6}  {'t(s)':>5}\n")
        fh.write(hdr)
        fh.write('-' * (len(hdr) - 1) + '\n')
        for r in _log_rows:
            rmsd_s = f"{r['rmsd']:.3f}" if r['rmsd']  is not None else '  ---'
            rbe_s  = f"{r['rbent']*100:.1f}%" if r['rbent'] is not None else '  N/A'
            rbd_s  = f"{r['rbend']*100:.1f}%" if r['rbend'] is not None else '  N/A'
            k_s    = f"{r['k']:.3f}"   if r['k'] is not None else '  N/A'
            B_s    = f"{r['B']:+6.2f}" if r['B'] is not None else '   N/A'
            fh.write(f"{r['label']:>7}  {rmsd_s:>6}  {r['n_active']:>6d}  "
                     f"{rbe_s:>6}  {rbd_s:>6}  "
                     f"{k_s:>6}  {B_s:>7}  "
                     f"{r['peak_sigma']:>+7.2f}σ  {r['peak_atom']:>20}  "
                     f"{r['peak_dist']:>5.2f}Å  {r['t']:>4.0f}s\n")
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
                   fill_fcalc=False, completeness_threshold=0.99):
    """Run refmac5 or phenix.refine on pdb_path+mtz_path; return
    (refined_pdb_path, refined_mtz_path).

    Tries refmac5 first ($CCP4/bin or PATH), then phenix.refine.
    Calls sys.exit(1) if both fail.

    Returning both paths lets the caller chain subsequent operations
    (altindex resolution, peak hunting, etc.) on the refined atom set —
    which can differ from the input atom set when refmac drops
    disordered residues during the run.

    fill_fcalc : if True, augment input MTZ with model-derived Fcalc at
                 SG-ASU HKLs missing from the deposited dataset before
                 calling refmac.  Required when completeness is below
                 `completeness_threshold` (default 99%) — refmac only
                 writes FWT for HKLs present in input, so an incomplete
                 input means incomplete FWT and incomplete downstream
                 bent.mtz coverage.  When False and input is below
                 threshold, sys.exit(1) with instructions.
    completeness_threshold : minimum SG-ASU completeness required when
                 fill_fcalc=False.
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
    # unless the caller opts into `fill_fcalc=True`, which augments the
    # input MTZ with Fcalc-derived rows at SG-ASU HKLs missing from the
    # deposited dataset.  This is a deliberate choice — silently filling
    # would risk biasing the refinement / downstream fit, so make it
    # explicit.
    n_obs, n_exp, frac = _mtz_completeness(mtz_path)
    print(f'  input completeness: {n_obs}/{n_exp} = {100*frac:.1f}%',
          flush=True)
    mtz_path_for_refmac = mtz_path
    if frac < completeness_threshold:
        if not fill_fcalc:
            print(f'ERROR: {os.path.basename(mtz_path)} is {100*frac:.1f}% '
                  f'complete (< {100*completeness_threshold:.0f}% threshold).\n'
                  f'  Re-run with fill_fcalc=True (or pass --fill-fcalc on CLI) '
                  f'to fill missing SG-ASU HKLs with model-derived Fcalc.\n'
                  f'  This is required so refmac writes FWT for every SG-ASU '
                  f'HKL and the downstream scan map is complete.',
                  file=sys.stderr)
            sys.exit(1)
        mtz_filled = os.path.join(outdir, f'{stem}_filled.mtz')
        try:
            _fill_missing_with_fcalc(mtz_path, pdb_path, mtz_filled)
            n2, e2, f2 = _mtz_completeness(mtz_filled)
            print(f'  filled with Fcalc → {n2}/{e2} = {100*f2:.1f}% complete',
                  flush=True)
            mtz_path_for_refmac = mtz_filled
        except Exception as e:
            print(f'  _fill_missing_with_fcalc failed ({e!r}); aborting',
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
        free   = next((l for l in ('FREE', 'RFREE', 'FreeR_flag')   if l in cols), None)
        if fp and sigfp:
            labin  = f'FP={fp} SIGFP={sigfp}' + (f' FREE={free}' if free else '')
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
        cmd = [phenix, pdb_path, mtz_path_for_refmac,
               f'refinement.main.cycles={n_cycles}',
               f'output.prefix={os.path.join(outdir, stem)}']
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
        candidate_mtz = os.path.join(outdir, f'{stem}_refine_001.mtz')
        candidate_pdb = os.path.join(outdir, f'{stem}_refine_001.pdb')
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
    """Apply F'(H) = exp(-2πi H·t) F(H) to every (F, PHI) column pair in mtz_in.
    All other columns (intensities, sigF, FREE flag) pass through unchanged.
    Pure phase rotation — preserves all column relationships, no re-refinement."""
    mtz = gemmi.read_mtz_file(mtz_in)
    miller = mtz.make_miller_array().astype(np.float64)        # (N, 3)
    delta_phi_deg = (-360.0 * (miller @ np.asarray(t_frac, dtype=np.float64))
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
        data[:, target] = ((data[:, target] + delta_phi_deg) % 360.0
                            ).astype(np.float32)
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


def _enum_alt_rot_origin_candidates(cell, sg, sg_name):
    """Yield (R_frac, t_frac) pairs covering (altindex × symop × origin) candidates
    for resolve_altindex.  R_frac is fractional, t_frac is fractional mod 1.

    sg_name is the H-M symbol from the PDB header (e.g. 'H 3'); used to look up
    origin choices in _ORIGINS_TABLE.  gemmi's xhm() returns 'R 3:H' for H3 etc.,
    which doesn't match _ORIGINS_TABLE keys.

    Altindex rotations come from `_get_altindex_ops(sg_name, cell)` —
    metric-tensor enumeration of integer proper rotations preserving the
    cell's lattice that are NOT already in the SG.  This finds the full
    set of valid altindex candidates (including 6-folds along c for
    Sohncke trig/hex groups, axis permutations for orthorhombic, etc.),
    not just a hand-coded subset."""
    alt_ops = [(np.eye(3), np.zeros(3))] + list(_get_altindex_ops(sg_name,
                                                                    cell=cell))
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


def resolve_altindex(mov_pdb, ref_pdb, mov_mtz, outdir,
                     refine_cycles=5, improve_threshold=0.7,
                     verbose=True, fill_fcalc=False):
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
    """
    os.makedirs(outdir, exist_ok=True)
    mov_st = gemmi.read_pdb(mov_pdb)
    ref_st = gemmi.read_pdb(ref_pdb)
    mov_st.setup_entities()
    ref_st.setup_entities()

    A_cart, B_cart = _collect_matched_ca(mov_st, ref_st)
    if len(A_cart) < 3:
        if verbose:
            print(f"resolve_altindex: <3 CA matches; skipping", flush=True)
        return _no_op_resolve_result(mov_pdb, mov_mtz)

    rmsd_base = float(np.sqrt(np.mean(np.sum((A_cart - B_cart) ** 2, axis=1))))
    R_lsq, _t_lsq, rmsd_lsq = _kabsch(A_cart, B_cart)

    if verbose:
        print(f"resolve_altindex: {len(A_cart)} CA pairs  "
              f"baseline RMSD {rmsd_base:.2f} A  LSQ {rmsd_lsq:.2f} A",
              flush=True)

    cell = mov_st.cell
    sg_name = mov_st.spacegroup_hm
    sg = gemmi.find_spacegroup_by_name(sg_name)
    if sg is None:
        if verbose:
            print(f"resolve_altindex: unknown SG {sg_name!r}; skipping", flush=True)
        return _no_op_resolve_result(mov_pdb, mov_mtz)

    O, Finv = _ortho_pair(cell)
    G = O.T @ O

    A_frac = (Finv @ A_cart.T).T
    B_frac = (Finv @ B_cart.T).T

    best = (rmsd_base, np.eye(3), np.zeros(3))   # (rmsd, R_frac, t_frac)
    for R_frac, t_frac in _enum_alt_rot_origin_candidates(cell, sg, sg_name):
        if not _is_metric_preserving(R_frac, G):
            continue
        Af2  = (R_frac @ A_frac.T).T + t_frac
        diff = Af2 - B_frac
        diff -= np.round(diff)
        cart_diff = (O @ diff.T).T
        rmsd_disc = float(np.sqrt(np.mean(np.sum(cart_diff**2, axis=1))))
        if rmsd_disc < best[0]:
            best = (rmsd_disc, R_frac, t_frac)

    rmsd_top, R_frac, t_frac = best
    R_cart = O @ R_frac @ Finv
    t_cart_opt = (B_cart.mean(axis=0) - (R_cart @ A_cart.T).T.mean(axis=0))
    drot_deg = (_rot_deviation_deg(R_lsq, R_cart)
                if rmsd_top < rmsd_base else 0.0)

    is_identity_R = np.allclose(R_frac, np.eye(3), atol=1e-6)
    t_wrapped = t_frac - np.round(t_frac)
    is_zero_t = float(np.linalg.norm(t_wrapped)) < 0.01

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
        if verbose:
            print(f"  best discrete op: rmsd={rmsd_top:.2f} A vs baseline "
                  f"{rmsd_base:.2f} A — no improvement, skipping", flush=True)
        return _no_op_resolve_result(mov_pdb, mov_mtz,
                                     drot=drot_deg, rmsd=rmsd_base)

    if is_identity_R and is_zero_t:
        if verbose:
            print(f"  best discrete op is identity — no transform needed",
                  flush=True)
        return _no_op_resolve_result(mov_pdb, mov_mtz,
                                     drot=drot_deg, rmsd=rmsd_top)

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
        if verbose:
            print(f"  action=sg_op_origin  R is an in-SG symop  "
                  f"t_frac={t_frac}  rmsd_after={rmsd_top:.2f} A", flush=True)
            print(f"  applied (R, t) to PDB (cartesian) and to MTZ "
                  f"(F-space, no re-refinement)", flush=True)
        return dict(mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
                    action=action, R_frac=R_frac, t_frac=t_frac,
                    drot_deg=drot_deg, rmsd_after=rmsd_top)

    if is_identity_R:
        # Origin-only: translate PDB, phase-shift MTZ.  No re-refinement.
        action = 'origin_only'
        out_pdb = os.path.join(outdir, f'{mov_stem}_origin.pdb')
        out_mtz = os.path.join(outdir, f'{mtz_stem}_origin.mtz')
        _apply_cart_transform_to_pdb(mov_pdb, R_cart, t_cart_opt, out_pdb,
                                     ref_cell=cell)
        _shift_mtz_origin(mov_mtz, t_frac, out_mtz)
        if verbose:
            print(f"  action=origin_only  t_frac={t_frac}  "
                  f"rmsd_after={rmsd_top:.2f} A", flush=True)
        return dict(mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
                    action=action, R_frac=R_frac, t_frac=t_frac,
                    drot_deg=drot_deg, rmsd_after=rmsd_top)

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

    if verbose:
        print(f"  action=altindex_refine  drot={drot_deg:.2f} deg  "
              f"rmsd_after={rmsd_top:.2f} A", flush=True)
        print(f"  R_frac =\n{R_frac_int}\n  t_frac = {t_frac}", flush=True)
        print(f"  re-refining {os.path.basename(pre_pdb)} against "
              f"{os.path.basename(pre_mtz)} ...", flush=True)

    # Re-refinement of the altindex-reindexed MTZ: reindex generally drops
    # ~half the reflections (the reindex is a SG-asymmetric reshuffle, and
    # any unmatched HKLs become missing).  Propagate fill_fcalc so the
    # completeness gate inside run_refinement engages here too.
    out_pdb, out_mtz = run_refinement(pre_pdb, pre_mtz, outdir=outdir,
                                       n_cycles=refine_cycles,
                                       fill_fcalc=fill_fcalc)

    return dict(mov_pdb_out=out_pdb, mov_mtz_out=out_mtz,
                action=action, R_frac=R_frac, t_frac=t_frac,
                drot_deg=drot_deg, rmsd_after=rmsd_top)


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
                     fill_fcalc=p['fill_fcalc'],
                     sample_rate=p['sample_rate'],
                     subtract=p['subtract'])
        sys.exit(0)

    # ── Single-fit path: refinement (if requested) then altindex resolution ──
    if mov_mtz and p['run_refinement']:
        pdb1, mov_mtz = run_refinement(pdb1, mov_mtz, outdir='refine_mov',
                                        n_cycles=p['refine_cycles'],
                                        fill_fcalc=p['fill_fcalc'])
    if ref_mtz and p['run_refinement']:
        pdb2, ref_mtz = run_refinement(pdb2, ref_mtz, outdir='refine_ref',
                                        n_cycles=p['refine_cycles'],
                                        fill_fcalc=p['fill_fcalc'])

    if mov_mtz:
        _res = resolve_altindex(pdb1, pdb2, mov_mtz,
                                 outdir='altindex_resolve',
                                 refine_cycles=p['refine_cycles'])
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
