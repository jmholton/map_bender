#!/programs/ccp4-8.0/bin/ccp4-python
"""bendfinder.py — fit a Fourier shift field between two crystal forms.

Usage:
    bendfinder.py moving.pdb reference.pdb [map.ccp4] [key=value ...]

Parameters:
    nhkls=30        max Fourier terms (default 30); reports effective fitreso
    fitreso=X       use all HKLs with d ≥ X Å (e.g. fitreso=12); overrides nhkls
    drop_snr=1      post-fit SNR threshold (default 1; 0 = disable)
    frac=1          scale factor for applied shifts (default 1)
    geotest=false   run refmac5 geometry check (default false)
    use_symm=true   enforce space-group symmetry constraint (default true)
    dimensions=xyz  which shifts to fit (subset of xyzoBà; default xyz)
    fitparams=file  load pre-computed psdvf.mtz (or other fitparams) and skip fitting
    deltamaps       write delta-x/y/z/r maps in addition to bent map

Public API:
    from bendfinder import bend_fit, bend_apply_pdb, bend_apply_map
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
from scipy.linalg import lstsq as sp_lstsq
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
                    bfacs=None, b_sigma=None):
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
    if n_drop_shift:
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
        if n_drop_b:
            print(f"  B-factor outlier rejection: dropped additional {n_drop_b} atoms "
                  f"(B>{b_tol:.1f} Å²; median={b_med:.1f}, σ_MAD={b_sig_eq:.1f} Å²)")
        keep = keep & b_keep

    return (fitme[keep], ca_mask[keep],
            [u for u, k in zip(uids, keep) if k],
            bfacs[keep] if bfacs is not None else None,
            keep)


def expand_to_p1(pdb_path, altloc_filter=False, altloc_fallback=1.0):
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
    sg_name = sg_obj.hm

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

    X[3n+α, 6m+β]   = Σ_k R_k[β,α] sin(2π H_m·(R_k x_n + t_k))
    X[3n+α, 6m+3+β] = Σ_k R_k[β,α] cos(2π H_m·(R_k x_n + t_k))

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
        x_t = frac_coords @ R.T + t                    # (N, 3)
        ph  = TWO_PI * (x_t @ canon_f.T)               # (N, M)
        s = np.sin(ph)                                  # (N, M)
        c = np.cos(ph)
        # Zero columns for canonicals where this op maps outside the basis
        s[:, ~valid_mk[:, ki]] = 0.0
        c[:, ~valid_mk[:, ki]] = 0.0
        for alpha in range(3):
            for beta in range(3):
                r = R[beta, alpha]                      # (R^T)[alpha,beta]
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

        U, s, Vt = np.linalg.svd(Xa, full_matrices=False)
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
            phi  = TWO_PI * float(np.dot(Hc, t))
            cph, sph = np.cos(phi), np.sin(phi)
            dA = R.T @ (Ac * cph - Bc * sph)
            dB = R.T @ (Ac * sph + Bc * cph)
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
        U, s, Vt = np.linalg.svd(Xa, full_matrices=False)
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

def write_bent_pdb(pdb1_path, pdb2_path, hkls, AB_xyz, outpath, frac=1.0):
    """Apply shift field to pdb1, convert to pdb2 orthogonal frame, write PDB.

    AB_xyz : (3, N_hkls, 2) — x, y, z shift coefficients.
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
    d_frac    = eval_shift_field(frac_arr, hkls, AB_xyz)     # (N, 3)
    bent_frac = frac_arr + frac * d_frac                     # (N, 3)
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


def write_ccp4(path, data, hdr_template, cell_override=None):
    """Write a CCP4 map, borrowing the header from hdr_template.

    Overwrites DMIN/DMAX/DMEAN (words 20-22) and RMS (word 54) with
    actual stats.  ISPG (word 23, byte 88) is preserved from the template.
    cell_override: optional (a,b,c,al,be,ga) to replace cell in header.
    """
    raw = bytearray(hdr_template['raw_header'])

    mn, mx = float(data.min()), float(data.max())
    mean   = float(data.mean())
    rms    = float(np.sqrt(np.mean(data**2)))
    for off, val in zip((76, 80, 84, 216), (mn, mx, mean, rms)):
        struct.pack_into('<f', raw, off, val)

    if cell_override is not None:
        for i, val in enumerate(cell_override):
            struct.pack_into('<f', raw, 40 + i*4, val)

    with open(path, 'wb') as f:
        f.write(raw)
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


def interpolate_map(data, hdr, probe_frac):
    """Trilinear/tricubic interpolation of map data at fractional positions.

    probe_frac : (N, 3) fractional coords to sample.
    Returns (N,) float32 density values.
    """
    idx = _frac_to_grid_indices(probe_frac, hdr)   # (3, N)
    return map_coordinates(data, idx, order=3, mode='wrap').astype(np.float32)


# ══════════════════════════════════════════════════════════════════════════════
# Bent map (replace mapman + floatgen pipeline)
# ══════════════════════════════════════════════════════════════════════════════

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
    nc, nr, ns = hdr['nc'], hdr['nr'], hdr['ns']
    nx, ny, nz = hdr['nx'], hdr['ny'], hdr['nz']
    ncstart, nrstart, nsstart = hdr['ncstart'], hdr['nrstart'], hdr['nsstart']
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']

    # Build fractional coordinates of every grid point
    sec_idx = np.arange(ns)
    row_idx = np.arange(nr)
    col_idx = np.arange(nc)
    sec_g, row_g, col_g = np.meshgrid(sec_idx, row_idx, col_idx, indexing='ij')

    # Map column/row/section axes to x/y/z based on MAPC/MAPR/MAPS.
    # axis_to_frac[ax] = grid-index array for crystallographic axis ax (1=X,2=Y,3=Z)
    # start[ax]        = ASU origin along ax (belongs to MAPC/MAPR/MAPS, not X/Y/Z)
    axis_to_frac = {mapc: col_g, mapr: row_g, maps: sec_g}
    start  = {mapc: ncstart, mapr: nrstart, maps: nsstart}
    Nxyz   = {1: nx, 2: ny, 3: nz}   # full unit-cell grid sizes per axis

    frac_x = (axis_to_frac[1].ravel() + start[1]) / Nxyz[1]
    frac_y = (axis_to_frac[2].ravel() + start[2]) / Nxyz[2]
    frac_z = (axis_to_frac[3].ravel() + start[3]) / Nxyz[3]
    frac_pts = np.stack([frac_x, frac_y, frac_z], axis=1)   # (N, 3)

    # Shift field at each grid point (fractional)
    delta = eval_shift_field(frac_pts, hkls, AB_xyz)   # (N, 3)

    # Sample the map at (grid_point - shift): negative because we invert
    source_frac = frac_pts - delta                      # (N, 3)
    new_vals = interpolate_map(data, hdr, source_frac)  # (N,)

    if use_jacobian:
        # Scale each voxel by (1 - div(δ)) to conserve density under the
        # coordinate transformation (first-order Jacobian correction).
        div = eval_divergence(frac_pts, hkls, AB_xyz)  # (N,)
        new_vals = new_vals * (1.0 - div)

    new_data = new_vals.reshape(ns, nr, nc)
    write_ccp4(outpath, new_data, hdr, cell_override=cell2)
    print(f"bent map written to {outpath}")
    return outpath


def make_delta_maps(map_path, hkls, AB_xyz, nhkls):
    """Write delta-x, delta-y, delta-z, delta-r maps."""
    data, hdr = read_ccp4(map_path)
    nc, nr, ns = hdr['nc'], hdr['nr'], hdr['ns']

    sec_idx = np.arange(ns)
    row_idx = np.arange(nr)
    col_idx = np.arange(nc)
    sec_g, row_g, col_g = np.meshgrid(sec_idx, row_idx, col_idx, indexing='ij')
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']
    start  = {mapc: hdr['ncstart'], mapr: hdr['nrstart'], maps: hdr['nsstart']}
    Nxyz   = {1: hdr['nx'], 2: hdr['ny'], 3: hdr['nz']}
    axis_to_frac = {mapc: col_g, mapr: row_g, maps: sec_g}

    frac_x = (axis_to_frac[1].ravel() + start[1]) / Nxyz[1]
    frac_y = (axis_to_frac[2].ravel() + start[2]) / Nxyz[2]
    frac_z = (axis_to_frac[3].ravel() + start[3]) / Nxyz[3]
    frac_pts = np.stack([frac_x, frac_y, frac_z], axis=1)

    delta = eval_shift_field(frac_pts, hkls, AB_xyz)  # (N, 3) fractional
    cell = hdr['cell']
    delta_orth = frac_to_orth(delta, cell)             # (N, 3) Å

    for i, label in enumerate('xyz'):
        d = delta_orth[:, i].reshape(ns, nr, nc).astype(np.float32)
        write_ccp4(f"delta_{label}{nhkls}.map", d, hdr)
    dr = np.sqrt(np.sum(delta_orth**2, axis=1)).reshape(ns, nr, nc).astype(np.float32)
    write_ccp4(f"delta_r{nhkls}.map", dr, hdr)
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
        mtz.add_column(f'd{dim.upper()}', 'F')   # amplitude of d{x,y,z} shift (fractional Å)
        mtz.add_column(f'PH{dim.upper()}', 'P')  # phase in degrees
    mtz.add_column('SNR', 'R')
    mtz.add_column('ACTIVE', 'R')

    # Convert (A, B) Cartesian coefficients to (amplitude, phase in degrees).
    # Shift field: dx = Σ a*cos(2π(hx+ky+lz) + phi)
    # where A = -a*sin(phi_rad), B = a*cos(phi_rad)
    # so a = sqrt(A²+B²), phi_rad = atan2(-A, B)
    n_cols = 3 + 2 * n_dims + 2
    data = np.zeros((n_hkls, n_cols), dtype=np.float32)
    data[:, 0] = hkls[:, 0]
    data[:, 1] = hkls[:, 1]
    data[:, 2] = hkls[:, 2]
    for d in range(n_dims):
        A, B = AB[d, :, 0], AB[d, :, 1]
        data[:, 3 + 2*d]     = np.sqrt(A**2 + B**2)
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

    data = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]

    H = data[:, labels.index('H')].astype(int)
    K = data[:, labels.index('K')].astype(int)
    L = data[:, labels.index('L')].astype(int)
    hkls = np.stack([H, K, L], axis=1)

    n_hkls = len(hkls)
    n_dims = len(dimensions)
    AB = np.zeros((n_dims, n_hkls, 2))
    for d, dim in enumerate(dimensions):
        DIM = dim.upper()
        amp_col = (f'd{DIM}' if f'd{DIM}' in labels else
                   f'D{DIM}' if f'D{DIM}' in labels else
                   f'F{DIM}' if f'F{DIM}' in labels else None)
        if amp_col is not None:
            # Amplitude + phase in degrees
            a       = data[:, labels.index(amp_col)]
            phi_rad = np.radians(data[:, labels.index(f'PH{DIM}')])
            AB[d, :, 0] = -a * np.sin(phi_rad)
            AB[d, :, 1] =  a * np.cos(phi_rad)
        else:
            # Legacy Cartesian A, B coefficients
            AB[d, :, 0] = data[:, labels.index(f'A{DIM}')]
            AB[d, :, 1] = data[:, labels.index(f'B{DIM}')]

    snr    = data[:, labels.index('SNR')]
    active = data[:, labels.index('ACTIVE')].astype(bool)

    return hkls, AB, active, snr, cell1, cell2, dimensions, rmsd


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

    # Report largest CA fractional shift
    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}")
    if ca_mags[max_i] > 0.1:
        raise ValueError(f"Largest fractional CA shift {ca_mags[max_i]:.4f} > 0.1 cells. "
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
                         fitreso_start=20.0, fitreso_end=None,
                         drop_snr=1.0, batch_hkls=20, od_margin=1.5,
                         outlier_sigma=3.0, b_sigma=None,
                         use_symm=True, dimensions='xyz',
                         frac=1.0, verbose=True,
                         altloc_filter=False, altloc_fallback=1.0):
    """Fit a Fourier shift field progressively, admitting HKLs coarsest-first.

    HKLs are admitted in batches from fitreso_start (coarsest, large d) down to
    fitreso_end (finest, small d).  After each batch the fit residuals are
    MAD-filtered to adaptively update the active atom set for the next iteration.
    Stops when the overdetermination ratio (N_asu_active × 3) / (N_canon × 6)
    drops below od_margin, or fitreso_end is reached.

    Parameters
    ----------
    fitreso_start : float — d-spacing (Å) of the initial coarse batch; default 20 Å
    fitreso_end   : float — finest d-spacing to reach; default 2 Å
    batch_hkls    : int   — canonical HKLs to add per iteration
    od_margin     : float — stop before OD ratio drops below this (default 1.5)

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

    # Permanent initial outlier rejection (B-factor + shift MAD)
    if verbose:
        print("rejecting outliers...", end=' ', flush=True)
    fitme, ca_mask, uids, bfacs, _ = reject_outliers(
        fitme, ca_mask, uids, cell1, mad_sigma=outlier_sigma,
        bfacs=bfacs, b_sigma=b_sigma)
    if verbose:
        print(f"{ca_mask.sum()} CA remaining", flush=True)

    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}", flush=True)
    if ca_mags[max_i] > 0.1:
        raise ValueError(f"Largest fractional CA shift {ca_mags[max_i]:.4f} > 0.1 cells. "
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
    """Parse bendfinder.com-style key=value arguments."""
    pdb1 = pdb2 = mapfile = None
    params = {
        'nhkls':      30,
        'fitreso':    None,
        'drop_snr':   1.0,
        'frac':       1.0,
        'geotest':    False,
        'use_symm':   True,
        'dimensions': 'xyz',
        'fitparams':  None,
        'deltamaps':  False,
    }
    for arg in argv:
        if arg.endswith('.pdb') or arg.endswith('.PDB'):
            if pdb1 is None:
                pdb1 = arg
            elif pdb2 is None:
                pdb2 = arg
        elif arg.endswith('.map') or arg.endswith('.ccp4') or arg.endswith('.mrc'):
            if mapfile is None:
                mapfile = arg
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
        elif arg in ('deltamaps', 'delta', 'nofit'):
            if arg in ('deltamaps', 'delta'):
                params['deltamaps'] = True
            # nofit: ignored — lstsq always fits
    return pdb1, pdb2, mapfile, params


def main():
    argv = sys.argv[1:]
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)

    pdb1, pdb2, mapfile, p = _parse_args(argv)

    if pdb1 is None or pdb2 is None:
        print("ERROR: provide two PDB files on the command line.", file=sys.stderr)
        sys.exit(1)

    # Optionally load pre-computed fitparams
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

        # Save fitparams
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

    # Bent map
    if mapfile:
        if not os.path.exists(mapfile):
            print(f"WARNING: map file not found: {mapfile}", file=sys.stderr)
        else:
            bent_map = f"bent{nhkls}.map"
            make_bent_map(mapfile, hkls, AB_xyz, bent_map, cell2=cell2)
            if p['deltamaps']:
                make_delta_maps(mapfile, hkls, AB_xyz, nhkls)


if __name__ == '__main__':
    main()
