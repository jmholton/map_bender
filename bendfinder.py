#!/programs/ccp4-8.0/bin/ccp4-python
"""bendfinder.py — fit a Fourier shift field between two crystal forms.

Usage:
    bendfinder.py moving.pdb reference.pdb [map.ccp4] [key=value ...]

Parameters:
    nhkls=30        max Fourier terms (default 30)
    reso=3          resolution cutoff in Angstroms (default 3)
    drop_snr=1      post-fit SNR threshold (default 1; 0 = disable)
    frac=1          scale factor for applied shifts (default 1)
    geotest=false   run refmac5 geometry check (default false)
    dimensions=xyz  which shifts to fit (subset of xyzoBà; default xyz)
    fitparams=file  load pre-computed fitparams.npy and skip fitting
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


def expand_to_p1(pdb_path):
    """Read PDB and expand ASU to P1 with all crystallographic symmetry ops.

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

    atoms = []
    for oi, op in enumerate(ops):
        for model in st:
            for chain in model:
                for res in chain:
                    if res.name in ('HOH', 'WAT', 'H2O'):
                        continue
                    for atom in res:
                        if atom.element == gemmi.Element('H'):
                            continue
                        if atom.altloc not in ('\x00', 'A', ' ', '\0'):
                            continue
                        orth = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                        frac = M_inv @ orth
                        nf = op.apply_to_xyz(frac.tolist())
                        nf = [f % 1.0 for f in nf]

                        aname = atom.name.strip()
                        rname = res.name.strip()
                        cname = chain.name.strip()
                        rnum  = res.seqid.num
                        icode = res.seqid.icode.strip()

                        uid = f"{cname}_{rnum}{icode}_{rname}_{aname}_op{oi}"
                        atoms.append({
                            'x': nf[0], 'y': nf[1], 'z': nf[2],
                            'bfac': atom.b_iso, 'occ': atom.occ,
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
    """
    idx2 = {a['uid']: a for a in atoms2}
    rows, ca, uids = [], [], []
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
    if not rows:
        raise ValueError("No atoms matched between the two structures — "
                         "check space groups match.")
    return np.array(rows, dtype=float), np.array(ca, dtype=bool), uids


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


def _fit_one_dim(Xa, b):
    """Fit one dimension: b ≈ Xa @ [A,B,...].

    Returns (AB array (M,2), snr array (M,)).
    """
    sol, res, rank, _ = sp_lstsq(Xa, b)
    n = len(b)
    dof = max(n - rank, 1)
    res_arr = np.atleast_1d(res)
    if res_arr.size > 0:
        sig2 = float(res_arr[0]) / dof
    else:
        sig2 = float(np.dot(b - Xa @ sol, b - Xa @ sol)) / dof

    # Covariance diagonal: sig2 * diag(V S⁻² Vᵀ) from thin SVD of Xa
    _, s_svd, Vt = np.linalg.svd(Xa, full_matrices=False)
    s_thresh = s_svd.max() * max(Xa.shape) * np.finfo(float).eps * 100
    s_inv = np.zeros_like(s_svd)
    nz = s_svd > s_thresh
    s_inv[nz] = 1.0 / s_svd[nz]
    cov_diag = sig2 * np.sum((Vt * s_inv[:, None])**2, axis=0)   # (p,)

    AB      = sol.reshape(-1, 2)           # (M, 2)
    cov_AB  = cov_diag.reshape(-1, 2)     # (M, 2)
    A, B    = AB[:, 0], AB[:, 1]
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
    active  = np.ones(N_hkls, dtype=bool)
    AB_full  = np.zeros((N_dims, N_hkls, 2))
    snr_full = np.zeros(N_hkls)

    for _ in range(max_rounds):
        col_mask = np.repeat(active, 2)
        Xa = X[:, col_mask]
        ai = np.where(active)[0]

        snr_sum = np.zeros(len(ai))
        for d in range(N_dims):
            AB_d, snr_d = _fit_one_dim(Xa, shifts[:, d])
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

    Overwrites DMIN/DMAX/DMEAN/RMS (words 20-23) with actual stats.
    cell_override: optional (a,b,c,al,be,ga) to replace cell in header.
    """
    raw = bytearray(hdr_template['raw_header'])

    mn, mx = float(data.min()), float(data.max())
    mean   = float(data.mean())
    rms    = float(np.sqrt(np.mean(data**2)))
    for off, val in zip((76, 80, 84, 88), (mn, mx, mean, rms)):
        struct.pack_into('<f', raw, off, val)

    if cell_override is not None:
        for i, val in enumerate(cell_override):
            struct.pack_into('<f', raw, 40 + i*4, val)

    with open(path, 'wb') as f:
        f.write(raw)
        f.write(data.astype(np.float32).tobytes())


def _frac_to_grid_indices(frac_xyz, hdr):
    """Convert fractional coords (N,3) → grid indices [sec, row, col] (3, N)."""
    nx, ny, nz = hdr['nx'], hdr['ny'], hdr['nz']
    ncstart, nrstart, nsstart = hdr['ncstart'], hdr['nrstart'], hdr['nsstart']
    mapc, mapr, maps = hdr['mapc'], hdr['mapr'], hdr['maps']

    # Grid coordinates for x, y, z axes
    gxyz = {
        1: frac_xyz[:, 0] * nx - ncstart,
        2: frac_xyz[:, 1] * ny - nrstart,
        3: frac_xyz[:, 2] * nz - nsstart,
    }
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

def make_bent_map(map_path, hkls, AB_xyz, outpath):
    """Resample map_path through the shift field and write outpath.

    For each grid point g in the OUTPUT map, we find where it came from
    in the INPUT map: sample at (g - δ(g)).  The negative sign is because
    the shift field maps the moving frame → reference frame; to resample
    we need the inverse.
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

    # Map column/row/section axes to x/y/z based on MAPC/MAPR/MAPS
    axis_to_frac = {mapc: col_g, mapr: row_g, maps: sec_g}
    N_grid = {mapc: nx, mapr: ny, maps: nz}
    start  = {mapc: ncstart, mapr: nrstart, maps: nsstart}

    frac_x = (axis_to_frac[1].ravel() + start[1]) / N_grid[1]
    frac_y = (axis_to_frac[2].ravel() + start[2]) / N_grid[2]
    frac_z = (axis_to_frac[3].ravel() + start[3]) / N_grid[3]
    frac_pts = np.stack([frac_x, frac_y, frac_z], axis=1)   # (N, 3)

    # Shift field at each grid point (fractional)
    delta = eval_shift_field(frac_pts, hkls, AB_xyz)   # (N, 3)

    # Sample the map at (grid_point - shift): negative because we invert
    source_frac = frac_pts - delta                      # (N, 3)
    new_vals = interpolate_map(data, hdr, source_frac)  # (N,)

    new_data = new_vals.reshape(ns, nr, nc)
    write_ccp4(outpath, new_data, hdr)
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
    N_grid = {hdr['mapc']: hdr['nx'], hdr['mapr']: hdr['ny'], hdr['maps']: hdr['nz']}
    start  = {hdr['mapc']: hdr['ncstart'], hdr['mapr']: hdr['nrstart'],
              hdr['maps']: hdr['nsstart']}
    axis_to_frac = {mapc: col_g, mapr: row_g, maps: sec_g}

    frac_x = (axis_to_frac[1].ravel() + start[1]) / N_grid[1]
    frac_y = (axis_to_frac[2].ravel() + start[2]) / N_grid[2]
    frac_z = (axis_to_frac[3].ravel() + start[3]) / N_grid[3]
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
# fitparams save / load (numpy format)
# ══════════════════════════════════════════════════════════════════════════════

def save_fitparams(path, hkls, AB, active, snr, cell1, cell2, dimensions, rmsd):
    # np.savez appends .npz automatically; strip .npz suffix first so the
    # caller's path is the actual filename on disk.
    path = str(path)
    if path.endswith('.npz'):
        path = path[:-4]
    np.savez(path,
             hkls=hkls, AB=AB, active=active, snr=snr,
             cell1=np.array(cell1), cell2=np.array(cell2),
             dimensions=np.array(list(dimensions)),
             rmsd=np.array(rmsd))
    return path + '.npz'


def load_fitparams(path):
    path = str(path)
    if not path.endswith('.npz'):
        path = path + '.npz'
    d = np.load(path, allow_pickle=False)
    return (d['hkls'], d['AB'], d['active'].astype(bool), d['snr'],
            tuple(d['cell1']), tuple(d['cell2']),
            ''.join(d['dimensions']), float(d['rmsd']))


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


def bend_fit(pdb1_path, pdb2_path, nhkls=30, reso=3.0, drop_snr=1.0,
             dimensions='xyz', geotest=False, frac=1.0, verbose=True):
    """Fit a Fourier shift field between pdb1 (moving) and pdb2 (reference).

    Returns BendResult with fitted coefficients and RMSD.
    """
    t0 = time.time()

    if verbose:
        print(f"expanding {pdb1_path} from", end=' ', flush=True)
    atoms1, cell1, sg1 = expand_to_p1(pdb1_path)
    if verbose:
        print(f"{sg1} to P1")
        print(f"expanding {pdb2_path} from", end=' ', flush=True)
    atoms2, cell2, sg2 = expand_to_p1(pdb2_path)
    if verbose:
        print(f"{sg2} to P1")

    fitme, ca_mask, uids = match_atoms(atoms1, atoms2)

    # Report largest CA fractional shift
    ca_shifts = fitme[ca_mask, 3:6]
    ca_mags   = np.sqrt(np.sum(ca_shifts**2, axis=1))
    ca_ids    = [u for u, m in zip(uids, ca_mask) if m]
    max_i     = int(np.argmax(ca_mags))
    print(f"largest fractional CA shift: {ca_mags[max_i]:.7f} for {ca_ids[max_i]}")
    if ca_mags[max_i] > 0.1:
        raise ValueError(f"Largest fractional CA shift {ca_mags[max_i]:.4f} > 0.1 cells. "
                         "Check structure alignment.")

    if verbose:
        print(f"generating basis functions out to {reso} A resolution")
    all_hkls = generate_hkls(cell1, reso)
    ntotal   = len(all_hkls)
    print(f"{ntotal - 1}  complex coefficients for each dimension: {dimensions}")

    nhkls_use = min(nhkls, ntotal)
    hkls      = all_hkls[:nhkls_use]

    # Dimension mapping: x=col3, y=col4, z=col5, o=col6, B=col7
    dim_col = {'x': 3, 'y': 4, 'z': 5, 'o': 6, 'B': 7}
    dim_list = [d for d in dimensions if d in dim_col]
    dim_indices = [dim_col[d] for d in dim_list]
    shifts = fitme[:, dim_indices]   # (N_atoms, N_dims)

    frac_coords = fitme[:, :3]
    X = build_design_matrix(frac_coords, hkls)

    dt = time.time() - t0
    print(f"\n=================  Starting nhkls {nhkls_use} fit at {dt:.0f} s\n")

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

    make_bent_map(map_path, result.hkls, AB_xyz, outpath)
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
        'reso':       3.0,
        'drop_snr':   1.0,
        'frac':       1.0,
        'geotest':    False,
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
            elif key == 'reso':
                params['reso'] = float(val)
            elif key == 'drop_snr':
                params['drop_snr'] = float(val)
            elif key == 'frac':
                params['frac'] = float(val)
            elif key == 'geotest':
                params['geotest'] = val.lower() in ('true', '1', 'yes')
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
            nhkls=p['nhkls'], reso=p['reso'], drop_snr=p['drop_snr'],
            dimensions=p['dimensions'], geotest=p['geotest'], frac=p['frac'],
        )
        nhkls = int(result.active.sum())

        # Save fitparams
        fp_path = save_fitparams(
            f"fitparams{nhkls}.npz", result.hkls, result.AB, result.active, result.snr,
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
            make_bent_map(mapfile, hkls, AB_xyz, bent_map)
            if p['deltamaps']:
                make_delta_maps(mapfile, hkls, AB_xyz, nhkls)


if __name__ == '__main__':
    main()
