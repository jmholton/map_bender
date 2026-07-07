#!/usr/bin/env ccp4-python
"""Embossing-artifact diagnostic for a bendfinder scan point.

For each strong diff-map peak, ask: is it sitting *on* the reference
2Fo-Fc peak, or on the *far side* of it with a matching negative peak
on the near side?  The latter is the classic dipole pattern of a
mis-registered atom — an alignment residual that survived the shift
field.

Inputs (in a fitreso_scan output dir like fitreso_scans/1-6-100kGy/best/):
  ref.mtz     — ref 2Fo-Fc (columns F, PHI)
  bent.mtz    — bent 2Fo-Fc (FDM,PHIDM) and diff (DELFWT,PHDELWT)
  ref.pdb     — reference atoms (for peak labelling)

Uses gemmi's NeighborSearch for symmetry-aware atom lookup, so peaks
near symmetry-equivalent atoms in a hexagonal / high-symmetry cell get
labelled correctly.

Usage:
  ccp4-python embossing_check.py <scan_point_dir> [n_top=20] [sigma_min=4.0]
"""
from __future__ import print_function
import sys
import os
import numpy as np

import gemmi
from scipy.ndimage import maximum_filter, minimum_filter, map_coordinates


# ── I/O ────────────────────────────────────────────────────────────────

def mtz_to_grid(mtz_path, F_col, P_col, sample_rate=3.0):
    """FFT F,PHI → numpy grid (NX,NY,NZ), σ-normalised."""
    mtz = gemmi.read_mtz_file(mtz_path)
    F = mtz.get_f_phi_on_grid(F_col, P_col,
                              size=mtz.get_size_for_hkl(sample_rate=sample_rate),
                              half_l=False)
    dens = gemmi.transform_f_phi_grid_to_map(F)
    arr  = np.array(dens.array, copy=True)
    rms  = float(np.sqrt(np.mean(arr**2)))
    if rms > 0:
        arr = arr / rms
    return arr, dens.unit_cell, (dens.nu, dens.nv, dens.nw)


def frac_to_orth(cell, frac):
    p = cell.orthogonalize(gemmi.Fractional(*frac))
    return np.array([p.x, p.y, p.z])


def orth_to_frac(cell, orth):
    f = cell.fractionalize(gemmi.Position(*orth))
    return np.array([f.x, f.y, f.z])


def index_to_frac(idx, dims):
    return np.array([idx[i] / dims[i] for i in range(3)])


def sample(arr, frac):
    """Cubic interpolation at fractional coord, crystal-wrap boundary."""
    dims = arr.shape
    idx  = np.array([frac[i] * dims[i] for i in range(3)])
    return float(map_coordinates(arr, idx.reshape(3, 1),
                                  order=3, mode='wrap')[0])


# ── Peak finding + ASU-dedup ──────────────────────────────────────────

def build_ns(structure, radius=6.0):
    """Symmetry-aware NeighborSearch — reads any (chain, resnum, atom)
    for a real-space position and its symmetry mates."""
    ns = gemmi.NeighborSearch(structure[0], structure.cell, radius).populate()
    return ns


def find_peaks_asu(arr, cell, dims, sigma_min, min_sep_A,
                    structure, positive=True):
    """Peak search, then keep only the ASU representative of each SG
    orbit.  Uses gemmi.find_asu_brick to filter by fractional bounds,
    or falls back to distance-based dedup vs a running list."""
    spacings = [cell.a / dims[0], cell.b / dims[1], cell.c / dims[2]]
    fp = tuple(max(3, int(np.ceil(min_sep_A / s)) | 1) for s in spacings)
    if positive:
        mx = maximum_filter(arr, size=fp, mode='wrap')
        mask = (arr == mx) & (arr >= sigma_min)
    else:
        mn = minimum_filter(arr, size=fp, mode='wrap')
        mask = (arr == mn) & (arr <= -sigma_min)
    idxs = np.argwhere(mask)
    vals = arr[mask]
    order = np.argsort(-np.abs(vals))
    # Build symop list for dedup
    sg  = structure.spacegroup_hm and gemmi.SpaceGroup(structure.spacegroup_hm)
    ops = sg.operations() if sg else None

    out = []
    seen_orth = []
    for k in order:
        i, j, m = idxs[k]
        frac  = index_to_frac((i, j, m), dims)
        orth  = frac_to_orth(cell, frac)
        # Skip if within 1.5 Å of a symmetry-equivalent already-kept peak
        is_dup = False
        if ops is not None:
            for prev_orth, prev_frac in seen_orth:
                for op in ops:
                    tf = op.apply_to_xyz(list(prev_frac))
                    tf = np.array([tf[0] % 1.0, tf[1] % 1.0, tf[2] % 1.0])
                    df = frac - tf
                    df -= np.round(df)
                    dd = frac_to_orth(cell, df)
                    if np.linalg.norm(dd) < 1.5:
                        is_dup = True
                        break
                if is_dup:
                    break
        if is_dup:
            continue
        seen_orth.append((orth, frac))
        out.append(dict(idx=(int(i), int(j), int(m)),
                        frac=frac, orth=orth, value=float(vals[k])))
    return out


def find_nearest_local_max(arr, seed_frac, cell, dims,
                            search_radius_A=3.0,
                            min_peak_sep_A=1.5):
    """Search around `seed_frac` for the nearest LOCAL MAXIMUM (a voxel
    whose value ≥ all neighbours) in `arr`.  Returns (frac_pos, value,
    dist_A) or (seed_frac, arr[seed], 0.0) if no local max in range.
    """
    spacings = [cell.a / dims[0], cell.b / dims[1], cell.c / dims[2]]
    r_vox    = [max(1, int(np.ceil(search_radius_A / s))) for s in spacings]
    seed_idx = np.array([seed_frac[i] * dims[i] for i in range(3)])

    # Extract a subgrid centred on the seed
    di = np.arange(-r_vox[0], r_vox[0] + 1)
    dj = np.arange(-r_vox[1], r_vox[1] + 1)
    dk = np.arange(-r_vox[2], r_vox[2] + 1)
    ii = (int(seed_idx[0]) + di) % dims[0]
    jj = (int(seed_idx[1]) + dj) % dims[1]
    kk = (int(seed_idx[2]) + dk) % dims[2]
    sub = arr[np.ix_(ii, jj, kk)]

    # Local maxima inside sub: value == 3x3x3 max-filter
    lm = (sub == maximum_filter(sub, size=3, mode='reflect'))
    lm_pts = np.argwhere(lm)
    if len(lm_pts) == 0:
        return seed_frac, float(arr[tuple(seed_idx.astype(int) % dims)]), 0.0

    # For each local-max candidate, compute distance to seed in Å
    best = None
    for lp in lm_pts:
        i2 = int(ii[lp[0]]); j2 = int(jj[lp[1]]); k2 = int(kk[lp[2]])
        lp_frac = np.array([i2/dims[0], j2/dims[1], k2/dims[2]])
        df = lp_frac - seed_frac
        df -= np.round(df)
        dd = frac_to_orth(cell, df)
        dist = float(np.linalg.norm(dd))
        val  = float(sub[tuple(lp)])
        if best is None or dist < best[2]:
            best = (lp_frac, val, dist)
    return best


# ── Symmetry-aware atom lookup ────────────────────────────────────────

def nearest_atom_symm(orth, ns, structure):
    """Nearest atom (any SG copy) to a cartesian position."""
    marks = ns.find_atoms(gemmi.Position(*orth), '\0', radius=8.0)
    if not marks:
        return None, None
    best = None
    for m in marks:
        image = m.pos                     # image = the symmetry copy position
        d = float(np.linalg.norm(np.array([image.x - orth[0],
                                            image.y - orth[1],
                                            image.z - orth[2]])))
        if best is None or d < best[1]:
            cra = m.to_cra(structure[0])
            lab = f"{cra.chain.name}/{cra.residue.seqid.num}{cra.residue.name}/{cra.atom.name}"
            best = (lab, d)
    return best


# ── Main ──────────────────────────────────────────────────────────────

def main(scan_dir, n_top=20, sigma_min=4.0):
    ref_mtz  = os.path.join(scan_dir, 'ref.mtz')
    bent_mtz = os.path.join(scan_dir, 'bent.mtz')
    ref_pdb  = os.path.join(scan_dir, 'ref.pdb')
    if not all(os.path.exists(p) for p in [ref_mtz, bent_mtz, ref_pdb]):
        print(f"missing input file in {scan_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"# embossing check for {scan_dir}")
    print(f"# top {n_top} SG-unique diff peaks by |σ|, sigma_min={sigma_min}")
    print()

    print("loading ref density from ref.mtz F/PHI ...")
    ref, cell, dims = mtz_to_grid(ref_mtz, 'F', 'PHI')
    print(f"  grid {dims}, cell {cell.a:.2f} {cell.b:.2f} {cell.c:.2f}")

    print("loading diff density from bent.mtz DELFWT/PHDELWT ...")
    diff, cell2, dims2 = mtz_to_grid(bent_mtz, 'DELFWT', 'PHDELWT')
    assert dims == dims2, f"grid mismatch: ref {dims} vs diff {dims2}"

    st = gemmi.read_structure(ref_pdb)
    st.setup_entities()
    ns = build_ns(st, radius=8.0)
    print(f"loaded {sum(1 for m in st for c in m for r in c for a in r)} atoms "
          f"in SG {st.spacegroup_hm}")

    print("finding peaks (SG-deduped) ...")
    pos = find_peaks_asu(diff, cell, dims, sigma_min, 2.0, st, positive=True)
    neg = find_peaks_asu(diff, cell, dims, sigma_min, 2.0, st, positive=False)
    all_peaks = sorted(pos + neg, key=lambda p: -abs(p['value']))[:n_top]
    print(f"kept {len(pos)} +peaks / {len(neg)} −peaks above ±{sigma_min}σ")
    print()

    print(f"{'rank':>4}  {'σ':>7}  {'atom':<24}  {'d_atom':>6}  "
          f"{'d_dens':>6}  {'mirror_σ':>9}  {'verdict':<11}  {'dip_ness':>8}")
    print('-' * 108)

    counts = dict(on_dens=0, dipole=0, weak_dipole=0, off=0)
    top_row = None
    # Continuous dipole-ness stats — σ²-weighted mean over all peaks in
    # the top list.  Per-peak formula:
    #   dip_ness = (-v_mirror / v_p) * min(1, d_dens / 0.5)
    # gate suppresses on-density peaks (d_dens ≈ 0) that would otherwise
    # trivially report mirror = self.
    d_score_num = 0.0          # Σ v² × dipole_ness
    d_score_den = 0.0          # Σ v²
    for i, p in enumerate(all_peaks, 1):
        lab_dist = nearest_atom_symm(p['orth'], ns, st)
        lab   = lab_dist[0] if lab_dist and lab_dist[0] else '(none)'
        d_at  = lab_dist[1] if lab_dist and lab_dist[1] is not None else 99.0

        q_frac, q_val, d_dens = find_nearest_local_max(
            ref, p['frac'], cell, dims, search_radius_A=3.0)

        p_prime = np.array([(2*q_frac[k] - p['frac'][k]) % 1.0
                             for k in range(3)])
        mirror_val = sample(diff, p_prime)
        ratio = -mirror_val / p['value']    # >0 if mirror has opposite sign

        if d_dens < 0.5:
            verdict = 'on_dens'
        elif ratio > 0.5 and 0.5 < d_dens < 2.5:
            verdict = 'dipole'
        elif ratio > 0.3 and 0.5 < d_dens < 2.5:
            verdict = 'weak_dipole'
        else:
            verdict = 'off'

        # Continuous dipole-ness (see comments above)
        gate = min(1.0, d_dens / 0.5)
        dip_ness = (-mirror_val / p['value']) * gate if p['value'] != 0 else 0.0
        w = p['value'] ** 2
        d_score_num += w * dip_ness
        d_score_den += w

        print(f"{i:4d}  {p['value']:+7.2f}  {lab:<24}  {d_at:6.2f}  "
              f"{d_dens:6.2f}  {mirror_val:+9.2f}  {verdict:<11}  {dip_ness:+.3f}")
        counts[verdict] += 1
        if i == 1:
            top_row = (p['value'], lab, d_at, d_dens, mirror_val, verdict)

    print()
    print(f"SUMMARY  n_on_dens={counts['on_dens']}  n_dipole={counts['dipole']}  "
          f"n_weak_dipole={counts['weak_dipole']}  n_off={counts['off']}  "
          f"n_top={sum(counts.values())}")
    d_score = d_score_num / d_score_den if d_score_den > 0 else 0.0
    print(f"DIPOLE   score={d_score:+.4f}  "
          f"(σ²-weighted mean of dipole_ness over all top peaks; "
          f"0=none, 1=perfect dipoles)")
    if top_row:
        print(f"TOP      σ={top_row[0]:+.2f}  atom={top_row[1]}  d_atom={top_row[2]:.2f}Å  "
              f"d_dens={top_row[3]:.2f}Å  mirror={top_row[4]:+.2f}σ  verdict={top_row[5]}")
    print()
    print("verdict key:")
    print("  on_dens      diff peak sits on the ref-density peak (< 0.5 Å) — real Δρ")
    print("  dipole       ref-density peak between + and − diff peaks (mirror ≥ 50% opposite)")
    print("  weak_dipole  30-50% opposite-sign mirror — suggestive")
    print("  off          no ref-density peak nearby, no partner — new density or noise")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    scan_dir  = sys.argv[1]
    n_top     = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    sigma_min = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
    main(scan_dir, n_top, sigma_min)
