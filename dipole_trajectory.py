#!/usr/bin/env ccp4-python
"""Single-process dipole-trajectory scan across all fitreso_scan points
of a chosen ref (default 1-6-100kGy).  Ref density is FFT'd once and
reused across every scan point; each point re-FFTs only its own diff
map.  Emits one row per scan point with:

    top_σ  n_on n_dip n_wk n_off  dipole_score  top_atom  top_verdict

The dipole_score is a σ²-weighted mean of per-peak dipole-ness in
~[-1,+1]:
    dipole_ness = (-v_mirror / v_p) * min(1, d_dens / 0.5)
where the gate suppresses on-density peaks (d_dens ≈ 0) that would
otherwise trivially report mirror = self.

Usage:
    ccp4-python dipole_trajectory.py [ref_name=1-6-100kGy] [top=30] [sigma=4.0]
"""
import os
import sys
import time
import numpy as np

import gemmi
from scipy.ndimage import maximum_filter, minimum_filter, map_coordinates


# ── helpers (copied so we can reuse ref grid without re-loading) ──────

def frac_to_orth(cell, frac):
    p = cell.orthogonalize(gemmi.Fractional(*frac))
    return np.array([p.x, p.y, p.z])


def index_to_frac(idx, dims):
    return np.array([idx[i] / dims[i] for i in range(3)])


def sample(arr, frac):
    dims = arr.shape
    idx  = np.array([frac[i] * dims[i] for i in range(3)])
    return float(map_coordinates(arr, idx.reshape(3, 1),
                                  order=3, mode='wrap')[0])


def _pick_f_phi(mtz_path, candidates):
    """Return the first (F_col, P_col) tuple whose labels are present
    in `mtz_path`.  Raises if none match."""
    mtz = gemmi.read_mtz_file(mtz_path)
    labels = {c.label for c in mtz.columns}
    for f, p in candidates:
        if f in labels and p in labels:
            return f, p
    raise RuntimeError(f"None of {candidates} found in {mtz_path} "
                       f"(have {sorted(labels)})")


def _mtz_to_grid(mtz_path, F_col, P_col, sample_rate=1.5):
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


def find_peaks_asu(arr, cell, dims, sigma_min, min_sep_A,
                    sg_ops, positive=True, max_candidates=200,
                    max_out=60):
    """SG-deduped local extrema.  Caps the candidate pool at
    `max_candidates` (highest |σ| first) so we don't walk thousands of
    raw peaks; caps kept peaks at `max_out` (still >> anything the
    caller will use)."""
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
    order = np.argsort(-np.abs(vals))[:max_candidates]   # cap here
    out = []
    seen = []           # list of (frac,) tuples of kept peaks
    for k in order:
        if len(out) >= max_out:
            break
        i, j, m = idxs[k]
        frac  = index_to_frac((i, j, m), dims)
        orth  = frac_to_orth(cell, frac)
        is_dup = False
        for prev_frac in seen:
            for op in sg_ops:
                tf = op.apply_to_xyz(list(prev_frac))
                tf = np.array([tf[0] % 1.0, tf[1] % 1.0, tf[2] % 1.0])
                df = frac - tf
                df -= np.round(df)
                if np.linalg.norm(frac_to_orth(cell, df)) < 1.5:
                    is_dup = True
                    break
            if is_dup:
                break
        if is_dup:
            continue
        seen.append(frac)
        out.append(dict(idx=(int(i), int(j), int(m)),
                        frac=frac, orth=orth, value=float(vals[k])))
    return out


def find_nearest_local_max(arr, seed_frac, cell, dims, search_radius_A=3.0):
    spacings = [cell.a / dims[0], cell.b / dims[1], cell.c / dims[2]]
    r_vox    = [max(1, int(np.ceil(search_radius_A / s))) for s in spacings]
    seed_idx = np.array([seed_frac[i] * dims[i] for i in range(3)])

    di = np.arange(-r_vox[0], r_vox[0] + 1)
    dj = np.arange(-r_vox[1], r_vox[1] + 1)
    dk = np.arange(-r_vox[2], r_vox[2] + 1)
    ii = (int(seed_idx[0]) + di) % dims[0]
    jj = (int(seed_idx[1]) + dj) % dims[1]
    kk = (int(seed_idx[2]) + dk) % dims[2]
    sub = arr[np.ix_(ii, jj, kk)]

    lm = (sub == maximum_filter(sub, size=3, mode='reflect'))
    lm_pts = np.argwhere(lm)
    if len(lm_pts) == 0:
        return seed_frac, float(arr[tuple(seed_idx.astype(int) % dims)]), 0.0

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


def nearest_atom_symm(orth, ns, structure):
    marks = ns.find_atoms(gemmi.Position(*orth), '\0', radius=8.0)
    if not marks:
        return None, None
    best = None
    for m in marks:
        image = m.pos
        d = float(np.linalg.norm(np.array([image.x - orth[0],
                                            image.y - orth[1],
                                            image.z - orth[2]])))
        if best is None or d < best[1]:
            cra = m.to_cra(structure[0])
            lab = f"{cra.chain.name}/{cra.residue.seqid.num}{cra.residue.name}/{cra.atom.name}"
            best = (lab, d)
    return best


# ── main ──────────────────────────────────────────────────────────────

def analyse_point(diff, cell, dims, sg_ops, ns, st, ref, n_top, sigma_min):
    """Return dict of per-point results — one call per scan point."""
    t0 = time.time()
    pos = find_peaks_asu(diff, cell, dims, sigma_min, 2.0, sg_ops, positive=True)
    t_pos = time.time() - t0
    t1 = time.time()
    neg = find_peaks_asu(diff, cell, dims, sigma_min, 2.0, sg_ops, positive=False)
    t_neg = time.time() - t1
    all_peaks = sorted(pos + neg, key=lambda p: -abs(p['value']))[:n_top]
    t_peaks = t_pos + t_neg
    print(f"  peak search: pos={t_pos:.0f}s (n={len(pos)}) neg={t_neg:.0f}s (n={len(neg)}) ", end='', flush=True)
    t_loop = time.time()

    counts = dict(on_dens=0, dipole=0, weak_dipole=0, off=0)
    d_num = d_den = 0.0
    top_row = None

    for i, p in enumerate(all_peaks, 1):
        _lab, _dat = nearest_atom_symm(p['orth'], ns, st) or (None, None)
        _lab = _lab if _lab else '(none)'
        _dat = _dat if _dat is not None else 99.0

        q_frac, _q_val, d_dens = find_nearest_local_max(
            ref, p['frac'], cell, dims, search_radius_A=3.0)
        p_prime = np.array([(2*q_frac[k] - p['frac'][k]) % 1.0
                             for k in range(3)])
        mirror_val = sample(diff, p_prime)
        ratio = -mirror_val / p['value'] if p['value'] != 0 else 0.0

        if d_dens < 0.5:
            verdict = 'on_dens'
        elif ratio > 0.5 and 0.5 < d_dens < 2.5:
            verdict = 'dipole'
        elif ratio > 0.3 and 0.5 < d_dens < 2.5:
            verdict = 'weak_dipole'
        else:
            verdict = 'off'
        counts[verdict] += 1

        gate = min(1.0, d_dens / 0.5)
        dip_ness = ratio * gate
        w = p['value'] ** 2
        d_num += w * dip_ness
        d_den += w

        if i == 1:
            top_row = dict(sig=p['value'], atom=_lab, d_atom=_dat,
                           d_dens=d_dens, mirror=mirror_val, verdict=verdict)

    print(f"loop={time.time()-t_loop:.0f}s ", end='', flush=True)
    score = d_num / d_den if d_den > 0 else 0.0
    return dict(counts=counts, d_score=score, top=top_row)


def main():
    # Positional: <scan_dir_or_ref_name> [n_top] [sigma_min] [out_log]
    # If arg1 is an existing directory, use it directly.
    # Otherwise treat as a name under JJD95 fitreso_scans/ (legacy default).
    arg1 = sys.argv[1] if len(sys.argv) > 1 else '1-6-100kGy'
    n_top    = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    sigma    = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
    out_log  = sys.argv[4] if len(sys.argv) > 4 else None

    if os.path.isdir(arg1):
        base = arg1
        ref_name = os.path.basename(base.rstrip('/'))
    else:
        ref_name = arg1
        base = ('/home/jamesh/projects/map_bender/arwen/JJD95_dose_series_mtzs/'
                f'fitreso_scans/{ref_name}')
    if out_log is None:
        out_log = ('/home/jamesh/projects/map_bender/arwen/JJD95_dose_series_mtzs/'
                   f'dipole_trajectory_{ref_name}.log')

    points = ['pre','hkl00','hkl01','hkl02','hkl03','hkl04','hkl05','hkl06',
              'hkl07','hkl08','hkl09','hkl10','fr20','fr15','fr12','fr10',
              'fr8','fr7','fr6','fr5','best']

    print(f"ref={ref_name}  base={base}")
    print(f"n_top={n_top}  sigma_min={sigma}")
    t0 = time.time()

    # Load ref density + structure once (shared across all points)
    print("loading ref density (once) ...", flush=True)
    ref_mtz_path = f'{base}/best/ref.mtz'
    ref_cols = _pick_f_phi(ref_mtz_path,
                            candidates=[('F', 'PHI'),
                                        ('FWT', 'PHWT'),
                                        ('FP', 'PHIC'),
                                        ('FC', 'PHIC')])
    print(f"  ref cols: {ref_cols[0]}/{ref_cols[1]}", flush=True)
    ref, cell, dims = _mtz_to_grid(ref_mtz_path, *ref_cols)
    print(f"  ref grid {dims}, cell {cell.a:.2f} {cell.b:.2f} {cell.c:.2f}",
          flush=True)

    print("loading ref structure + symops (once) ...", flush=True)
    st = gemmi.read_structure(f'{base}/best/ref.pdb')
    st.setup_entities()
    ns = gemmi.NeighborSearch(st[0], cell, 8.0).populate()
    sg_ops = gemmi.SpaceGroup(st.spacegroup_hm).operations()
    print(f"  {sum(1 for m in st for c in m for r in c for a in r)} atoms, "
          f"SG={st.spacegroup_hm}, {len(list(sg_ops))} operations", flush=True)

    rows = []
    for pt in points:
        pdir = f'{base}/{pt}'
        if not os.path.isdir(pdir):
            continue
        bent_mtz = f'{pdir}/bent.mtz'
        if not os.path.exists(bent_mtz):
            continue
        t = time.time()
        print(f">> {pt} ... ", end='', flush=True)
        try:
            t_fft = time.time()
            diff, _cell2, dims2 = _mtz_to_grid(bent_mtz, 'DELFWT', 'PHDELWT')
            dt_fft = time.time() - t_fft
            if dims2 != dims:
                # This shouldn't happen — same ref cell — but bail cleanly
                print(f"grid mismatch {dims2} vs {dims}, skipping")
                continue
            t_ana = time.time()
            res = analyse_point(diff, cell, dims, list(sg_ops), ns, st, ref,
                                 n_top, sigma)
            dt_ana = time.time() - t_ana
            rows.append((pt, res))
            top = res['top']
            c = res['counts']
            print(f"top {top['sig']:+.2f}σ  score={res['d_score']:+.4f}  "
                  f"n_on/dip/wk/off = {c['on_dens']}/{c['dipole']}/"
                  f"{c['weak_dipole']}/{c['off']}  "
                  f"(fft={dt_fft:.0f}s  ana={dt_ana:.0f}s  tot={time.time()-t:.0f}s)",
                  flush=True)
        except Exception as e:
            print(f"FAILED: {e}", flush=True)

    # Write final table
    with open(out_log, 'w') as fh:
        fh.write(f"# dipole trajectory for {ref_name}\n")
        fh.write(f"# top {n_top} peaks per scan point, |σ| ≥ {sigma}\n")
        fh.write(f"# dipole_ness = (-v_mirror/v_peak) * min(1, d_dens/0.5)\n")
        fh.write(f"# score       = Σ(v_peak² · dipole_ness) / Σ(v_peak²)\n")
        fh.write(f"# 0 = no dipole content;  1 = perfect dipoles;  "
                 f"< 0 = same-sign mirrors (rare)\n\n")
        fh.write(f"{'point':<7} {'top_σ':>7} {'dipole_score':>13}  "
                 f"{'on':>3} {'dip':>4} {'wk':>3} {'off':>4}  "
                 f"{'top_atom':<20} {'top_verdict':<11}\n")
        fh.write('-' * 90 + '\n')
        for pt, res in rows:
            top = res['top']
            c = res['counts']
            fh.write(f"{pt:<7} {top['sig']:>+7.2f} {res['d_score']:>+13.4f}  "
                     f"{c['on_dens']:>3} {c['dipole']:>4} "
                     f"{c['weak_dipole']:>3} {c['off']:>4}  "
                     f"{top['atom']:<20} {top['verdict']:<11}\n")
    print(f"\nWrote {out_log}")
    print(f"Total wall time: {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
