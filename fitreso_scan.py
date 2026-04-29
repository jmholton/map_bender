#!/programs/ccp4-8.0/bin/ccp4-python
"""Standard bendfinder fitreso-scan protocol.

Usage:
  fitreso_scan.py  MOV.pdb  MOV.map  REF.pdb  REF.map  [OPTIONS]
  fitreso_scan.py  MOV.pdb  MOV.mtz  REF.pdb  REF.mtz  [OPTIONS]

Map inputs may be CCP4 maps (.map/.ccp4) or MTZ files (.mtz); MTZ is
converted to a map grid via the F/PHI column pair (auto-detected or
specified with --fcol/--phicol).

Positional arguments:
  MOV.pdb   Moving PDB (will be bent onto REF frame)
  MOV.map   Moving CCP4 map or MTZ
  REF.pdb   Reference PDB
  REF.map   Reference CCP4 map or MTZ

Options:
  --outdir DIR          Output directory (default: ./fitreso_scan)
  --fullcell-mov        Read MOV map as full unit cell; avoids ASU-boundary
                        artefacts when both datasets share the same space group
                        and unit cell (e.g. humidity/temperature series)
  --fitreso LIST        Comma-separated fitreso_end values in Å
                        (default: 20,15,12,10,8,7,6,5)
  --chunk N             Voxels per chunk for shift-field evaluation;
                        0 = no chunking  (default: 50000)
  --fcol NAME           MTZ amplitude column (default: auto-detect)
  --phicol NAME         MTZ phase column    (default: auto-detect)
  --sample-rate FLOAT   Map oversampling when converting from MTZ (default: 3.0)

Outputs (per scan point under OUTDIR/hkl00..hkl10 and fr20..fr5):
  bent.map       Resampled moving map on reference grid
  diff_norm.map  Z-scored difference map (ref_norm − bent_norm)
  bent.mtz       Fbent/PHIbent + DELFWT/PHDELWT (for Coot DELFWT display)
  PSDVF.mtz      Fitted shift-field coefficients
  bent.pdb       Bent moving PDB
  ref.pdb        Reference PDB (copy)
  unbent.pdb     Unshifted moving PDB resampled on reference grid (copy)
  riso.log       Riso log (requires diff.com + CCP4 scaleit)

Summary line per scan point:
  label  RMSD0  RMSD  active  RFAC  peak+  atom+  peak-  atom-  [time]
"""
import sys, os, argparse, time, subprocess, contextlib, shutil, tempfile
import numpy as np, gemmi, math
from scipy import ndimage

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from bendfinder import (bend_fit_progressive, save_fitparams,
                        read_ccp4, read_ccp4_fullcell, write_ccp4,
                        eval_shift_field, interpolate_map,
                        expand_to_p1, match_atoms, reject_outliers,
                        rmsd_ca, write_bent_pdb)

# column name search order for auto-detection
_FCOL_TRY   = ['FWT',  '2FOFCWT', 'FC',   'F']
_PHICOL_TRY = ['PHWT', 'PH2FOFCWT', 'PHIC', 'PHI']

# ── argument parsing ───────────────────────────────────────────────────────────
def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('mov_pdb')
    p.add_argument('mov_map',  help='Moving CCP4 map or MTZ')
    p.add_argument('ref_pdb')
    p.add_argument('ref_map',  help='Reference CCP4 map or MTZ')
    p.add_argument('--outdir',       default='fitreso_scan')
    p.add_argument('--fullcell-mov', action='store_true')
    p.add_argument('--fitreso',      default='20,15,12,10,8,7,6,5')
    p.add_argument('--chunk',        type=int, default=50000)
    p.add_argument('--fcol',         default=None)
    p.add_argument('--phicol',       default=None)
    p.add_argument('--sample-rate',  type=float, default=3.0)
    return p.parse_args()

# ── map I/O ────────────────────────────────────────────────────────────────────
def _mtz_cols(mtz_path, fcol_hint, phicol_hint):
    mtz    = gemmi.read_mtz_file(mtz_path)
    labels = {c.label for c in mtz.columns}
    if fcol_hint and phicol_hint:
        if fcol_hint not in labels or phicol_hint not in labels:
            raise ValueError(f'MTZ columns {fcol_hint}/{phicol_hint} not found; '
                             f'available: {sorted(labels)}')
        return fcol_hint, phicol_hint
    for f, ph in zip(_FCOL_TRY, _PHICOL_TRY):
        if f in labels and ph in labels:
            return f, ph
    raise ValueError(f'Cannot find F/PHI columns in MTZ {mtz_path}; '
                     f'tried {list(zip(_FCOL_TRY, _PHICOL_TRY))}; '
                     f'available: {sorted(labels)}')

def _mtz_to_map(mtz_path, out_map, fcol, phicol, sample_rate):
    mtz   = gemmi.read_mtz_file(mtz_path)
    fc, ph = _mtz_cols(mtz_path, fcol, phicol)
    print(f'  MTZ columns: {fc}/{ph}', flush=True)
    grid  = mtz.transform_f_phi_to_map(fc, ph, sample_rate=sample_rate)
    ccp4  = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header(2, True)
    ccp4.write_ccp4_map(out_map)

def _load_map(path, outdir_cache, tag, fcol, phicol, sample_rate, fullcell):
    """Load a CCP4 map or MTZ into (data, header).

    If path is an MTZ, convert to a CCP4 map cached under outdir_cache so it
    can be re-read later (e.g. for map2mtz).  Returns (data, header, ccp4_path).
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == '.mtz':
        os.makedirs(outdir_cache, exist_ok=True)
        ccp4_path = os.path.join(outdir_cache, f'{tag}.map')
        if not os.path.exists(ccp4_path):
            print(f'  converting {tag} MTZ → map ({ccp4_path})...', flush=True)
            _mtz_to_map(path, ccp4_path, fcol, phicol, sample_rate)
        else:
            print(f'  re-using cached {tag} map', flush=True)
        path = ccp4_path
    else:
        ccp4_path = path

    if fullcell:
        d, h = read_ccp4_fullcell(path)
    else:
        d, h = read_ccp4(path)
    return d, h, ccp4_path

# ── map → MTZ (for ref.mtz used by diff.com) ──────────────────────────────────
def _map2mtz(mapfile, mtzfile):
    try:
        ccp4 = gemmi.read_ccp4_map(mapfile); ccp4.setup(float('nan'))
        grid = ccp4.grid
        sf   = gemmi.transform_map_to_f_phi(grid, half_l=True)
        data = sf.prepare_asu_data(dmin=0.0)
        mtz  = gemmi.Mtz(with_base=True)
        mtz.spacegroup = grid.spacegroup
        mtz.cell       = grid.unit_cell
        mtz.add_dataset('dataset')
        mtz.add_column('F',   'F')
        mtz.add_column('PHI', 'P')
        n   = len(data)
        arr = np.zeros((n, 5), dtype=np.float32)
        arr[:, :3] = data.miller_array
        arr[:,  3] = np.abs(data.value_array)
        arr[:,  4] = np.degrees(np.angle(data.value_array))
        mtz.set_data(arr)
        mtz.write_to_file(mtzfile)
        return True
    except Exception:
        return False

def _write_scan_mtz(bent_arr, diff_arr, ref_h, mtzfile):
    try:
        tmp_b = tempfile.NamedTemporaryFile(suffix='.map', delete=False).name
        tmp_d = tempfile.NamedTemporaryFile(suffix='.map', delete=False).name
        write_ccp4(tmp_b, bent_arr, ref_h)
        write_ccp4(tmp_d, diff_arr, ref_h)
        ccp4_b = gemmi.read_ccp4_map(tmp_b); ccp4_b.setup(float('nan'))
        ccp4_d = gemmi.read_ccp4_map(tmp_d); ccp4_d.setup(float('nan'))
        os.unlink(tmp_b); os.unlink(tmp_d)
        grid_b = ccp4_b.grid
        sf_b   = gemmi.transform_map_to_f_phi(grid_b,      half_l=True)
        data_b = sf_b.prepare_asu_data(dmin=0.0)
        sf_d   = gemmi.transform_map_to_f_phi(ccp4_d.grid, half_l=True)
        data_d = sf_d.prepare_asu_data(dmin=0.0)
        mtz = gemmi.Mtz(with_base=True)
        mtz.spacegroup = grid_b.spacegroup
        mtz.cell       = grid_b.unit_cell
        mtz.add_dataset('dataset')
        mtz.add_column('Fbent',   'F')
        mtz.add_column('PHIbent', 'P')
        mtz.add_column('DELFWT',  'F')
        mtz.add_column('PHDELWT', 'P')
        n   = len(data_b)
        arr = np.zeros((n, 7), dtype=np.float32)
        arr[:, :3] = data_b.miller_array
        arr[:,  3] = np.abs(data_b.value_array)
        arr[:,  4] = np.degrees(np.angle(data_b.value_array))
        arr[:,  5] = np.abs(data_d.value_array)
        arr[:,  6] = np.degrees(np.angle(data_d.value_array))
        mtz.set_data(arr)
        mtz.write_to_file(mtzfile)
        return True
    except Exception:
        return False

# ── Riso via diff.com ──────────────────────────────────────────────────────────
def _run_riso(outdir, ref_mtz, diff_com):
    test_mtz = os.path.join(outdir, 'bent.mtz')
    if not os.path.exists(test_mtz) or not os.path.exists(diff_com):
        return None
    logfile = os.path.join(outdir, 'riso.log')
    subprocess.run(
        ['tcsh', '-c',
         f'{diff_com} {ref_mtz} F {test_mtz} Fbent >& {logfile}'],
        cwd=outdir, env=dict(os.environ, DEBUG='1'))
    if os.path.exists(logfile):
        for line in open(logfile):
            if 'THE TOTALS' in line:
                parts = line.split()
                try:
                    return float(parts[6])
                except Exception:
                    pass
    return None

# ── geometry helpers ───────────────────────────────────────────────────────────
def _cell_matrix(hdr):
    a, b, c, al, be, ga = hdr['cell']
    al, be, ga = math.radians(al), math.radians(be), math.radians(ga)
    cos_al, cos_be, cos_ga = math.cos(al), math.cos(be), math.cos(ga)
    sin_ga = math.sin(ga)
    v = math.sqrt(1 - cos_al**2 - cos_be**2 - cos_ga**2 + 2*cos_al*cos_be*cos_ga)
    return np.array([
        [a,  b*cos_ga,  c*cos_be],
        [0,  b*sin_ga,  c*(cos_al - cos_be*cos_ga)/sin_ga],
        [0,  0,         c*v/sin_ga],
    ])

def _read_pdb_atoms_p1(pdb_path):
    st  = gemmi.read_pdb(pdb_path)
    sg  = st.find_spacegroup()
    ops = sg.operations() if sg else gemmi.SpaceGroup('P1').operations()
    out = []
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.has_altloc() and atom.altloc != 'A':
                        continue
                    f0 = st.cell.fractionalize(atom.pos)
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
                                    'frac': (R @ f0 + t) % 1.0})
    return out

def _nearest_atom(frac_target, atoms, M):
    arr  = np.array([a['frac'] for a in atoms])
    diff = arr - frac_target[np.newaxis, :]
    diff -= np.round(diff)
    dists = np.linalg.norm(diff @ M.T, axis=1)
    i = np.argmin(dists)
    return atoms[i], dists[i]

def _voxel_to_frac(idx_3d, hdr):
    ms, mr, mc = hdr['maps'], hdr['mapr'], hdr['mapc']
    nx, ny, nz = hdr['nx'],   hdr['ny'],   hdr['nz']
    ss, sr, sc = hdr['nsstart'], hdr['nrstart'], hdr['ncstart']
    is_, ir, ic = idx_3d
    grid  = {mc: ic, mr: ir, ms: is_}
    start = {mc: sc, mr: sr, ms: ss}
    N     = {1: nx, 2: ny, 3: nz}
    return np.array([(grid[ax] + start[ax]) / N[ax] for ax in [1, 2, 3]])

def _find_peaks(data, hdr, atoms, M, n=1, footprint=5):
    fp, sigma = np.ones((footprint,)*3), data.std()
    out = {}
    for sign, lbl in [(1, 'pos'), (-1, 'neg')]:
        d        = sign * data
        lmax     = d == ndimage.maximum_filter(d, footprint=fp)
        labeled, nf = ndimage.label(lmax & (d > 3*sigma))
        peaks = []
        for k in range(1, nf+1):
            mask = labeled == k
            idx  = np.unravel_index(np.argmax(d * mask), d.shape)
            peaks.append((sign * d[idx] / sigma, np.array(idx)))
        peaks.sort(key=lambda x: -abs(x[0]))
        top = []
        for val, idx in peaks[:n]:
            frac = _voxel_to_frac(idx, hdr)
            a, dist = _nearest_atom(frac, atoms, M)
            top.append((val, frac, a, dist))
        out[lbl] = top
    return out

def _atom_str(a):
    return f"{a['chain']}/{a['resseq']}{a['resname']}/{a['atomname']}"

_DUMMY = {'chain': '?', 'resseq': 0, 'resname': '?', 'atomname': '?'}

def _eval_chunked(ref_pts, hkls, AB, chunk):
    if chunk <= 0 or len(ref_pts) <= chunk:
        return eval_shift_field(ref_pts, hkls, AB)
    delta = np.zeros_like(ref_pts)
    for i in range(0, len(ref_pts), chunk):
        delta[i:i+chunk] = eval_shift_field(ref_pts[i:i+chunk], hkls, AB)
    return delta

# ── main ───────────────────────────────────────────────────────────────────────
def main():
    args = _parse_args()
    FITRESO_LIST = [int(x) for x in args.fitreso.split(',')]
    SCAN_DIR = os.path.abspath(args.outdir)
    MOV_PDB  = os.path.abspath(args.mov_pdb)
    REF_PDB  = os.path.abspath(args.ref_pdb)
    MOV_MAP  = os.path.abspath(args.mov_map)
    REF_MAP  = os.path.abspath(args.ref_map)
    CHUNK    = args.chunk
    DIFF_COM = os.path.join(_HERE, 'diff.com')

    print(f'MOV: {args.mov_pdb}  +  {args.mov_map}', flush=True)
    print(f'REF: {args.ref_pdb}  +  {args.ref_map}', flush=True)
    print(f'OUT: {SCAN_DIR}', flush=True)

    # ── load maps ──────────────────────────────────────────────────────────────
    os.makedirs(SCAN_DIR, exist_ok=True)
    print('Loading maps...', flush=True)
    mov_d, mov_h, _ = _load_map(MOV_MAP, SCAN_DIR, 'mov',
                                 args.fcol, args.phicol, args.sample_rate,
                                 fullcell=args.fullcell_mov)
    ref_d, ref_h, ref_ccp4_path = _load_map(REF_MAP, SCAN_DIR, 'ref',
                                              args.fcol, args.phicol, args.sample_rate,
                                              fullcell=False)
    if args.fullcell_mov:
        print('  moving map read as full unit cell', flush=True)

    # ── build reference grid fractional coordinates ────────────────────────────
    nc, nr, ns = ref_h['nc'], ref_h['nr'], ref_h['ns']
    nx, ny, nz = ref_h['nx'], ref_h['ny'], ref_h['nz']
    mc, mr, ms = ref_h['mapc'], ref_h['mapr'], ref_h['maps']
    g_s, g_r, g_c = np.meshgrid(np.arange(ns), np.arange(nr), np.arange(nc),
                                  indexing='ij')
    amap  = {mc: g_c, mr: g_r, ms: g_s}
    start = {mc: ref_h['ncstart'], mr: ref_h['nrstart'], ms: ref_h['nsstart']}
    Nxyz  = {1: nx, 2: ny, 3: nz}
    frac_x = (amap[1].ravel() + start[1]) / Nxyz[1]
    frac_y = (amap[2].ravel() + start[2]) / Nxyz[2]
    frac_z = (amap[3].ravel() + start[3]) / Nxyz[3]
    ref_pts = np.stack([frac_x, frac_y, frac_z], axis=1)
    n_vox = len(ref_pts)
    print(f'  ref grid: {ns}×{nr}×{nc} = {n_vox:,} voxels  '
          f'(chunk size: {CHUNK or "off"})', flush=True)

    M     = _cell_matrix(ref_h)
    atoms = _read_pdb_atoms_p1(REF_PDB)

    # ── ref.mtz for Riso (FT of reference density map) ────────────────────────
    ref_mtz_path = os.path.join(SCAN_DIR, 'ref.mtz')
    if not os.path.exists(ref_mtz_path):
        print('Converting ref map → ref.mtz...', flush=True)
        _map2mtz(ref_ccp4_path, ref_mtz_path)
    shutil.copy(REF_PDB, os.path.join(SCAN_DIR, 'ref.pdb'))
    _h0  = np.zeros((0, 3), dtype=int)
    _AB0 = np.zeros((3, 0, 2))
    write_bent_pdb(MOV_PDB, REF_PDB, _h0, _AB0,
                   os.path.join(SCAN_DIR, 'unbent.pdb'))

    # ── scan-point helper ──────────────────────────────────────────────────────
    def save_scan_point(label, outdir, bent_map, rmsd0, rmsd1,
                        n_active, t_elapsed, hkls=None, AB=None):
        write_ccp4(f'{outdir}/bent.map', bent_map, ref_h)
        if hkls is not None and AB is not None:
            write_bent_pdb(MOV_PDB, REF_PDB, hkls, AB, f'{outdir}/bent.pdb')
        shutil.copy(REF_PDB, f'{outdir}/ref.pdb')
        shutil.copy(os.path.join(SCAN_DIR, 'unbent.pdb'), f'{outdir}/unbent.pdb')
        ref_n  = (ref_d   - ref_d.mean())    / ref_d.std()
        bent_n = (bent_map - bent_map.mean()) / bent_map.std()
        diff   = ref_n - bent_n
        dnorm  = diff / diff.std()
        write_ccp4(f'{outdir}/diff_norm.map', dnorm, ref_h)
        _write_scan_mtz(bent_map, diff, ref_h, f'{outdir}/bent.mtz')
        rfac  = _run_riso(outdir, ref_mtz_path, DIFF_COM)
        peaks = _find_peaks(dnorm, ref_h, atoms, M, n=1)
        pp = peaks['pos'][0] if peaks['pos'] else (0, None, _DUMMY, 0)
        pn = peaks['neg'][0] if peaks['neg'] else (0, None, _DUMMY, 0)
        rfac_s  = f'{rfac*100:.1f}%' if rfac is not None else '  N/A'
        rmsd0_s = f'{rmsd0:.3f}'     if rmsd0 is not None else '  ---'
        rmsd1_s = f'{rmsd1:.3f}'     if rmsd1 is not None else '  ---'
        print(
            f"{label:>7}  {rmsd0_s:>6}  {rmsd1_s:>6}  {n_active:>6d}  {rfac_s:>6}  "
            f"{pp[0]:>+7.2f}σ  {_atom_str(pp[2]):>20} {pp[3]:.2f}Å  "
            f"{pn[0]:>+7.2f}σ  {_atom_str(pn[2]):>20} {pn[3]:.2f}Å  "
            f"[{t_elapsed:.0f}s]", flush=True)

    # ── column header ──────────────────────────────────────────────────────────
    print(f"\n{'label':>7}  {'RMSD0':>6}  {'RMSD':>6}  {'active':>6}  {'RFAC':>6}  "
          f"{'peak+':>7}  {'atom+':>20}  {'peak-':>7}  {'atom-':>20}", flush=True)
    print('-' * 122, flush=True)

    # ── Section 1: hkl00 — zero shift ─────────────────────────────────────────
    outdir0 = os.path.join(SCAN_DIR, 'hkl00')
    os.makedirs(outdir0, exist_ok=True)
    t0   = time.time()
    bent0 = interpolate_map(mov_d, mov_h, ref_pts).reshape(ns, nr, nc)
    with open(os.devnull, 'w') as _dev, contextlib.redirect_stdout(_dev):
        _a1, _cell1, _ = expand_to_p1(MOV_PDB)
        _a2, _cell2, _ = expand_to_p1(REF_PDB)
        _fm, _ca, _uid, _bf = match_atoms(_a1, _a2)
        _fm, _ca, *_ = reject_outliers(_fm, _ca, _uid, _cell1,
                                       mad_sigma=2.5, b_sigma=3.0, bfacs=_bf)
    raw_rmsd = rmsd_ca(_fm, _ca, _h0, _AB0, _cell2)
    save_scan_point('hkl00', outdir0, bent0, raw_rmsd, raw_rmsd,
                    0, time.time() - t0, hkls=_h0, AB=_AB0)

    # ── Section 2: hkl01..10 — one canonical HKL at a time ────────────────────
    def hkl_callback(iter_i, nhkls, n_canon, rmsd, hkls_now, AB_now, active, snr):
        n_non_dc = n_canon - 1
        if n_non_dc < 1 or n_non_dc > 10:
            return
        t0    = time.time()
        label = f'hkl{n_non_dc:02d}'
        odir  = os.path.join(SCAN_DIR, label)
        os.makedirs(odir, exist_ok=True)
        delta = _eval_chunked(ref_pts, hkls_now, AB_now, CHUNK)
        bent  = interpolate_map(mov_d, mov_h, ref_pts - delta).reshape(ns, nr, nc)
        save_scan_point(label, odir, bent, raw_rmsd, rmsd,
                        int(active.sum()), time.time() - t0,
                        hkls=hkls_now, AB=AB_now)
        save_fitparams(os.path.join(odir, 'PSDVF.mtz'),
                       hkls_now, AB_now, active, snr,
                       ref_h['cell'], ref_h['cell'], 'xyz', rmsd)

    bend_fit_progressive(
        MOV_PDB, REF_PDB,
        fitreso_start=100.0, fitreso_end=2.0,
        max_canon=11, batch_hkls=1, od_margin=1.5,
        drop_snr=0.0, outlier_sigma=2.5, b_sigma=3.0,
        verbose=False,
        iter_callback=hkl_callback,
    )

    # ── Section 3: fitreso scan ────────────────────────────────────────────────
    print('-' * 122, flush=True)
    for fitreso in FITRESO_LIST:
        odir = os.path.join(SCAN_DIR, f'fr{fitreso}')
        os.makedirs(odir, exist_ok=True)
        t0 = time.time()
        result = bend_fit_progressive(
            MOV_PDB, REF_PDB,
            fitreso_start=20.0, fitreso_end=float(fitreso),
            drop_snr=0.0, batch_hkls=100, od_margin=1.5,
            outlier_sigma=2.5, b_sigma=3.0,
            verbose=False,
        )
        t_fit = time.time() - t0
        delta = _eval_chunked(ref_pts, result.hkls, result.AB, CHUNK)
        bent  = interpolate_map(mov_d, mov_h, ref_pts - delta).reshape(ns, nr, nc)
        save_scan_point(f'fr{fitreso}', odir, bent, raw_rmsd, result.rmsd,
                        int(result.active.sum()), t_fit,
                        hkls=result.hkls, AB=result.AB)
        save_fitparams(f'{odir}/PSDVF.mtz',
                       result.hkls, result.AB, result.active, result.snr,
                       result.cell1, result.cell2, result.dimensions, result.rmsd)

    print('\nDone.', flush=True)


if __name__ == '__main__':
    main()
