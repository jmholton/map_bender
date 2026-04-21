#!/programs/ccp4-8.0/bin/ccp4-python
"""Standard protocol: zero-cycle phenix.refine → 2Fo-Fc maps → bendfinder → CC comparison.

Usage:
  cd dhfr/
  srun /programs/ccp4-8.0/bin/ccp4-python run_test.py [pdb1 mtz1 pdb2 mtz2] [nhkls=N|fitreso=X]

Defaults to 1rx1 (moving) vs 1rx2 (reference).
Output: bent<N>.pdb, 1rx1_2fofc.map, 1rx2_2fofc.map, bent<N>.map, psdvf.mtz
"""

import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'map_bender'))
from bendfinder import bend_fit, read_ccp4, write_ccp4, save_fitparams, eval_shift_field, interpolate_map

PHENIX   = '/programs/phenix-2.0-5936/bin/phenix.refine'
CCP4_FFT = '/programs/ccp4-8.0/bin/fft'
CCP4_PY  = '/programs/ccp4-8.0/bin/ccp4-python'
USE_SRUN = True   # set False to run inline (for debugging)

# ── Parse arguments ────────────────────────────────────────────────────────────
pdb1 = '1rx1.pdb'; mtz1 = '1rx1.mtz'
pdb2 = '1rx2.pdb'; mtz2 = '1rx2.mtz'
nhkls = 100
fitreso = None
for arg in sys.argv[1:]:
    if arg.endswith('.pdb'):
        if '1' in arg or pdb1 == '1rx1.pdb':
            pdb1 = arg
        else:
            pdb2 = arg
    elif arg.startswith('nhkls='):
        nhkls = int(arg.split('=')[1])
    elif arg.startswith('fitreso='):
        fitreso = float(arg.split('=')[1])
    elif arg.startswith('pdb1='):
        pdb1 = arg.split('=', 1)[1]
    elif arg.startswith('pdb2='):
        pdb2 = arg.split('=', 1)[1]
    elif arg.startswith('mtz1='):
        mtz1 = arg.split('=', 1)[1]
    elif arg.startswith('mtz2='):
        mtz2 = arg.split('=', 1)[1]

stem1 = pdb1.replace('.pdb', '').replace('.PDB', '')
stem2 = pdb2.replace('.pdb', '').replace('.PDB', '')


def run(cmd, stdin=None, logfile=None):
    """Run a command, optionally via srun, and check return code."""
    if USE_SRUN:
        cmd = ['srun', '--ntasks=1'] + cmd
    print(f"  running: {' '.join(cmd[:4])} ...")
    result = subprocess.run(cmd, input=stdin, capture_output=True, text=True)
    if logfile:
        with open(logfile, 'w') as f:
            f.write(result.stdout)
            f.write(result.stderr)
    if result.returncode != 0:
        print(f"  FAILED (exit {result.returncode})")
        print(result.stderr[-2000:])
        sys.exit(1)
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: zero-cycle phenix.refine
# ══════════════════════════════════════════════════════════════════════════════
print(f'\n── Step 1: zero-cycle phenix.refine ──────────────────────────────────')

for pdb, mtz, stem in [(pdb1, mtz1, stem1), (pdb2, mtz2, stem2)]:
    refined_mtz = f'{stem}_refine_001.mtz'
    refined_pdb = f'{stem}_refine_001.pdb'
    if os.path.exists(refined_mtz) and os.path.exists(refined_pdb):
        print(f"  {refined_mtz} exists — skipping refinement")
        continue
    print(f"  refining {pdb} vs {mtz}")
    run([PHENIX, pdb, mtz,
         'refinement.main.number_of_macro_cycles=0',
         'refinement.main.nqh_flips=False',
         f'output.prefix={stem}_refine',
         'output.write_maps=False',
         ],
        logfile=f'{stem}_refine.log')
    if not os.path.exists(refined_mtz):
        print(f"ERROR: phenix did not produce {refined_mtz}")
        sys.exit(1)
    print(f"  → {refined_pdb}, {refined_mtz}")


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: make 2Fo-Fc maps via CCP4 FFT
# ══════════════════════════════════════════════════════════════════════════════
print(f'\n── Step 2: 2Fo-Fc maps (CCP4 FFT) ───────────────────────────────────')

CCP4_MAPMASK = '/programs/ccp4-8.0/bin/mapmask'

def make_2fofc_map(mtz_in, map_out):
    if os.path.exists(map_out):
        print(f"  {map_out} exists — skipping")
        return
    print(f"  FFT {mtz_in} → {map_out}")
    tmp = map_out.replace('.map', '_full.map')
    run([CCP4_FFT, 'hklin', mtz_in, 'mapout', tmp],
        stdin='labin F1=2FOFCWT PHI=PH2FOFCWT\ngrid sample 4.0\nend\n',
        logfile=map_out.replace('.map', '_fft.log'))
    # Canonical axis order and trim to ASU
    run([CCP4_MAPMASK, 'mapin', tmp, 'mapout', map_out],
        stdin='XYZLIM ASU\nAXIS X Y Z\nEND\n',
        logfile=map_out.replace('.map', '_mapmask.log'))
    os.remove(tmp)
    print(f"  → {map_out}")

make_2fofc_map(f'{stem1}_refine_001.mtz', f'{stem1}_2fofc.map')
make_2fofc_map(f'{stem2}_refine_001.mtz', f'{stem2}_2fofc.map')


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: fit shift field with bendfinder
# ══════════════════════════════════════════════════════════════════════════════
print(f'\n── Step 3: bendfinder shift-field fit ────────────────────────────────')
print(f"  moving={pdb1}  reference={pdb2}  nhkls={nhkls}")

result = bend_fit(pdb1, pdb2, nhkls=nhkls, fitreso=fitreso, reso=3.0, drop_snr=1.0, use_symm=True)
n_active = int(result.active.sum())
print(f"  RMSD(CA) = {result.rmsd:.3f} Å   active HKLs = {n_active}")

fp_path = save_fitparams('psdvf', result.hkls, result.AB, result.active, result.snr,
                         result.cell1, result.cell2, result.dimensions, result.rmsd)
print(f"  fitparams → {fp_path}")

# Build AB_xyz (3, N, 2) for xyz dims
dim_list = list(result.dimensions)
xyz_dims = [dim_list.index(d) for d in 'xyz' if d in dim_list]
AB_xyz = np.zeros((3, len(result.hkls), 2))
for new_i, old_i in enumerate(xyz_dims):
    AB_xyz[new_i] = result.AB[old_i]


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 & 5: apply shift on the REFERENCE grid, then compare
# ══════════════════════════════════════════════════════════════════════════════
# The correct way: for each reference grid point g_ref, find the source point
# in the moving frame (x_mov = g_ref - δ(g_ref)), interpolate the moving map
# there, and compare to the reference value at g_ref.
print(f'\n── Step 4: sample shifted {stem1} map on {stem2} grid ───────────────')

map_mov = f'{stem1}_2fofc.map'
map_ref = f'{stem2}_2fofc.map'

mov_d, mov_h = read_ccp4(map_mov)
ref_d, ref_h = read_ccp4(map_ref)

# Build fractional coordinates for every voxel in the REFERENCE grid
nc2, nr2, ns2 = ref_h['nc'], ref_h['nr'], ref_h['ns']
nx2, ny2, nz2 = ref_h['nx'], ref_h['ny'], ref_h['nz']
mapc2, mapr2, maps2 = ref_h['mapc'], ref_h['mapr'], ref_h['maps']
ncstart2, nrstart2, nsstart2 = ref_h['ncstart'], ref_h['nrstart'], ref_h['nsstart']

sec2, row2, col2 = np.meshgrid(np.arange(ns2), np.arange(nr2), np.arange(nc2), indexing='ij')
axis_to_frac2 = {mapc2: col2, mapr2: row2, maps2: sec2}
start2  = {mapc2: ncstart2, mapr2: nrstart2, maps2: nsstart2}
Nxyz2   = {1: nx2, 2: ny2, 3: nz2}

frac_x2 = (axis_to_frac2[1].ravel() + start2[1]) / Nxyz2[1]
frac_y2 = (axis_to_frac2[2].ravel() + start2[2]) / Nxyz2[2]
frac_z2 = (axis_to_frac2[3].ravel() + start2[3]) / Nxyz2[3]
ref_pts = np.stack([frac_x2, frac_y2, frac_z2], axis=1)   # (N, 3) in ref-frame fractional

# Evaluate shift field at reference grid points, find source in moving frame
delta = eval_shift_field(ref_pts, result.hkls, AB_xyz)   # (N, 3)
source_pts = ref_pts - delta                              # subtract δ to go backwards

# Interpolate the moving map at source points
bent_vals = interpolate_map(mov_d, mov_h, source_pts)     # (N,)
bent_data = bent_vals.reshape(ns2, nr2, nc2)

map_out = f'bent{n_active}.map'
write_ccp4(map_out, bent_data, ref_h)   # output on reference grid/cell
print(f"  → {map_out}  shape=({ns2},{nr2},{nc2})")


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: compare maps
# ══════════════════════════════════════════════════════════════════════════════
print(f'\n── Step 5: map comparison ────────────────────────────────────────────')

# Also compute the baseline (no-transform) CC: interpolate moving map at ref points with δ=0
noshift_vals = interpolate_map(mov_d, mov_h, ref_pts)
noshift_data = noshift_vals.reshape(ns2, nr2, nc2)

def cc(a, b):
    a = a.ravel() - a.mean()
    b = b.ravel() - b.mean()
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

ref_flat  = ref_d.ravel()
bent_flat = bent_vals
movp_flat = noshift_vals

cc_bent = cc(bent_flat.reshape(ref_d.shape), ref_d)
cc_mov  = cc(noshift_data, ref_d)

diff = bent_flat - ref_d.ravel()
sig  = ref_d.std()

print(f"  shapes:  ref={ref_d.shape}  bent={bent_data.shape}")
print(f"\n  CC(bent_{stem1}, {stem2}) = {cc_bent:.6f}   [after shift field]")
print(f"  CC({stem1},      {stem2}) = {cc_mov:.6f}   [baseline: no transform]")
print(f"\n  Diff sigma / ref sigma:  {diff.std()/sig*100:.2f}%")
print(f"  ref sigma = {sig:.5f}")
print(f"\nDone.  RMSD(CA)={result.rmsd:.3f} Å  CC gain = {cc_bent - cc_mov:+.6f}")
