#!/usr/bin/env ccp4-python
"""Combined d_opt selector for a fitreso_scan output.

For each fr-row (and best/, pre/, hkl00/): pull Rbent + CA-RMSD from
scan_fitreso.log, bondZ from check_geometry on bent.pdb, dipole_score
from an accompanying dipole_trajectory_<ref>.log.  Combine into

    S = R_bent + α·RMSD + β·max(0, dipole) + γ·max(0, bondZ − 1)

with α=0.1/Å, β=0.5, γ=0.05.  Report per-row scores and the argmin.
"""
import os, re, sys
sys.path.insert(0, '/home/jamesh/projects/map_bender/claude')
from bendfinder import check_geometry

REF = sys.argv[1] if len(sys.argv) > 1 else '1-6-100kGy'
BASE = f'/home/jamesh/projects/map_bender/arwen/JJD95_dose_series_mtzs/fitreso_scans/{REF}'
DIP  = f'/home/jamesh/projects/map_bender/arwen/JJD95_dose_series_mtzs/dipole_trajectory_{REF}.log'

# Weights
ALPHA = 0.1
BETA  = 0.5
GAMMA = 0.05

# Parse scan_fitreso.log for label → (rmsd, rbent)
scan_rows = {}
with open(f'{BASE}/scan_fitreso.log') as fh:
    for line in fh:
        if line.startswith('#') or line.startswith('-') or not line.strip():
            continue
        toks = line.split()
        if len(toks) < 5: continue
        label = toks[0]
        try:
            rmsd = float(toks[1])
        except ValueError:
            continue
        # Rbent is the token that ends in '%'
        rbent = None
        for t in toks[2:8]:
            if t.endswith('%'):
                try:
                    rbent = float(t.rstrip('%')) / 100.0
                    break
                except ValueError:
                    pass
        if rbent is None: continue
        scan_rows[label] = dict(rmsd=rmsd, rbent=rbent)

# Parse dipole trajectory
dip_rows = {}
with open(DIP) as fh:
    for line in fh:
        if line.startswith('#') or line.startswith('point') or line.startswith('---'):
            continue
        toks = line.split()
        if len(toks) < 3: continue
        try:
            dip_rows[toks[0]] = float(toks[2])
        except ValueError:
            continue

# Rows to include in the ranking — exclude pre/hkl0x since those don't
# have bondZ (no bent.pdb) and are baselines anyway.
FR_ROWS = ['fr20','fr15','fr12','fr10','fr8','fr7','fr6','fr5','best']

print(f"# combined d_opt score for {REF}")
print(f"# S = Rbent + {ALPHA}·RMSD(Å) + {BETA}·max(0,dipole) + "
      f"{GAMMA}·max(0,bondZ−1)")
print()
print(f"{'row':<6} {'RMSD':>6} {'Rbent':>6} {'dipole':>7} {'bondZ':>6}  "
      f"{'S':>6}  {'terms (bent  rmsd  dip  bond)':<30}")
print('-' * 82)

results = []
for lbl in FR_ROWS:
    r = scan_rows.get(lbl)
    if r is None: continue
    dip = dip_rows.get(lbl, 0.0)
    pdb = f'{BASE}/{lbl}/bent.pdb'
    if not os.path.exists(pdb):
        bondZ = None
    else:
        try:
            g = check_geometry(pdb)
            bondZ = g['bond_rmsz'] if g else None
        except Exception:
            bondZ = None

    if bondZ is None:
        bondZ = 1.0     # treat missing as neutral (no penalty)

    t_bent = r['rbent']
    t_rmsd = ALPHA * r['rmsd']
    t_dip  = BETA  * max(0.0, dip)
    t_bond = GAMMA * max(0.0, bondZ - 1.0)
    S = t_bent + t_rmsd + t_dip + t_bond
    results.append((lbl, r['rmsd'], r['rbent'], dip, bondZ, S,
                     (t_bent, t_rmsd, t_dip, t_bond)))

    print(f"{lbl:<6} {r['rmsd']:>6.3f} {r['rbent']*100:>5.1f}% "
          f"{dip:>+7.3f} {bondZ:>6.2f}  "
          f"{S:>6.3f}  "
          f"({t_bent:.3f} {t_rmsd:.3f} {t_dip:.3f} {t_bond:.3f})")

# Pick winner
best = min(results, key=lambda x: x[5])
print()
print(f"BEST by combined score: {best[0]}  S={best[5]:.3f}")
print(f"  RMSD={best[1]:.3f}Å  Rbent={best[2]*100:.1f}%  "
      f"dipole={best[3]:+.3f}  bondZ={best[4]:.2f}")
