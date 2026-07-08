#!/programs/ccp4-8.0/bin/ccp4-python
"""test_magdoff.py — Magdoff-type synthetic deformation tests on ribonuclease A (7rsa).

7rsa: ribonuclease A, P2₁, a=30.2 b=38.4 c=53.3 Å β=105.9°, 124 residues, 248 CA in P1.

Test 1 — Cell scaling (0.5% isomorphous cell change):
  Scale a, b, c in CRYST1 by 1.005; leave atom Cartesian coordinates unchanged.
  Simulates a pure unit-cell expansion with no intramolecular rearrangement.
  The fractional coordinates in the scaled cell are x/1.005 ≈ x·(1−0.005), so the
  fractional shift Δx ≈ −0.005·x_frac and the Cartesian shift Δr ≈ −0.005·r_cart.

Test 2 — Rigid-body rotation (0.5° about a fixed random axis):
  Rotate all atoms by 0.5° in Cartesian space.  Expected shift: Δr ≈ ω × r_cart.
  Not SG-consistent for a general axis (the 2₁ screw mixes the rotation).

Each test is run twice — use_symm=True (P2₁-constrained, 2× fewer parameters) and
use_symm=False (unconstrained, can fit non-SG-consistent components) — to show the
tradeoff.
"""
import sys, os, math, numpy as np, gemmi
sys.path.insert(0, '/home/jamesh/projects/map_bender/claude')
from bendfinder import (_ortho_matrix, bend_fit_progressive, eval_shift_field,
                        bend_apply_pdb)

PDB_REF = '/home/jamesh/projects/map_bender/claude/magdoff/7rsa.pdb'
WORKDIR = '/home/jamesh/projects/map_bender/claude/magdoff'

RNG = np.random.default_rng(42)

DMIN_RISO = 1.5
RATE_RISO = 1.5


# ── helpers ───────────────────────────────────────────────────────────────────

def rotation_matrix(axis, angle_deg):
    """Rodrigues rotation: axis is a 3-vector (need not be unit), angle in degrees."""
    n = np.asarray(axis, float); n /= np.linalg.norm(n)
    th = math.radians(angle_deg)
    c, s = math.cos(th), math.sin(th)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3)*c + s*K + (1-c)*np.outer(n, n)


def modify_pdb(src, dst, fn):
    """Copy PDB, applying fn(xyz) → xyz to each ATOM/HETATM record."""
    with open(src) as f, open(dst, 'w') as g:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                xyz = np.array([float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])])
                xp, yp, zp = fn(xyz)
                line = (line[:30] + f"{xp:8.3f}{yp:8.3f}{zp:8.3f}" + line[54:])
            g.write(line)


def scale_cell_pdb(src, dst, scale):
    """Copy PDB, scaling a, b, c in CRYST1 by scale; strip SCALE/ORIGX cards;
    atom coordinates unchanged.  Stripping forces gemmi to derive the fractional
    coordinate matrix from the (modified) CRYST1 rather than the stale cards."""
    with open(src) as f, open(dst, 'w') as g:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15]) * scale
                b = float(line[15:24]) * scale
                c = float(line[24:33]) * scale
                line = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{line[33:]}"
            elif line.startswith(('SCALE', 'ORIGX')):
                continue
            g.write(line)


def pdb_to_f(pdb_path, dmin=DMIN_RISO, rate=RATE_RISO):
    """Fcalc structure factors from PDB via DensityCalculatorX → dict {hkl: |F|}."""
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()
    dc = gemmi.DensityCalculatorX()
    dc.d_min = dmin
    dc.rate  = rate
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    sf = gemmi.transform_map_to_f_phi(dc.grid, half_l=True)
    data = sf.prepare_asu_data(dmin=dmin)
    hkl = data.miller_array; fabs = np.abs(data.value_array)
    return {tuple(hkl[i].tolist()): fabs[i] for i in range(len(fabs))}


def riso(f_ref, f_other):
    """Isomorphous R-factor (%) between two {hkl: |F|} dicts, common HKLs only."""
    common = set(f_ref.keys()) & set(f_other.keys())
    num = sum(abs(f_ref[h] - f_other[h]) for h in common)
    den = sum(f_ref[h] for h in common)
    return 100.*num/den, len(common)


def one_fit(pdb_mod, pdb_ref, use_symm):
    """Return (rmsd_before_Å, result) for one fit."""
    result = bend_fit_progressive(
        pdb_mod, pdb_ref,
        fitreso_start=20.0, fitreso_end=7.0,
        batch_hkls=100, od_margin=1.5,
        outlier_sigma=2.5, b_sigma=3.0,
        use_symm=use_symm, verbose=True,
    )
    fitme = result.fitme
    ca    = result.ca_mask
    M     = _ortho_matrix(result.cell1)
    cart  = fitme[ca, 3:6] @ M.T
    rmsd_before = float(np.sqrt(np.mean(np.sum(cart**2, axis=1))))
    return rmsd_before, result


def run_test(label, pdb_mod, pdb_ref):
    print(f"\n{'═'*66}")
    print(f"  {label}")
    print(f"{'═'*66}")

    # Constrained (P2₁-symmetric) fit
    print("\n  --- use_symm=True (P2₁-constrained) ---")
    rmsd_before, r_symm = one_fit(pdb_mod, pdb_ref, use_symm=True)
    recovery_symm = (1.0 - r_symm.rmsd / rmsd_before) * 100.0

    # Unconstrained fit
    print("\n  --- use_symm=False (unconstrained) ---")
    _, r_unsymm = one_fit(pdb_mod, pdb_ref, use_symm=False)
    recovery_unsymm = (1.0 - r_unsymm.rmsd / rmsd_before) * 100.0

    # Summary
    print(f"\n  ── Summary ──────────────────────────────────────────────")
    print(f"  Imposed deformation RMSD(CA)    : {rmsd_before:.4f} Å")
    print(f"  Constrained   RMSD(CA) after fit: {r_symm.rmsd:.4f} Å  "
          f"({recovery_symm:.1f}% recovery)  "
          f"{r_symm.active.sum()} active Friedel HKLs [P2₁-constrained]")
    print(f"  Unconstrained RMSD(CA) after fit: {r_unsymm.rmsd:.4f} Å  "
          f"({recovery_unsymm:.1f}% recovery)  "
          f"{r_unsymm.active.sum()} active Friedel HKLs [unconstrained]")

    # Show first 3 CA shifts (from the constrained fit)
    fitme = r_symm.fitme;  ca = r_symm.ca_mask
    M     = _ortho_matrix(r_symm.cell1)
    true_cart  = fitme[ca, 3:6] @ M.T
    fit_frac   = eval_shift_field(fitme[ca, :3], r_symm.hkls, r_symm.AB)
    fit_cart   = fit_frac @ M.T
    err_mag    = np.linalg.norm(true_cart - fit_cart, axis=1)

    print(f"\n  First 3 CA atoms (constrained fit):")
    for i in range(min(3, int(ca.sum()))):
        ts = true_cart[i]
        rs = err_mag[i]
        print(f"    CA{i+1}: shift=[{ts[0]:+.3f},{ts[1]:+.3f},{ts[2]:+.3f}] Å  "
              f"|shift|={np.linalg.norm(ts):.4f} Å  |resid|={rs:.4f} Å")

    # ── Riso ──────────────────────────────────────────────────────────────────
    print(f"\n  ── Riso at {DMIN_RISO} Å ────────────────────────────────────────────")
    f_ref = pdb_to_f(pdb_ref)
    f_mod = pdb_to_f(pdb_mod)
    r_before, n_hkl = riso(f_ref, f_mod)
    print(f"  Before bending       : Riso={r_before:.2f}%  ({n_hkl} HKLs)")

    base = os.path.basename(pdb_mod).replace('.pdb', '')
    bent_symm_path = os.path.join(WORKDIR, f"{base}_bent_symm.pdb")
    bend_apply_pdb(pdb_mod, pdb_ref, r_symm, outpath=bent_symm_path)
    r_rs, _ = riso(f_ref, pdb_to_f(bent_symm_path))
    print(f"  Constrained   (symm) : Riso={r_rs:.2f}%  ({n_hkl} HKLs)")

    bent_un_path = os.path.join(WORKDIR, f"{base}_bent_un.pdb")
    bend_apply_pdb(pdb_mod, pdb_ref, r_unsymm, outpath=bent_un_path)
    r_ru, _ = riso(f_ref, pdb_to_f(bent_un_path))
    print(f"  Unconstrained        : Riso={r_ru:.2f}%  ({n_hkl} HKLs)")


# ── Test 1: uniform cell scaling 0.5% (CRYST1 only, atoms unchanged) ─────────
SCALE      = 1.005
pdb_scaled = os.path.join(WORKDIR, '7rsa_scaled.pdb')
scale_cell_pdb(PDB_REF, pdb_scaled, SCALE)
run_test(f"Test 1 — Cell scaling ×{SCALE}  (CRYST1 only, atoms unchanged)", pdb_scaled, PDB_REF)


# ── Test 2: rigid-body rotation 0.5° about a random axis ─────────────────────
ROT_AXIS = RNG.standard_normal(3); ROT_AXIS /= np.linalg.norm(ROT_AXIS)
ROT_DEG  = 0.5
R_ROT    = rotation_matrix(ROT_AXIS, ROT_DEG)

print(f"\nRotation axis (Cartesian): [{ROT_AXIS[0]:.4f}, {ROT_AXIS[1]:.4f}, {ROT_AXIS[2]:.4f}]")
print(f"Rotation angle           : {ROT_DEG}°")

pdb_rotated = os.path.join(WORKDIR, '7rsa_rotated.pdb')
modify_pdb(PDB_REF, pdb_rotated, lambda r: R_ROT @ r)
run_test(
    f"Test 2 — Rigid-body rotation {ROT_DEG}° "
    f"[{ROT_AXIS[0]:.3f},{ROT_AXIS[1]:.3f},{ROT_AXIS[2]:.3f}]",
    pdb_rotated, PDB_REF,
)

print("\nDone.")
