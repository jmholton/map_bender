#!/programs/ccp4-8.0/bin/ccp4-python
"""test_symm_all_sgs.py — symmetry constraint test across all 65 Sohncke SGs.

For each protein-compatible (Sohncke) space group:
  1. Place N_ASU_ATOMS atoms at random general positions in the ASU.
  2. Apply small random displacements δ_i to each ASU atom.
  3. Expand BOTH states to P1 via all SG operators:
       ref copy k of atom i  → R_k · x_i  + t_k
       moved copy k of atom i → R_k · (x_i + δ_i) + t_k
     ∴ displacement of copy k of atom i = R_k · δ_i   (SG-consistent!)
     This is the correct crystallographic scenario: both structures in the SG.
  4. Fit PSDVF two ways:
       unconstrained — build_design_matrix on ALL N_ops × N_ASU atoms
       constrained   — build_design_matrix_symm on ASU atoms only
  5. Evaluate both fields on a grid and measure the SG symmetry violation:
       max_{k,x} |Δr(R_k·x+t_k) − R_k·Δr(x)|  in Å
  6. Report: SG, n_ops, n_atoms_P1, n_canon, viol(uncon), viol(con).

Key question: does the unconstrained fit (all P1 atoms, no SG constraint) naturally
produce an SG-symmetric field when the input data is SG-symmetric?
This depends on whether the Friedel-unique HKL set is closed under the point group.
The constrained fit enforces SG symmetry by construction → violation should be ~0.
"""

import sys, os
import numpy as np
import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from bendfinder import (
    generate_hkls, get_proper_symops, assign_canonical,
    build_design_matrix, fit_lstsq,
    build_design_matrix_symm, fit_lstsq_symm, expand_ab_canon,
    eval_shift_field, _ortho_matrix,
)

# ── parameters ────────────────────────────────────────────────────────────────
N_ASU_ATOMS = 8       # atoms in ASU (general positions)
SHIFT_AMP   = 0.008   # fractional displacement amplitude (~0.3 Å for 40 Å cell)
GRID_N      = 15      # violation check grid: 15³ = 3375 points
FITRESO     = 20.0    # Å — use ALL HKLs at this resolution (no truncation)
                      # The resolution sphere is point-group-invariant, so the
                      # full set is orbit-complete; truncation would break closure.


# ── cell parameters for each crystal system ───────────────────────────────────
def standard_cell(sg_number):
    if sg_number <= 2:
        return (30., 35., 40., 70., 80., 75.)          # Triclinic
    elif sg_number <= 15:
        return (30., 40., 35., 90., 100., 90.)          # Monoclinic
    elif sg_number <= 74:
        return (30., 40., 50., 90., 90., 90.)           # Orthorhombic
    elif sg_number <= 142:
        return (40., 40., 55., 90., 90., 90.)           # Tetragonal
    elif sg_number <= 194:
        return (40., 40., 65., 90., 90., 120.)          # Trigonal / Hexagonal
    else:
        return (55., 55., 55., 90., 90., 90.)           # Cubic


def crystal_system(n):
    if n <= 2:   return 'Triclinic'
    if n <= 15:  return 'Monoclinic'
    if n <= 74:  return 'Orthorhombic'
    if n <= 142: return 'Tetragonal'
    if n <= 167: return 'Trigonal'
    if n <= 194: return 'Hexagonal'
    return 'Cubic'


# ── enumerate Sohncke (protein-compatible) space groups ───────────────────────
def sohncke_sgs():
    result = []
    for i in range(1, 231):
        sg = gemmi.find_spacegroup_by_number(i)
        if all(round(np.linalg.det(np.array(op.rot, dtype=float) / 24.)) == 1
               for op in sg.operations()):
            result.append(sg)
    return result


# ── symmetry violation metric ─────────────────────────────────────────────────
def compute_violation(AB_xyz, hkls, proper_ops, M):
    """Max over all non-identity ops of |Δr(Rx+t) − R·Δr(x)| in Å on grid."""
    xi  = np.linspace(0, 1, GRID_N, endpoint=False)
    gx, gy, gz = np.meshgrid(xi, xi, xi, indexing='ij')
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)    # (N,3)

    dr       = eval_shift_field(pts, hkls, AB_xyz)                   # (N,3) frac
    max_viol = 0.
    for R, t in proper_ops[1:]:                                      # skip identity
        pts_sym  = (pts @ R.T + t) % 1.0
        dr_sym   = eval_shift_field(pts_sym, hkls, AB_xyz)           # Δr(Rx+t)
        dr_R     = dr @ R.T                                          # R·Δr(x)
        viol_ang = np.linalg.norm((dr_sym - dr_R) @ M.T, axis=1).max()
        max_viol = max(max_viol, viol_ang)
    return max_viol


# ── per-SG test ───────────────────────────────────────────────────────────────
def run_sg(sg):
    cell = standard_cell(sg.number)
    M    = _ortho_matrix(cell)

    proper_ops = get_proper_symops(sg.hm)         # list of (R, t) in fractional
    n_ops      = len(proper_ops)

    # Random ASU atom positions (general, interior of [0.1, 0.4]³ avoids specials)
    rng      = np.random.default_rng(sg.number)   # reproducible per SG
    asu_ref  = rng.uniform(0.1, 0.4, (N_ASU_ATOMS, 3))     # (N_ASU, 3) fractional
    delta_asu = rng.standard_normal((N_ASU_ATOMS, 3)) * SHIFT_AMP  # (N_ASU, 3) frac

    # Expand BOTH states to P1; displacement of op-k copy of atom i = R_k · δ_i
    ref_pts_list, shifts_list = [], []
    for R, t in proper_ops:
        for i in range(N_ASU_ATOMS):
            ref_pts_list.append((R @ asu_ref[i] + t) % 1.0)
            shifts_list.append(R @ delta_asu[i])
    ref_pts = np.array(ref_pts_list)   # (N_ops * N_ASU, 3)
    shifts  = np.array(shifts_list)    # (N_ops * N_ASU, 3) — SG-consistent

    # HKL pool — ALL HKLs at FITRESO; no truncation (orbit-complete set)
    hkls = generate_hkls(cell, FITRESO)
    n_h  = len(hkls)

    # ── unconstrained fit (all N_ops × N_ASU P1 atoms) ───────────────────────
    X_un       = build_design_matrix(ref_pts, hkls)          # (N_P1, 2*n_h)
    AB_un, _, _ = fit_lstsq(X_un, shifts, drop_snr=0.)       # (3, n_h, 2)

    # ── constrained fit (ASU atoms only, SG symmetry enforced) ───────────────
    canon_hkls, hkl_to_canon, hkl_all_ops = assign_canonical(hkls, proper_ops)
    n_canon = len(canon_hkls)

    X_con    = build_design_matrix_symm(
                   asu_ref, canon_hkls, proper_ops, hkl_to_canon, hkl_all_ops)
    b_con    = delta_asu.ravel()                              # (3*N_ASU,) interleaved
    params, active_c, _ = fit_lstsq_symm(X_con, b_con, drop_snr=0.)
    AB_con   = expand_ab_canon(params, active_c, canon_hkls, hkls,
                               hkl_to_canon, hkl_all_ops, proper_ops)  # (3,n_h,2)

    # ── symmetry violation ────────────────────────────────────────────────────
    viol_un  = compute_violation(AB_un,  hkls, proper_ops, M)
    viol_con = compute_violation(AB_con, hkls, proper_ops, M)

    return {
        'sg':       sg.hm,
        'number':   sg.number,
        'cs':       crystal_system(sg.number),
        'n_ops':    n_ops,
        'n_p1':     n_ops * N_ASU_ATOMS,
        'n_h':      n_h,
        'n_canon':  n_canon,
        'viol_un':  viol_un,
        'viol_con': viol_con,
    }


# ── main ──────────────────────────────────────────────────────────────────────
sgs = sohncke_sgs()
print(f"Testing {len(sgs)} protein-compatible (Sohncke) space groups\n"
      f"  ASU atoms:   {N_ASU_ATOMS} (random general positions)\n"
      f"  Shift amp:   {SHIFT_AMP} fractional ({SHIFT_AMP*40:.2f} Å for 40 Å cell)\n"
      f"  HKLs/SG:     all Friedel-unique at d ≥ {FITRESO} Å (orbit-complete)\n"
      f"  Grid:        {GRID_N}³ = {GRID_N**3} points\n")

hdr = (f"{'SG':12s}  {'#':4s}  {'system':12s}  {'ops':4s}  {'P1atm':6s}  "
       f"{'hkls':5s}  {'canon':5s}  {'redc':5s}  "
       f"{'viol(uncon) Å':14s}  {'viol(con) Å':12s}")
print(hdr)
print('-' * len(hdr))

results = []
for sg in sgs:
    r = run_sg(sg)
    results.append(r)
    redc = r['n_h'] / max(r['n_canon'], 1)
    ok   = 'OK' if r['viol_con'] < 1e-6 else 'FAIL'
    print(f"{r['sg']:12s}  {r['number']:4d}  {r['cs']:12s}  {r['n_ops']:4d}  "
          f"{r['n_p1']:6d}  {r['n_h']:5d}  {r['n_canon']:5d}  {redc:5.1f}x  "
          f"{r['viol_un']:14.6f}  {r['viol_con']:12.2e}  {ok}")

# ── summary ───────────────────────────────────────────────────────────────────
n_pass = sum(1 for r in results if r['viol_con'] < 1e-6)
n_nonz = sum(1 for r in results if r['viol_un']  > 1e-6)
max_con = max(r['viol_con'] for r in results)
mean_un = np.mean([r['viol_un'] for r in results if r['n_ops'] > 1])

print(f"\n{'='*70}")
print(f"Constrained fit:   {n_pass}/{len(results)} SGs have violation < 1e-6 Å  "
      f"(max = {max_con:.2e} Å)")
print(f"Unconstrained fit: {n_nonz}/{len(results)} SGs have violation > 1e-6 Å  "
      f"(mean for n_ops>1 = {mean_un:.4f} Å)")
print("Done.")
