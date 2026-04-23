# map_bender / bendfinder

Finds a smooth, periodic coordinate-shift field that maps the coordinates of one crystal onto a reference crystal of the same protein. The shift field is expressed as a truncated Fourier series using real Miller-index sine/cosine functions — physically motivated basis functions that respect the periodicity of the crystal lattice.

Primary author: James Holton

## What it does

Given two PDB files of the same protein (same or different crystal forms), `bendfinder.py` computes a smooth 3D vector field **Δr(x,y,z)** such that applying that field to the coordinates of `bendme.pdb` minimises the all-atom RMSD to `reference.pdb`. Optionally, any CCP4 map in the frame of `bendme.pdb` can be spline-interpolated into the reference frame.

The shift field is parameterised as:

```
Δd(x,y,z) = Σ_{hkl}  A_{hkl} · sin(2π(hx+ky+lz))  +  B_{hkl} · cos(2π(hx+ky+lz))
```

where the sum runs over (h,k,l) triplets sorted by resolution (low resolution first), so the most physically meaningful large-scale deformations are captured first. Symmetry constraints are applied: in space group P4₃2₁2 (8 proper rotations), 247 canonical parameters control 1401 Friedel-unique coefficients.

## Quick start

```python
from bendfinder import bend_fit_progressive, bend_apply_pdb

result = bend_fit_progressive('bendme.pdb', 'reference.pdb')
bend_apply_pdb('bendme.pdb', 'reference.pdb', result, outpath='bent.pdb')
print(f"CA RMSD: {result.rmsd:.3f} Å")
```

## Requirements

- Python 3.8+ with numpy, scipy, gemmi
- CCP4 Python environment (`ccp4-python`) recommended for map I/O

## Options (`bendfinder.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fitreso_start` | 20.0 | Starting resolution (Å) — coarsest HKLs first |
| `fitreso_end` | 7.0 | Stopping resolution (Å) — limited by overdetermination ratio |
| `batch_hkls` | 100 | HKLs admitted per progressive iteration |
| `od_margin` | 1.5 | Stop when (N_atoms × 3) / (N_canon × 6) drops below this |
| `outlier_sigma` | 2.5 | Reject CA pairs > this many robust σ from median shift |
| `b_sigma` | 3.0 | Reject CA pairs with B-factor > median + b_sigma × σ_MAD |
| `drop_snr` | 0.0 | Drop HKLs with \|A,B\| / σ < drop_snr (0 = keep all) |
| `use_symm` | True | Apply space-group proper-rotation constraints |
| `dimensions` | 'xyz' | Which coordinate axes to fit |
| `frac` | 1.0 | How far along the bend path to place the output (0–1) |
| `verbose` | True | Print per-iteration progress |

## How it works

1. Both PDB files are expanded from their crystallographic space group to P1 (via gemmi).
2. CA atom pairs are matched by residue and atom name. Outliers are rejected by shift magnitude and B-factor using robust statistics (median ± `outlier_sigma` × σ_MAD).
3. All (h,k,l) indices out to `fitreso_end` Å are enumerated, sorted by resolution (low-res first), and deduplicated for Friedel symmetry.
4. Proper-rotation symmetry of the space group is applied: each Friedel-unique HKL is assigned a canonical representative, reducing the free parameter count by the order of the point group (e.g. ×8 for P4₃2₁2).
5. A joint design matrix **X** (shape 3N_atoms × 6M_canon) is built from the symmetry-expanded sine/cosine basis. All three coordinate axes are fitted simultaneously.
6. **X @ params = shifts** is solved by SVD (scipy.linalg), giving the globally optimal (A, B) coefficients and their uncertainties from the covariance matrix.
7. Steps 2–6 are repeated in a coarse-to-fine loop (`fitreso_start` → `fitreso_end`), admitting `batch_hkls` new HKLs per iteration. Atoms whose residuals remain large after each fit are down-weighted, so the model is not distorted by conformational outliers.
8. The loop stops when the overdetermination ratio (atoms / canonical parameters) drops below `od_margin` — naturally limiting resolution to where the data can support the model.
9. The fitted shift field is evaluated at all atom positions; bent PDB and optional map are written.

## Benchmarks

All runs use default parameters (`fitreso_end=7.0 Å`, `batch_hkls=100`, `outlier_sigma=2.5`, `b_sigma=3.0`, `drop_snr=0`).

### Python version (`bendfinder.py`)

| System | Crystal forms | Space group | CA pairs | CA RMSD | Time |
|--------|--------------|-------------|----------|---------|------|
| Lysozyme | 3aw6 → 3aw7 | P4₃2₁2 | 1008 | **0.033 Å** | 118 s |
| DHFR | 1rx1 → 1rx2 | P2₁2₁2₁ | 636 | **0.071 Å** | 191 s |
| Myoglobin | 1mbo → 1a6m | P2₁ | 294 | **0.060 Å** | 106 s |
| Lysozyme raddam | 5kxk → 5kxl | P4₃2₁2 | 976 | **0.082 Å** | ~550 s |
| Lysozyme raddam | 5kxk → 5kxm | P4₃2₁2 | 984 | **0.046 Å** | ~575 s |
| Lysozyme raddam | 5kxk → 5kxn | P4₃2₁2 | 992 | **0.055 Å** | ~580 s |

### vs prototype (tcsh + gnuplot)

| System | Prototype RMSD | Prototype time | Python RMSD | Python time | Speedup |
|--------|---------------|----------------|-------------|-------------|---------|
| Lysozyme 3aw6/3aw7 | 0.209 Å | 2938 s | 0.033 Å | 118 s | **25×** |
| Myoglobin 1mbo/1a6m | 0.229 Å | 237 s | 0.060 Å | 106 s | **2×** |
| Raddam 5kxk→5kxl | 0.177 Å | 643 s | 0.082 Å | ~550 s | **same speed, 2× better** |
| Raddam 5kxk→5kxm | 0.165 Å | 506 s | 0.046 Å | ~575 s | **3.6× better** |
| Raddam 5kxk→5kxn | 0.220 Å | 558 s | 0.055 Å | ~580 s | **4× better** |

The lysozyme comparison is against the gold-standard prototype run (order 5, 91 HKLs, 71.9% vs 84.2% relative humidity causing ~2.5% cell contraction). The Python version achieves 6× better RMSD in 1/25th the time, primarily because:
- Linear (A, B) parameterisation allows all HKLs to be fitted simultaneously via SVD
- Space-group symmetry constraints reduce free parameters ~8× for P4₃2₁2
- No incremental gnuplot fitting loop required

The myoglobin prototype used nhkls=100; the Python version naturally fits more HKLs (401) before the overdetermination ratio stops it, contributing to the lower RMSD.

### Prototype convergence (lysozyme 3aw6/3aw7, for reference)

| HKLs fitted | RMSD(CA) | Time  | Old approach equivalent        |
|-------------|----------|-------|-------------------------------|
| 5           | 0.327 Å  | 9 s   | Order 2 (19 HKLs)             |
| 20          | 0.245 Å  | 580 s | Order 3 (37 HKLs)             |
| 26          | 0.215 Å  | 1070 s| Order 4 (61 HKLs)             |
| 30          | 0.211 Å  | 1503 s| Order 5 (91 HKLs, 2938 s)    |

## Output

`bend_fit_progressive` returns a result object with:

| Attribute | Description |
|-----------|-------------|
| `result.rmsd` | CA RMSD after bending (Å) |
| `result.hkls` | Array of active (h,k,l) indices, shape (N,3) |
| `result.AB` | Fourier coefficients, shape (3, N, 2) — [dim, hkl, sin/cos] |
| `result.active` | Boolean mask of active HKLs |
| `result.snr` | Signal-to-noise ratio per HKL |
| `result.cell1`, `result.cell2` | Unit cell parameters for both crystals |
| `result.dimensions` | Which coordinate axes were fitted |

Key functions:

| Function | Description |
|----------|-------------|
| `bend_fit_progressive(pdb1, pdb2, ...)` | Fit the shift field; returns result object |
| `bend_apply_pdb(pdb1, pdb2, result, outpath=...)` | Write bent PDB |
| `save_fitparams(path.mtz, ...)` | Save fitted coefficients to MTZ |
| `load_fitparams(path.mtz)` | Load previously fitted coefficients |
| `eval_shift_field(frac_pts, hkls, AB)` | Evaluate shift at arbitrary fractional coordinates |
| `interpolate_map(map_data, header, frac_pts)` | Trilinear interpolation into a CCP4 map |
| `read_ccp4(path)` / `write_ccp4(path, data, header)` | CCP4 map I/O |

## Files in this repository

| File | Description |
|------|-------------|
| `bendfinder.py` | Main Python module |
| `bendfinder.com` | Prototype tcsh script (historical reference) |
| `origins.com` | Helper: test symmetry origin choices |
| `examples/3aw6_3aw7/` | Lysozyme canonical example data |
| `LICENSE` | License |
