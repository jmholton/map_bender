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
- CCP4 programs (scaleit, cad) required only for Riso computation via `diff.com`; map→MTZ conversion uses gemmi directly

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
| `iter_callback` | None | Called as `f(iter_i, nhkls, n_canon, rmsd, hkls, AB_xyz, active, snr)` after each progressive iteration |
| `max_canon` | None | Stop after first iteration where n_canon ≥ this value; useful for stepping through canonical HKL checkpoints |

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

| System | Crystal forms | Space group | CA pairs | CA RMSD | Riso (fr5) | Time |
|--------|--------------|-------------|----------|---------|-----------|------|
| Lysozyme | 3aw6 → 3aw7 | P4₃2₁2 | 1008 | **0.034 Å** | 33.2% | 335 s |
| DHFR | 1rx1 → 1rx2 | P2₁2₁2₁ | 636 | **0.071 Å** | 41.9% | 287 s |
| Myoglobin | 1mbo → 1a6m | P2₁ | 294 | **0.063 Å** | — | 150 s |
| Lysozyme raddam | 5kxk → 5kxl | P4₃2₁2 | 976 | **0.082 Å** | 20.6% | 112 s |
| Lysozyme raddam | 5kxk → 5kxm | P4₃2₁2 | 984 | **0.046 Å** | — | — |
| Lysozyme raddam | 5kxk → 5kxn | P4₃2₁2 | ~992 | — | — | — |

Times are for the fr5 fitreso scan run (fitreso_start=20, fitreso_end=5 Å). Raddam times are faster because the Fc maps have smaller structural differences.

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
| `interpolate_map(map_data, header, frac_pts)` | Tricubic spline interpolation into a CCP4 map; 5-voxel periodic padding prevents b-spline overshoot at cell boundaries |
| `read_ccp4(path)` / `write_ccp4(path, data, header)` | CCP4 map I/O |
| `read_ccp4_fullcell(path)` | Read CCP4 map and expand to full P1 unit cell via gemmi symmetry |

## Files in this repository

| File | Description |
|------|-------------|
| `bendfinder.py` | Main Python module |
| `bendfinder.com` | Prototype tcsh script (historical reference) |
| `origins.com` | Helper: test symmetry origin choices |
| `examples/3aw6_3aw7/` | Lysozyme canonical example data |
| `LICENSE` | License |

## Fitreso scan workflow

The scan scripts (`lyso/lyso_fitreso_scan.py`, `dhfr/dhfr_fitreso_scan.py`, `myoglobin/myoglobin_fitreso_scan.py`, `raddam/raddam_fitreso_scan.py`) evaluate bend quality as a function of resolution. Each script runs three sections:

1. **hkl00** — no bending; resample moving map onto reference grid with zero shift. Gives the unregistered baseline Rfac.
2. **hkl01..hkl10** — add canonical HKLs one at a time using a single `bend_fit_progressive` call with `batch_hkls=1, max_canon=11` and an `iter_callback` that snapshots the bent map at each n_non_DC checkpoint.
3. **fr20..fr5** — separate `bend_fit_progressive` calls with `fitreso_end` stepped from 20 to 5 Å.

For each scan point the script writes `bent.map`, `diff_norm.map` (z-scored difference), and `fitparams.mtz`, then computes Riso via `diff.com` and reports the strongest positive and negative difference-map peaks with their nearest atoms.

Map→MTZ conversion uses `gemmi.transform_map_to_f_phi` (no sfall required).
