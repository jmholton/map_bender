# map_bender / bendfinder

Finds a smooth, periodic coordinate-shift field that maps the coordinates of one crystal onto a reference crystal of the same protein. The shift field is expressed as a truncated Fourier series using real Miller-index sine/cosine functions — physically motivated basis functions that respect the periodicity of the crystal lattice.

Primary author: James Holton

## What it does

Given two PDB files of the same protein — same crystal form at different conditions (humidity, temperature, ligand), or genuinely different crystal forms — `bendfinder.py` computes a smooth 3D vector field **Δr(x,y,z)** such that applying that field to the coordinates of `bendme.pdb` minimises the all-atom RMSD to `reference.pdb`. Optionally, any CCP4 map in the frame of `bendme.pdb` can be spline-interpolated into the reference frame.

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
- No CCP4 programs required; Riso computation and map→MTZ conversion are pure Python via gemmi

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

| System | Datasets | Space group | CA pairs | fr5 RMSD | fr5 Rbent | best Rbent | d_opt |
|--------|----------|-------------|----------|----------|-----------|------------|-------|
| Lysozyme (humidity-driven cell change) | 3aw6 → 3aw7 | P4₃2₁2 | 1008 | **0.033 Å** | 33.2% | 30.3% | 9.7 Å |
| DHFR | 1rx2 → 1rx1 | P2₁2₁2₁ | 592 | **0.070 Å** | 41.5% | 37.7% | 12.6 Å |
| Myoglobin | 1mbo → 1a6m | P2₁ | 294 | **0.063 Å** | 49.8% | 49.8% | 5.5 Å |
| Lysozyme raddam | 5kxk → 5kxl | P4₃2₁2 | 976 | **0.087 Å** | 21.9% | 11.5% | 20 Å (clamped) |
| Lysozyme raddam | 5kxk → 5kxm | P4₃2₁2 | 984 | **0.047 Å** | 19.2% |  9.9% | 20 Å (clamped) |
| Lysozyme raddam | 5kxk → 5kxn | P4₃2₁2 | 992 | **0.048 Å** | 24.1% | 17.8% | 20 Å (clamped) |
| Insulin T→R | 4fg3 → 4e7u | H3 | 801 | **0.510 Å** | 60.3% | 60.5% | 5.8 Å |
| Porin | 3poq → 3pou | H 3 2 | 340 | **0.366 Å** | 57.7% | 57.8% | 9.5 Å |

`fr5 RMSD` and `fr5 Rbent` are the values at the finest resolution
admitted (5 Å); `best Rbent` and `d_opt` come from the parabola-fit
re-run (see [Best d_opt parabola fit](#best-d_opt-parabola-fit) below).
Insulin requires `fill_fcalc=True` because both deposited MTZs are
below 99% SG-ASU complete; porin runs end-to-end with the
in-bendfinder altindex resolution (the obverse/reverse 2-fold for the
H 3 2 pair) — see [Obverse/reverse and altalign.py](#obversereverse-and-altalignpy)
for the dual-solution writer when you need an R 3 2:R MTZ for
downstream refinement.

The radiation-damage systems gain 8–13 Rbent points from the parabola
because the damage signal is purely low-frequency — `d_opt` clamps to
the coarsest scan point (20 Å) where the smooth shift field captures
the bulk swelling and the high-frequency HKLs only add noise.
Conversely, insulin's T→R conformational change exceeds the smooth-
PSDVF model (LEU B6 shifts ~8 Å), so no fitreso choice helps and the
parabola vertex collapses near fr5.

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

## Synthetic validation (Magdoff tests)

`magdoff/test_magdoff.py` imposes two controlled deformations on ribonuclease A (7rsa, P2₁, 248 CA in P1) and measures how well bendfinder recovers them. Riso is computed at 1.5 Å (full crystallographic resolution) using gemmi Fcalc; PSDVF coefficients are fitted to 7 Å spatial resolution — these are distinct.

**Test 1 — Isomorphous cell change (Magdoff):** Scale a, b, c in CRYST1 by 1.005; leave atom Cartesian coordinates unchanged. This is the classic Magdoff scenario: same protein, slightly different unit cell. Riso 15.0% → 2.7%; RMSD(CA) 0.197 → 0.018 Å (91% recovery).

**Test 2 — Rigid-body rotation 0.5°:** Rotate all atom Cartesian coordinates by 0.5° about a random axis. Riso 24.0% → 4.3%; RMSD(CA) 0.335 → 0.028 Å (92% recovery).

The residual Riso after bending (~2.7% and ~4.3%) is set by the 7 Å PSDVF bandwidth: structure factor content at 1.5–7 Å that the shift field cannot represent.

> **Implementation note:** When modifying only the CRYST1 record in a PDB file (cell dimensions), the SCALE and ORIGX records must also be stripped. If left, gemmi reads the fractional coordinate matrix from those stale cards rather than deriving it from the modified CRYST1, so the fractional coordinates are unchanged and no deformation is seen.

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
| `compute_riso(ref_mtz, ref_col, test_mtz, test_col)` | Wilson isotropic B scaling + Riso; returns (riso, kF, B_iso); pure Python, no CCP4 |
| `fitreso_scan(mov_pdb, ref_pdb, mov_mtz, ref_mtz, scan_dir, ...)` | Run full hkl00..hkl10 + fr20..fr5 + best scan; write bent maps, diff maps, bent.mtz per point. Accepts MTZ or CCP4-map inputs; `subtract` and `fill_fcalc` arguments control diff sign and missing-HKL filling. |

## Space-group generality

The shift field must respect the crystallographic symmetry of the space group:

```
Δr(R·x + t) = R · Δr(x)   for all {R, t} in the space group
```

`use_symm=True` (default) enforces this by fitting only canonical HKL representatives and expanding to the full Friedel-unique set via the point-group expansion formula. The free-parameter count is reduced by the order of the proper point group: ×2 for monoclinic, ×4 for orthorhombic, ×8 for P4₃2₁2, ×12 for cubic, etc.

**All 65 Sohncke (protein-compatible) space groups are tested**, including all screw-axis and non-symmorphic settings, all Bravais lattice types (P, C, I, F, R), and all trigonal and hexagonal groups. The symmetry constraint is satisfied to machine precision (<10⁻¹² Å) in every case. The test (`test_symm_all_sgs.py`) places 8 atoms in a general ASU position, displaces each independently, expands both states to P1 via the SG operators, fits with and without the symmetry constraint, and checks the residual violation on a 15³ grid.

**Implementation note.** In fractional coordinates, proper rotation matrices R are not in general orthogonal — R^T ≠ R^{−1} for trigonal and hexagonal crystal systems. The vectorial part of the symmetry expansion requires R^{−1}, not R^T. For all other crystal systems the two are equal, so the distinction only matters for trigonal/hexagonal.

## Files in this repository

| File | Description |
|------|-------------|
| `bendfinder.py` | Main Python module |
| `altalign.py` | Standalone LSQ altindex+origin search; dual-solution writer for R-lattice non-normalizing ops (H32+SYMM and R32:R outputs) |
| `bendfinder.com` | Prototype tcsh script (historical reference) |
| `origins.com` | Helper: test symmetry origin choices |
| `examples/3aw6_3aw7/` | Lysozyme canonical example data |
| `LICENSE` | License |

## Fitreso scan workflow

`fitreso_scan()` in `bendfinder.py` evaluates bend quality as a function of resolution. Call it directly:

```python
import sys; sys.path.insert(0, '/path/to/map_bender')
from bendfinder import fitreso_scan
fitreso_scan(mov_pdb='5kxk.pdb', ref_pdb='5kxl.pdb',
             mov_map='5kxk_fc.map', ref_map='5kxl_fc.map',
             scan_dir='fitreso_scan_5kxl')
```

Three sections:

1. **hkl00** — no bending; resample moving map onto reference grid with zero shift. Gives the unregistered baseline Rfac.
2. **hkl01..hkl10** — add canonical HKLs one at a time using a single `bend_fit_progressive` call with `batch_hkls=1, max_canon=11` and an `iter_callback` that snapshots the bent map at each n_non_DC checkpoint.
3. **fr20..fr5** — separate `bend_fit_progressive` calls with `fitreso_end` stepped from 20 to 5 Å.

For each scan point the function writes:

- `bent.mtz` — `FDM`/`PHIDM` (bent 2Fo-Fc as structure factors) plus `DELFWT`/`PHDELWT` (the difference map). **Load in Coot** via *File → Open MTZ*; `FDM`/`PHIDM` are recognised by default. Display `DELFWT`/`PHDELWT` as a difference map contoured at ±3σ.
- `diff_norm.map` — z-scored real-space difference map. Sign controlled by the `subtract` parameter: `subtract='ref'` (default) → diff = bent − ref, positive peaks = density present in bent but absent in ref; `subtract='bent'` flips it (use for radiation-damage runs where you want positive = features appearing with dose).
- `bent.map` — bent moving-crystal map resampled on the reference grid (CCP4 format)
- `bent.pdb`, `ref.pdb`, `unbent.pdb` — coordinates for atom context in the viewer
- `PSDVF.mtz` (fr\* + best points) — fitted shift-field (h,k,l) coefficients. Amplitudes are in Å (= √(A²+B²) × cell_axis) — `load_fitparams` reads the unit tag from the MTZ history and converts back to fractional on load.

Riso is computed by `compute_riso()`: Wilson isotropic B scaling of `F_test` to match `F_ref`, then Riso = Σ|F_ref − F_scaled| / Σ|F_ref|. No CCP4 programs needed. Because both maps include model bias and water molecules, reported Riso values should not be compared directly to crystallographic Riso from experimental data. The dominant peaks at high resolution (fr7 and finer) are typically displaced water molecules.

**Raddam sign convention**: 5kxk (undamaged) is the moving model; 5kxl/5kxm/5kxn (increasing dose) are references. Run with `subtract='bent'` so positive diff peaks = features appearing with radiation dose; negative peaks = features disappearing.

### Best d_opt parabola fit

Across every example system the Rbent-vs-fitreso curve is a clear U: it drops as the smooth PSDVF absorbs structural detail with finer resolution, then climbs again at the highest resolutions where the high-frequency HKLs add noise outside the shift field's natural bandwidth. `fitreso_scan` exploits this by fitting a parabola in **x = 1/d²** (the natural axis for R-factor-vs-resolution behaviour) to the 3–5 fr-rows centred on the empirical argmin, locating the vertex `d_opt = sqrt(1 / x_vert)` (clamping if the vertex falls outside the scan bracket), and re-running `bend_fit_progressive(fitreso_end=d_opt)` once more. The result lands in `scan_dir/best/` with the same per-point outputs as a normal fr-row. A footer in `scan_fitreso.log` records which rows were used, `d_opt`, and the parabola-predicted Rbent.

Empirically (May 2026 reference runs) `d_opt` sits in the 8–13 Å band for the bend-friendly systems and collapses to the coarsest scan point for radiation-damage signal (low-frequency only).

### Obverse/reverse and `altalign.py`

For pairs in R-lattice space groups (R3, R32, R3m, etc., in their hexagonal H setting), the aligning altindex operation can be metric-preserving but not a normalizer of the space group — it flips obverse↔reverse centering. `bendfinder.py` handles this transparently for the scan itself (the conjugate setting is benign as long as both moving and reference end up on the same grid), but downstream refmac on the reindexed moving MTZ may need the data in rhombohedral primitive (R 3 2 :R) form.

`altalign.py` is a standalone tool that does the LSQ-based altindex+origin search and, for non-normalizing R-lattice ops, writes *both* output settings: `<stem>_H32.{pdb,mtz}` (hexagonal axes preserved) and `<stem>_R32R.{pdb,mtz}` (rhombohedral primitive, gemmi+refmac native).

```
ccp4-python altalign.py mov.pdb ref.pdb [out.pdb] [mov.mtz] [out.mtz]
```

Verified on the porin 3poq→3pou pair: refmac rigid completes at R=0.37 on the R 3 2 :R output.
