# map_bender / bendfinder — Claude context

## Project overview

`bendfinder.py` fits a smooth, periodic 3D coordinate-shift field between two non-isomorphous crystal forms of the same protein. The shift field is a truncated Fourier series indexed by Miller indices (h,k,l).

The git repo lives at `./map_bender/` (relative to `../`). The working development copy is `../bendfinder.py`. **Keep both in sync when making changes — they must stay identical.**

`bendfinder.com` (tcsh prototype) is historical reference only; development is Python-only.

## Directory structure

```
../                             working area, not in git
  bendfinder.py                 current development copy (keep in sync with map_bender/)
  bendfinder.com                prototype tcsh script (historical reference)
  map_bender/                   git repository
    bendfinder.py
    CLAUDE.md                   this file
    README.md
    origins.com
    LICENSE
    examples/3aw6_3aw7/         canonical lysozyme example data (checked into git)
  lyso/                         lysozyme 3aw6→3aw7 (P4₃2₁2, 1008 CA pairs)
    3aw6.pdb, 3aw7.pdb
    3aw6_2fofc.map, 3aw7_2fofc.map
    fitreso_scan/               hkl00..hkl10, fr5..fr20 output subdirs
  dhfr/                         DHFR 1rx1→1rx2 (P2₁2₁2₁, 636 CA pairs)
    fitreso_scan/
  myoglobin/                    myoglobin 1mbo→1a6m (P2₁, 294 CA pairs)
    fitreso_scan/
  raddam/                       radiation damage 5kxk→5kxl/m/n (P4₃2₁2, ~980 CA pairs)
    fitreso_scan_5kxl/, fitreso_scan_5kxm/, fitreso_scan_5kxn/
  insulin/                      insulin hexamer 4fg3→4e7u (H3, T→R transition)
    4fg3.pdb, 4e7u.pdb
  lyso_test_031419/             gold-standard reference run (old prototype, RMSD=0.209 Å)
  magdoff/                      Magdoff synthetic deformation validation tests
    test_magdoff.py             test script (7rsa, P2₁, 248 CA)
    7rsa.pdb                    reference structure (ribonuclease A, 1.26 Å)
```

To run a scan, call `fitreso_scan()` directly from `bendfinder.py`, e.g.:
```python
cd raddam && ccp4-python -c "
import sys; sys.path.insert(0, '..')
from bendfinder import fitreso_scan
fitreso_scan(mov_pdb='5kxk.pdb', ref_pdb='5kxl.pdb',
             mov_map='5kxk_fc.map', ref_map='5kxl_fc.map',
             scan_dir='fitreso_scan_5kxl')
"
```
Per-system `*_fitreso_scan.py` wrapper scripts have been removed.

## Algorithm (Python version)

1. Both PDBs are expanded from their crystallographic space group to P1 via gemmi.
2. CA atom pairs are matched by residue+atom name. Outliers rejected by shift magnitude (`outlier_sigma`) and B-factor (`b_sigma`) using robust statistics (MAD-based).
3. All (h,k,l) out to `fitreso_end` Å are enumerated, sorted low-resolution first, deduplicated for Friedel symmetry.
4. Proper-rotation symmetry of the space group is applied (`use_symm=True`): each Friedel-unique HKL is assigned a canonical representative, reducing free parameters by the point-group order (×8 for P4₃2₁2).
5. A joint design matrix **X** (3N_atoms × 6M_canon) is built from symmetry-expanded sine/cosine basis. All three axes are fitted simultaneously.
6. **X @ params = shifts** is solved by SVD (scipy.linalg), giving globally optimal (A, B) coefficients and covariance-matrix uncertainties.
7. Steps 3–6 repeat coarse-to-fine (`fitreso_start` → `fitreso_end`), admitting `batch_hkls` new HKLs per iteration. The loop stops when the overdetermination ratio (atoms / canonical params) drops below `od_margin`.
8. Optional `iter_callback(iter_i, nhkls, n_canon, rmsd, hkls, AB_xyz, active, snr)` is called after each iteration. `max_canon` stops the loop once n_canon ≥ that value.

## Key parameters

| Parameter | Default | Role |
|-----------|---------|------|
| `fitreso_start` | 20.0 | Coarsest resolution to admit (Å) |
| `fitreso_end` | 7.0 | Finest resolution to admit (Å) |
| `batch_hkls` | 100 | New HKLs admitted per iteration |
| `od_margin` | 1.5 | Stop when atoms/(canon params) < od_margin |
| `outlier_sigma` | 2.5 | CA pair rejection threshold (robust σ of shift magnitude) |
| `b_sigma` | 3.0 | CA pair rejection threshold (B-factor; median + b_sigma×σ_MAD) |
| `drop_snr` | 0.0 | Drop HKLs with \|A,B\|/σ < drop_snr (0 = keep all) |
| `use_symm` | True | Apply space-group proper-rotation constraints |
| `dimensions` | 'xyz' | Which coordinate axes to fit |
| `iter_callback` | None | Called after each progressive iteration (see above) |
| `max_canon` | None | Stop after first iteration where n_canon ≥ this value |

## Fitreso scan

The scan logic is implemented as `fitreso_scan()` directly in `bendfinder.py`. Each system's `*_fitreso_scan.py` script is a thin wrapper that sets paths and calls `fitreso_scan()`.

`fitreso_scan(mov_pdb, ref_pdb, mov_map, ref_map, scan_dir, ...)` runs three sections and reports per-point: label, RMSD(CA) before/after, active HKLs, Rfac (Riso vs reference map), strongest positive and negative difference map peaks with nearest atom.

**Section 1 — hkl00**: zero shift; just resample the moving map onto the reference grid via `interpolate_map`. Establishes the unregistered baseline.

**Section 2 — hkl01..hkl10**: single `bend_fit_progressive` call with `batch_hkls=1, max_canon=11, fitreso_start=100`. The `iter_callback` fires once per canonical HKL added (n_non_dc = n_canon − 1). Efficient: only one initialization (atom matching, outlier rejection) for all 10 checkpoints.

**Section 3 — fr20..fr5**: separate `bend_fit_progressive` calls with `fitreso_end` in [20,15,12,10,8,7,6,5] Å.

**Riso calculation**: `compute_riso(ref_mtz, ref_col, test_mtz, test_col)` in `bendfinder.py`. Uses Wilson isotropic B scaling: fit `log(F_ref/F_test) = log(kF) − (B/4)·(1/d²)` by OLS, then Riso = Σ|F_ref − kF·exp(−B/4d²)·F_test| / Σ|F_ref|. Returns (riso, kF, B_iso). No CCP4 programs needed; replaces the old `diff.com` subprocess.

**Output files per scan point**:
- `cootme.mtz` — columns `FDM`/`PHIDM` (bent map) + `DELFWT`/`PHDELWT` (diff map). Load in Coot; FDM/PHIDM are recognised by default.
- `diff_norm.map` — z-scored real-space difference map (ref − bent) / σ; positive peaks = more density in reference.
- `bent.map` — bent moving-crystal map resampled on the reference grid.
- `PSDVF.mtz` (fr\* points only) — fitted shift-field (h,k,l) coefficients.

**Diff map sign convention**: `diff = ref_normalised − bent_normalised`. Positive peaks = more density in reference. For the raddam series where the damaged structure is the reference, positive peaks indicate features appearing with radiation dose.

**Raddam sign convention**: 5kxk (undamaged, lowest dose) is the **moving** model; 5kxl/5kxm/5kxn (increasingly damaged) are the **references**. Dose ordering is alphabetical: 5kxk < 5kxl < 5kxm < 5kxn. With this convention, positive diff peaks = features appearing with dose, negative peaks = features disappearing.

**Large-map memory**: `eval_shift_field` allocates an (N_voxels × N_hkls) phase matrix. For large maps (raddam: 3.7M voxels × 900+ HKLs ≈ 26 GB) this must be chunked. `fitreso_scan` uses an internal `_eval_chunked` helper in 50k-voxel batches (configurable via `chunk_size` parameter).

## Origin alignment (`_find_best_origin`)

`_find_best_origin(atoms1_op0, atoms2, sg_name, n_polar=12)` returns a 4-tuple `(shift, symop_k, R_alt, improved)`.

**Standard search**: For each space-group symop k and each polar-axis candidate shift (an `n_polar`-interval grid along each allowed origin-shift direction from `_ORIGINS_TABLE`), compute the median fractional CA shift between the op0 atoms of structure 1 and the op-k atoms of structure 2 after applying the trial shift. Take the minimum-median (shift, k) pair; call it `improved` if it reduces the unshifted score by >10%.

**Altindex fallback**: Only triggered when the standard search fails to improve (score ≥ 0.9 × ref_score). Tries proper rotation matrices from the Laue-group holohedry that are NOT in the space group's point group (defined in `_ALTINDEX_CANDIDATES` per crystal system). For each candidate R_alt, applies it to the fractional coordinates of atoms2 and reruns the standard search. If any altindex+shift combination yields a better score, that R_alt is applied permanently to atoms2 before fitting.

```python
_ALTINDEX_CANDIDATES = {
    'trigonal':    [2-folds about [100]/[010]/[110] in hex fractional],
    'hexagonal':   [same],
    'tetragonal':  [[0,1,0],[1,0,0],[0,0,-1]],
    'orthorhombic': [axis permutations],
    'cubic':       [3-fold permutations],
    'monoclinic':  [],
    'triclinic':   [],
}
```

Callers (`bend_fit`, `bend_fit_progressive`) unpack all 4 values and print: `origin shift: (x, y, z) for SG  +altindex  +symop_idx=k`.

**CA shift sanity checks** (both in `bend_fit` and `bend_fit_progressive`): After outlier rejection, two checks run:
1. **CA RMSD > 5 Å** — RMS orthogonal shift over all remaining CA pairs. Printed as `CA RMSD after origin fix: X.XXX Å`. Raises error if > 5 Å; a correctly aligned pair should always be below this even for large conformational changes. (Insulin T→R: 2.417 Å.)
2. **Largest fractional CA shift > 0.35 cells** (~16 Å for a 47 Å cell). Raised from the old 0.1 threshold to accommodate genuine large differences (e.g. insulin LEU B6, 0.165 fractional ≈ 8 Å). Values >0.35 indicate gross misalignment.

## Cubic b-spline boundary fix (interpolate_map)

`scipy.ndimage.map_coordinates` with `order=3, mode='wrap'` can produce overshoot at cell boundaries if adjacent grid rows have large density gradients. This creates spikes in difference maps at x=0, y=0, z=0.

All test-system maps are ASU maps (not full-cell): lysozyme covers ~½ cell in each axis; DHFR covers full X and Z but ¼ of Y; myoglobin is full-cell. For ASU maps, `mode='wrap'` wraps the boundary (e.g. x≈0.5) to x=0, which is a physically unrelated density — creating a flat artefact in the difference map near the ASU boundary.

**Fix**: `interpolate_map` pads by 5 voxels on each side with `np.pad(data, 5, mode='reflect')` before calling `map_coordinates` with `mode='nearest'`. Reflection gives a smooth continuation at any boundary (ASU edge or unit-cell edge), so the IIR prefilter sees no discontinuity. The prefilter decays as 0.268^k, so 5-voxel padding gives <0.15% boundary error at any interior query point. Do not use `pad=2` (insufficient) or `order=1` (eliminates overshoot at the cost of cubic quality everywhere).

## gemmi map→MTZ conversion

The `map2mtz` function in scan scripts uses:

```python
ccp4 = gemmi.read_ccp4_map(mapfile)
ccp4.setup(float('nan'))
sf   = gemmi.transform_map_to_f_phi(ccp4.grid, half_l=True)
data = sf.prepare_asu_data(dmin=0.0)
mtz  = gemmi.Mtz(with_base=True)
mtz.spacegroup = ccp4.grid.spacegroup
mtz.cell       = ccp4.grid.unit_cell
mtz.add_dataset('dataset')
mtz.add_column('F', 'F')
mtz.add_column('PHI', 'P')
...
mtz.write_to_file(mtzfile)
```

This replaces the old sfall subprocess approach (which failed for P2₁ maps). No CCP4 programs needed for this step; gemmi handles axis ordering from the CCP4 header automatically.

## Space-group generality and testing

### Symmetry constraint

The PSDVF must satisfy the crystallographic symmetry constraint for all proper operators {R_k, t_k} of the space group:

```
Δr(R_k x + t_k) = R_k · Δr(x)
```

where R_k and t_k are in fractional coordinates. `use_symm=True` enforces this by construction via `build_design_matrix_symm` and `expand_ab_canon`.

**Critical implementation detail — R^{−1} vs R^T:** In fractional coordinates, proper rotation matrices R_k are NOT in general orthogonal (R_k^T ≠ R_k^{−1}). The correct vectorial transform for the shift field is R_k^{−1}, not R_k^T. For orthogonal crystal systems (cubic, tetragonal, orthorhombic, monoclinic) R_k is orthogonal so R_k^T = R_k^{−1} and the distinction doesn't matter. For **trigonal and hexagonal** systems the 3- and 6-fold fractional rotation matrices are not orthogonal, making this distinction essential. An earlier implementation used R_k^T throughout and gave violations of 1–100 Å for all trigonal/hexagonal SGs. The fix: `np.linalg.inv(R)` in both `build_design_matrix_symm` and `expand_ab_canon`.

### test_symm_all_sgs.py

`claude/test_symm_all_sgs.py` verifies the symmetry constraint across all 65 Sohncke (protein-compatible) space groups. Protocol per SG:

1. Place 8 atoms at random general positions in the ASU.
2. Apply small random fractional displacements δ_i to each ASU atom.
3. Expand both reference and displaced states to P1 via all SG operators — displacement of op-k copy of atom i = R_k · δ_i (the correct crystallographic scenario: both structures in the same SG).
4. Fit unconstrained (all P1 atoms, `build_design_matrix`) and constrained (ASU atoms only, `build_design_matrix_symm`).
5. Evaluate the symmetry violation max|Δr(R·x+t) − R·Δr(x)| in Å on a 15³ grid.

**HKL set must be orbit-complete.** The full `generate_hkls(cell, FITRESO)` set (no truncation) is used. The resolution sphere is preserved under any point-group rotation (d-spacing is invariant), so the full set at a given resolution is automatically orbit-complete. Truncating to a fixed number of HKLs breaks this closure for trigonal/hexagonal and produces spurious violations.

**Results (after fix):**
- Constrained fit: **65/65** SGs have violation < 10⁻⁶ Å (max ~3 × 10⁻¹³ Å — machine precision)
- Unconstrained fit: **0/65** SGs have violation > 10⁻⁶ Å when the HKL set is orbit-complete and the data is SG-consistent

The unconstrained fit also gives a SG-symmetric field in this controlled test because: (a) the input data is SG-consistent by construction, and (b) the orbit-complete HKL set makes the SG-symmetric solution the unique minimum-norm SVD result.

### Centering translations (I, F, C, R)

`get_proper_symops` accepts all operators with det(R) = +1, including centering translations (R = I, t = centering vector). For I 2 (SG 5, I-centered monoclinic) this correctly adds operators with t = (½,½,½), imposing Δr(x + ½,½,½) = Δr(x) — which in Fourier space restricts to h+k+l even (the I-centering systematic absence). Tested explicitly: I 1 2 1 passes at 4.6 × 10⁻¹⁵ Å.

## Magdoff synthetic deformation tests (`magdoff/test_magdoff.py`)

Controlled validation on ribonuclease A (7rsa, P2₁, a=30.18 b=38.40 c=53.32 Å β=105.85°, 124 residues, 248 CA in P1). Each test imposes a known deformation and measures how well bendfinder recovers it. Riso is computed at 1.5 Å using gemmi `DensityCalculatorX` → `transform_map_to_f_phi`.

**Important distinction:** Riso is evaluated at **1.5 Å data resolution**; the PSDVF is fitted to **7 Å spatial resolution**. These are separate concepts. Residual Riso after bending reflects high-frequency content (1.5–7 Å band) that lies beyond the PSDVF bandwidth — not a failure of the shift field within its design range.

### Test 1 — Isomorphous cell change (true Magdoff)

Scale a, b, c in CRYST1 by 1.005; leave atom Cartesian coordinates unchanged. Simulates a 0.5% crystal expansion: same protein structure, slightly different unit cell. The fractional coordinates shift by ~0.5% because the same Cartesian positions map to different fractions under the new metric.

**Critical:** SCALE and ORIGX records must be stripped from the modified PDB. If left, gemmi reads the stale fractional matrix from those cards rather than deriving it from the modified CRYST1, so the fractional coordinates are unchanged and no deformation is seen.

| | RMSD(CA) | Riso (1.5 Å) |
|---|---|---|
| Before bending | 0.197 Å | 15.0% |
| Constrained (P2₁, 301 HKLs) | 0.018 Å (91.1% recovery) | 2.7% |
| Unconstrained (201 HKLs) | 0.015 Å (92.5% recovery) | 2.7% |

### Test 2 — Rigid-body rotation 0.5°

Rotate all atom Cartesian coordinates by 0.5° about a fixed random axis (seed 42 → axis [0.231, −0.789, 0.569]). Not SG-consistent for a general axis.

| | RMSD(CA) | Riso (1.5 Å) |
|---|---|---|
| Before bending | 0.335 Å | 24.0% |
| Constrained (P2₁, 301 HKLs) | 0.028 Å (91.8% recovery) | 4.3% |
| Unconstrained (201 HKLs) | 0.026 Å (92.2% recovery) | 4.3% |

Both tests fit 3 progressive iterations (20→7 Å), completing in ~25 s per fit on a single CPU. Constrained and unconstrained give nearly identical recovery for P2₁ (order 2), though constrained reaches 50% more HKLs before hitting the overdetermination limit.

## Empirical results (fitreso scans)

All systems use default parameters (`outlier_sigma=2.5`, `b_sigma=3.0`, `drop_snr=0`, `batch_hkls=100`).

| System | Space group | CA pairs | fr5 RMSD | fr5 Rfac | hkl00 Rfac |
|--------|------------|----------|----------|----------|------------|
| Lyso 3aw6→3aw7 | P4₃2₁2 | 1008 | 0.034 Å | 33.2% | 59.5% |
| DHFR 1rx1→1rx2 | P2₁2₁2₁ | 636 | 0.071 Å | 41.9% | 44.1% |
| Myoglobin 1mbo→1a6m | P2₁ | 294 | 0.063 Å | — | — |
| Raddam 5kxk→5kxl | P4₃2₁2 | 976 | 0.082 Å | 20.6% | 13.4% |
| Raddam 5kxk→5kxm | P4₃2₁2 | 984 | — | — | 11.0% |
| Raddam 5kxk→5kxn | P4₃2₁2 | ~992 | — | — | — |
| Insulin 4fg3→4e7u | H3 | 534 | 1.063 Å | 68.7% | 83.7% |

Notes:
- Lyso Rfac plateaus at ~33% by fr20 and barely changes with higher resolution — real structural differences remain in the diff map (A/74ASN/O is the persistent −10σ peak).
- DHFR Rfac barely improves (44% → 42%) — these crystal forms are more dissimilar; the FOL ligand and Mn ion dominate the diff map at all resolutions.
- Raddam Rfac *starts* low (13–11%) because the fc maps are nearly identical; huge water peaks (±30–65σ) reflect water molecules appearing/disappearing with radiation dose.
- Myoglobin Rfac pending (gemmi map2mtz re-run in progress); heme iron dominates diff map throughout (±35σ at A/154HEM/FE).
- Insulin: high Rfac (68.7%) and high residual RMSD (1.063 Å) reflect genuine T→R conformational change — LEU B6 shifts ~8 Å between T-state (4fg3) and R-state (4e7u), which is outside the smooth shift-field model. RMSD best at fr10 (0.988 Å, 401 HKLs); OD limit hit at 501 HKLs for fr8–fr5. Dominant diff peaks: +34σ at D/101ZN/ZN (zinc position differs), −10σ at waters/SCN.
