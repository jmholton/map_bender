# map_bender / bendfinder

Finds a smooth, periodic coordinate-shift field that maps the coordinates of one crystal onto a reference crystal of the same protein. The shift field is expressed as a truncated Fourier series using real Miller-index sine/cosine functions — physically motivated basis functions that respect the periodicity of the crystal lattice.

Primary author: James Holton

## What it does

Given two PDB files of the same protein — same crystal form at different conditions (humidity, temperature, ligand), or crystals so non-isomorphous as to appear to be different forms but aren't — `bendfinder.py` computes a smooth 3D vector field **Δr(x,y,z)** such that applying that field to the coordinates of `bendme.pdb` minimises the all-atom RMSD to `reference.pdb`. Optionally, any CCP4 map in the frame of `bendme.pdb` can be spline-interpolated into the reference frame.

The shift field is parameterised as:

```
Δd(x,y,z) = Σ_{hkl}  A_{hkl} · sin(2π(hx+ky+lz))  +  B_{hkl} · cos(2π(hx+ky+lz))
```

where the sum runs over (h,k,l) triplets sorted by resolution (low resolution first), so the most physically meaningful large-scale deformations are captured first. Symmetry constraints are applied. I.E. in space group P4₃2₁2 (8 proper rotations), 247 canonical parameters control 1401 Friedel-unique coefficients.

## Quick start

```python
from bendfinder import bend_fit_progressive, bend_apply_pdb, bend_apply_map

result = bend_fit_progressive('bendme.pdb', 'reference.pdb')
bend_apply_pdb('bendme.pdb', 'reference.pdb', result, outpath='bent.pdb')
bend_apply_map('bendme.map', result, outpath='bent.map')
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

### Latest version (`bendfinder.py`)

Numbers below come from the 2026-06-16 full-gamut run on one octamus1 node
(`./run_all_tests.com` with 64-CPU pthreaded OpenBLAS for the SVDs).
Pre-bend RMSD/Rfac are the `hkl00` baseline (no shift field applied,
moving map resampled onto reference grid only). Best RMSD/Rbent are
from the parabola-vertex re-fit at `d_opt`. Top peak revealed is the
largest |σ| feature in the post-bend difference map that appears
consistently across the fr-row scan — the structural finding the
shift field uncovers.

| System | Datasets | Pre-bend RMSD | Pre-bend Rfac | Best RMSD | Best Rbent | Top peak revealed | Wall |
|--------|----------|---------------|---------------|-----------|------------|-------------------|------|
| Lysozyme (humidity) | 3aw6 → 3aw7 | 0.76 Å | 53.3% | **0.071 Å** | **30.2%** | −5.2σ A/15HIS/CD2 (side-chain rearrangement) | 26 min |
| DHFR | 1rx2 → 1rx1 | 0.45 Å | 43.1% | **0.176 Å** | **37.8%** | +9.1σ A/161FOL/C14 (folate ligand) | 4 min |
| Myoglobin | 1mbo → 1a6m | 0.31 Å | 55.7% | **0.062 Å** | **50.0%** | −41.5σ A/155HEM/FE (heme iron) | 6 min |
| Lysozyme raddam | 5kxk → 5kxl | 0.12 Å | 11.0% | **0.114 Å** | **11.5%** | −11.1σ A/115CYS/SG (disulfide damage) | 18 min |
| Lysozyme raddam | 5kxk → 5kxm | 0.08 Å | 9.4%  | **0.080 Å** | **9.9%**  | +5.5σ A/105MET/SD (Met oxidation) | 19 min |
| Lysozyme raddam | 5kxk → 5kxn | 0.11 Å | 18.1% | **0.099 Å** | **17.7%** | −14.7σ A/94CYS/SG (disulfide damage) | 26 min |
| Insulin T→R | 4fg3 → 4e7u | 2.31 Å | 85.5% | **0.678 Å** | **64.4%** | −46.7σ D/101ZN/ZN (T→R Zn shift) | 16 min |
| Lipoxygenase | 9o4s → 9o4t | 1.08 Å | 63.6% | **0.206 Å** | **55.0%** | +8.0σ A/91MET/SD (Met sidechain) | ~1 h |
| Porin | 3poq → 3pou | 3.14 Å | 73.6% | **0.369 Å** | **57.8%** | +10.3σ A/244PHE/CB (Phe sidechain) | ~72 min |

What the columns mean:

- **Pre-bend**: the moving map resampled onto the reference grid with
  zero shift (`hkl00` row of the scan). This is the baseline that
  bendfinder has to improve on.
- **Best**: the result of the parabola-vertex re-fit at the optimal
  fitting resolution `d_opt` (see [Best d_opt parabola fit](#best-d_opt-parabola-fit)
  below). Inside the scan log this is the `best` row, computed by
  fitting `Rbent` vs `1/d²` across the `fr20`…`fr5` rows and re-running
  `bend_fit_progressive` at the vertex.
- **Top peak revealed**: dominant |σ| feature visible in the post-bend
  difference map across the scan — the structurally meaningful peak
  the shift field uncovers. Sign convention is `subtract='ref'` by
  default (positive σ = density present in bent but absent from ref);
  raddam runs use `subtract='bent'` (positive = density appearing
  with dose).
- **Wall**: total time for that test's slot in the gamut.

Read each row as "after bending, how close did we get and what was the
first thing the difference map could no longer hide". The radiation-
damage systems clamp at `d_opt = 20 Å` (the coarsest scan point) — the
damage signal is purely low-frequency, so finer HKLs only add noise.
Insulin's T→R transition exceeds the smooth-PSDVF model (LEU B6 shifts
~8 Å between T and R), so no fitreso choice helps and the residual is
honestly the T→R Zn site. Lipoxygenase is the same crystal habit as
the reference (same SG, same general cell) at extreme non-isomorphism
— so much so that the deposited cells look like a form change but
aren't; the pipeline stretches the moving cell into the reference,
picks up an alt-index 180°-about-z, and rigid-body re-refines before
bending. Porin is the obverse↔reverse R-lattice case; the in-bendfinder
altindex pass picks the right 2-fold and the scan runs cleanly
afterwards.

Insulin, lipoxygenase, and porin all need `fill_fcalc=True` because
their deposited MTZs are below 99 % SG-ASU complete; without it
refmac inherits the gaps and `bent.mtz` shows missing-HKL chunks in
Coot. For porin you can also run `altalign.py` directly to get a
refmac-ready R 3 2 :R output — see
[Obverse/reverse and altalign.py](#obversereverse-and-altalignpy)
below.

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
| `run_all_tests.com` | Full 11-test gamut runner (test_symm + magdoff + 8 example scans + porin altalign+refmac) |
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

## Altindex and origin search

Before the shift field can do anything useful, the two crystals have to
be put in the same setting: same fractional origin, same point-group
orbit, and (if the cells differ even slightly) the same metric. PDB
deposits routinely violate all three. The same protein, redetermined in
the same crystal habit, can land at a different cell-origin choice, a
chain-letter swap that's actually just an in-SG symop, an alternate
indexing choice that picks one Laue mate over another, or a few-percent
cell rescaling that makes two near-isomorphous datasets look like
different forms. Bendfinder ships an internal `resolve_altindex` step
that handles all four; `altalign.py` is the standalone equivalent.

### Algorithm

Both paths share the same enumeration kernel
(`_enum_alt_rot_origin_candidates` in `bendfinder.py`):

1. Match CA atoms by (chain, residue, atom name) → arrays `A`
   (moving) and `B` (reference); compute pre-fit RMSD and the Kabsch
   LSQ rigid-body floor as cross-checks.
2. Enumerate every integer basis-change matrix `M` (entries in
   `{−1, 0, 1}`, `|det M| ≤ 1`), build the alt-cell metric
   `G_alt = M G Mᵀ`, look up catalog space groups of that crystal
   system, and conjugate each catalog op back to the original frame as
   `R = Mᵀ R_alt M⁻ᵀ`, `t = Mᵀ t_alt`. Keep only candidates whose
   `R` is integer in `{−1, 0, 1}` with denominators in `{1, 2, 3}` for
   `t`, and whose `R` preserves the cell metric.
3. Rank each candidate by post-transform CA RMSD.

`resolve_altindex` then picks one of four actions on the rank-1 op:

| Action | Trigger | What it does |
|--------|---------|--------------|
| `none` | no candidate beats `0.7 × baseline` | leaves mov as-is |
| `origin_only` | `R = I` and \|t\| > 0.01 | translates PDB; phase-shifts MTZ by `exp(−2πi h·t)` (no re-refinement) |
| `sg_op_origin` | `R` is already a symop of the SG | applies `(R, t)` to PDB cartesian and to MTZ via the SF theorem `F'(h) = F(R^T h)·exp(2πi h·t)` (no re-refinement; \|F\| is SG-symmetric so the lookup is a no-op on amplitudes) |
| `altindex_refine` | `R ≠ I` and `R ∉ SG` | applies `R, t_cart_opt` to PDB cartesian; reindexes MTZ Fobs by `R`; re-refines with refmac to regenerate map coefficients |

Two extensions matter for non-trivial pairs:

- **Cell-stretch pre-step.** When the moving and reference cells
  don't match (`_cells_match` returns False), `resolve_altindex` first
  re-orthogonalizes the moving model into the reference cell —
  fractional coordinates preserved — and relabels the moving MTZ's
  cell record. The isomorphous distortion is absorbed as a uniform
  elastic stretch, leaving only the genuine alt-indexing rotation +
  origin difference for the enumerator. The moving experimental Fobs
  is preserved through this step.
- **Loose metric tolerance for stretched pairs.** Same-cell pairs
  enumerate with `metric_tol_rel = 1e-6`. After a cell stretch the
  tolerance loosens to 5 %, surfacing operations that are
  metric-preserving in a nearby higher-symmetry holohedry the
  deposited cell only approximates. For lipoxygenase this is exactly
  the alt-cell 180°-about-z that the strict 1e-6 tolerance would
  reject.

### What the gamut pairs actually need

For most pairs the answer is trivial — the deposits are already
aligned. For four of the nine they aren't:

| Pair | SG | Action | Op chosen | What was hiding |
|------|----|--------|-----------|-----------------|
| Lysozyme 3aw6 → 3aw7 | P 4₃2₁2 | `origin_only` | `t = (−0.007, −0.003, +0.007)` | Sub-cell origin offset from a humidity-driven 2.5 % cell contraction. |
| Insulin 4fg3 → 4e7u | H 3 | `sg_op_origin` | `R ∈ H 3 symops`, `t ≈ (0.02, 0.006, ⅓)` | The two deposits chose different H 3-equivalent ASUs along the 3-fold axis. SF theorem applies the (R, t) cleanly without re-refining. |
| Porin 3poq → 3pou | H 3 2 | `altindex_refine` | `R = [[0,−1,0],[−1,0,0],[0,0,−1]]`, `t = (⅓, −⅓, 0)` | The obverse↔reverse 2-fold. Metric-preserving but **not** a normalizer of H 3 2, so the reindexed moving crystal lands in the conjugate (reverse-H) setting — benign for the scan but awkward for refmac. See "altalign.py and the dual-solution writer" below. |
| Lipoxygenase 9o4s → 9o4t | P 2₁ | cell-stretch + `altindex_refine` | `R = diag(−1, −1, +1)`, `drot ≈ 0.00°` from LSQ | Same crystal, just very distorted — same SG and same general cell, but the metric drifts by enough (a 92.0→96.0, b 93.0→94.5, c 49.0→50.5, β 92.7→91.2°) that the deposited cells look like a form change but aren't. The monoclinic cell is nearly orthorhombic, so the alt-cell 180°-about-z is metric-preserving only at 5 % tolerance; once the cell is stretched into the reference metric that op aligns mov to ref at drot=0.00°. Moving Fobs is preserved through the reindex + refmac rigid-body step. |

For DHFR, the three radiation-damage pairs, and myoglobin the search
ran but found no candidate that beats 0.7 × baseline (DHFR / raddam:
already aligned, action `none`), or returned identity (myoglobin: cells
differ slightly but no rotation or translation needed beyond the
stretch).

### `altalign.py` — standalone diagnostic and dual-solution writer

`altalign.py` runs the same enumeration kernel as a standalone tool
and prints the top-N ranked candidate list, the LSQ floor, the
chosen op, and (for non-normalizing R-lattice ops) writes dual outputs.

```
ccp4-python altalign.py mov.pdb ref.pdb [out.pdb] [mov.mtz] [out.mtz]
```

Use it when you want a diagnostic view of the candidate ranking
(altalign emits the full top-15 with discrete RMSD, drot vs LSQ, and
comfit RMSD per candidate) or when you need the dual-solution writer
for R-lattice non-normalizing ops. For porin and similar
obverse/reverse cases, altalign emits both `<stem>_H32.{pdb,mtz}`
(hexagonal axes preserved; data is in the conjugate reverse-H setting,
needs `mtzutils` to swap SYMM records before refmac) and
`<stem>_R32R.{pdb,mtz}` (rhombohedral primitive, gemmi + refmac
native — setting-unambiguous). Refmac rigid completes at R=0.37 on
the porin R 3 2 :R output.

For the same-crystal-just-distorted case (lipox-style: same SG and
general cell, but with enough metric drift that strict 1e-6
metric-preservation tolerance rejects the alt-cell op), `altalign.py`
does **not** do the cell-stretch pre-step or relax that tolerance, so
it will miss the op that `resolve_altindex` finds. Run `fitreso_scan`
directly and read the `resolve_altindex:` block in `scan_fitreso.log`.
