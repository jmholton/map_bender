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
  dhfr/                         DHFR 1rx2→1rx1 (P2₁2₁2₁, 636 CA pairs).
                                 Canonical direction: mov=1rx2 (has FOL +
                                 Mn²⁺ + BME), ref=1rx1 (just NAP + Ca²⁺).
                                 With subtract=ref this makes FOL appear as
                                 positive (green) density and Ca²⁺ as
                                 negative (red).  1rx1→1rx2 also works but
                                 flips the colours.
    fitreso_scan/
  myoglobin/                    myoglobin 1mbo→1a6m (P2₁, 294 CA pairs)
    fitreso_scan/
  raddam/                       radiation damage 5kxk→5kxl/m/n (P4₃2₁2, ~980 CA pairs)
    fitreso_scan_5kxl/, fitreso_scan_5kxm/, fitreso_scan_5kxn/
  insulin/                      insulin hexamer 4fg3→4e7u (H3, T→R transition)
    4fg3.pdb, 4e7u.pdb          NB: both deposited MTZs are < 99% SG-ASU
                                 complete (4e7u 93%, 4fg3 97.6%) — pass
                                 `fill_fcalc=True` (CLI `--fill-fcalc`)
                                 to fitreso_scan or run_refinement will
                                 sys.exit(1).
  porin/                        porin 3poq→3pou (H 3 2)
    3poq.pdb, 3pou.pdb          NB: 3poq.mtz is 89.9% complete
                                 (3pou.mtz 99.6%) — needs
                                 `fill_fcalc=True`.  3poq/3pou are an
                                 obverse/reverse pair — see the
                                 "altalign.py and the porin
                                 obverse/reverse problem" section
                                 below.  The PDB side is solved
                                 (altalign aligns to ~2 Å); the moving
                                 MTZ needs the dual-solution writer.
  lipox/                        soybean lipoxygenase-1 9o4s→9o4t (P2₁,
    9o4s.pdb, 9o4t.pdb           ~800 CA, ~4% non-isomorphous cell
                                 expansion — same SG, different cells
                                 a 92.0→96.0, b 93.0→94.5, c 49.0→50.5,
                                 β 92.7→91.2°).  Raw cif2mtz outputs
                                 are I-only and < 99% SG-ASU complete,
                                 so `fill_fcalc=True` is required and
                                 `run_refinement` auto-ctruncates I→F
                                 before refmac.  The mov→ref relation
                                 is a non-crystallographic ~180°
                                 rotation about Cartesian z, but the
                                 monoclinic cell is *nearly*
                                 orthorhombic so a 5%-loose metric
                                 tolerance surfaces the alt-cell 2-fold
                                 as a discrete altindex op — see
                                 "Cross-cell pairs" section below.
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

**Section 4 — best**: parabolic fit of Rbent vs 1/d² across the fr-rows; one final `bend_fit_progressive` at the vertex `d_opt`. See [Best d_opt parabola fit](#best-d_opt-parabola-fit) below.

**Riso calculation**: `compute_riso(ref_mtz, ref_col, test_mtz, test_col)` in `bendfinder.py`. Uses Wilson isotropic B scaling: fit `log(F_ref/F_test) = log(kF) − (B/4)·(1/d²)` by OLS, then Riso = Σ|F_ref − kF·exp(−B/4d²)·F_test| / Σ|F_ref|. Returns (riso, kF, B_iso). No CCP4 programs needed; replaces the old `diff.com` subprocess.

**Output files per scan point**:
- `bent.mtz` — columns `FDM`/`PHIDM` (bent map) + `DELFWT`/`PHDELWT` (diff map). Load in Coot; FDM/PHIDM are recognised by default.
- `diff_norm.map` — z-scored real-space difference map (sign controlled by `subtract`); default `subtract=ref` → diff = bent − ref, positive peaks = density present in bent but absent (or weaker) in ref.
- `bent.map` — bent moving-crystal map resampled on the reference grid.
- `PSDVF.mtz` (fr\* + best points) — fitted shift-field (h,k,l) coefficients.

### Best d_opt parabola fit

Across every example system the Rbent-vs-fitreso curve is a clear U: it
drops as the smooth PSDVF absorbs more structural detail with finer
resolution, then climbs again at the highest resolutions as
high-frequency HKLs add noise outside the shift field's natural
bandwidth.  `fitreso_scan` exploits this by fitting a parabola in
**x = 1/d²** (the natural axis for R-factor-vs-resolution behaviour) to
the 3–5 fr-rows centred on the empirical argmin, locating the vertex
`d_opt = sqrt(1 / x_vert)`, clamping if the vertex falls outside the
bracket, and then re-running `bend_fit_progressive(fitreso_end=d_opt)`
once more.  The result lands in `scan_dir/best/` with the same per-point
outputs as a normal fr-row.  A footer in `scan_fitreso.log` records
which rows were used, `d_opt`, and the parabola-predicted Rbent.

Empirically (May 2026 reference runs) d_opt sits in the 8–11 Å band for
the bend-friendly systems (lyso d_opt = 9.7 Å; raddam similar) and
collapses to whatever the empirical argmin already was for systems
where Rbent saturates quickly (insulin, where the T→R shift exceeds the
PSDVF's smooth-deformation assumption — the curve barely bends, so the
parabola vertex is close to the coarsest scan point).

If fewer than 3 fr-rows have a valid Rbent, or the parabola opens
downward (no interior minimum), the best section silently falls back to
the empirical argmin or is skipped.

**Diff map sign convention**: controlled by the `subtract` parameter on `fitreso_scan` (CLI: `subtract=ref|bent`).
- `subtract='ref'` (default) → diff = bent − ref. Positive peaks = density present in bent (the moving structure resampled into ref's frame) but absent (or weaker) in ref. This is the "what new density does the moving structure bring in?" view.
- `subtract='bent'` → diff = ref − bent. Positive peaks = density present in ref but absent (or weaker) in bent. (Old default before 2026-05-16; flip with this option if you have prior figures/notes using that convention.)

The diff is computed in F-space after k+B scaling of bent → ref (`_fspace_scale_and_diff`), then inverse-FFT'd to `diff.map`/`diff_norm.map` and stored as DELFWT/PHDELWT in `bent.mtz`. The sign flip is exact (negation of every voxel; DELFWT amplitude unchanged, PHDELWT shifts by 180°).

**Raddam sign convention**: 5kxk (undamaged, lowest dose) is the **moving** model; 5kxl/5kxm/5kxn (increasingly damaged) are the **references**. Dose ordering is alphabetical: 5kxk < 5kxl < 5kxm < 5kxn. With the new default `subtract='ref'`: positive diff peaks = features in the **undamaged** model absent from the damaged reference (features disappearing with dose); negative peaks = features appearing with dose. (Sign flipped from the pre-2026-05-16 default — pass `subtract=bent` to recover the old convention.)

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

**CA shift sanity check**: After outlier rejection, if the largest fractional CA shift among remaining atoms exceeds **0.35 cells** (~16 Å for a 47 Å cell), an error is raised. This threshold was raised from 0.1 to accommodate genuine large conformational differences (e.g. insulin T→R, LEU B6 shifts ~8 Å ≈ 0.165 fractional). Values >0.35 indicate gross misalignment or a wrong origin.

## Basis-change altindex enumeration (`_get_altindex_ops`)

`_get_altindex_ops(sg_name, cell, max_M=1, det_max=1)` returns the
altindex candidate ops as `(R_frac, t_frac)` tuples (previously it
returned bare `R` arrays — callers `_find_best_origin` and
`_enum_alt_rot_origin_candidates` were updated to unpack the pair).

The enumeration is ported from `~/projects/origins/claude/gemmi_altindex.py`:
for every integer basis-change matrix `M` (entries in `{-max_M..max_M}`,
`|det M| ≤ det_max`), build the alt-cell metric `G_alt = M G Mᵀ`,
identify its crystal system, look up catalog SGs of that system, and
conjugate each catalog symop back to the original frame as
`R_old = Mᵀ R_alt M⁻ᵀ`, `t_old = Mᵀ t_alt`.  Keep ops where `R_old` is
integer with entries in `{-1,0,1}`, metric-preserving, and `t_old` has
denominator in `{1,2,3}`.  This finds altindex ops the old
`{-1,0,1}`-only direct enumeration missed (e.g. the porin obverse↔reverse
2-fold).  Vectorised + `lru_cache`d; `_catalog_ops_by_cs()` pre-extracts
the catalog once.  ~30 s cold for a large trigonal cell, then free.

A **normalizing** altindex op (one in the normalizer of the SG —
preserves both rotations *and* centering) reindexes the moving MTZ
cleanly: `R·SG·R⁻¹ = SG`.  A **non-normalizing** metric-preserving op
does not — it lands the moving data in a conjugate setting (see porin
below).  For H 3 2, 44 of the 75 enumerated ops normalize.

## Cross-cell pairs (`resolve_altindex` cell stretch + loose-tolerance altindex)

When mov.cell ≠ ref.cell (same SG, different cells — e.g. lipox
9o4s/9o4t with ~4% isomorphous expansion), the strict discrete
`(R_frac, t_frac)` enumeration in `resolve_altindex` can only ever
return lattice translations as the "best" op: no purely-fractional
rotation compensates for a metric change.  Three coupled mechanisms
recover the cross-cell case while keeping mov's experimental Fobs
intact through the entire pipeline:

1. **Cell stretch pre-alignment.**  When `_cells_match(mov, ref)` is
   False, `_stretch_pdb_to_cell` re-orthogonalizes the moving model
   into ref's cell preserving each atom's fractional coordinates, and
   `_relabel_mtz_cell` updates both the file-level and per-dataset
   cell records on the moving MTZ (HKLs and column values untouched).
   The isomorphous distortion is absorbed as a uniform elastic
   stretch; the remaining real-space difference is the genuine
   alt-indexing rotation + translation that the discrete enumeration
   knows how to find.  Port of `~/projects/origins/claude/stretch_to_cell.py`.

2. **COM-optimal continuous translation.**  Discrete origin-table
   entries don't generally contain the continuous COM offset between
   two cross-cell deposits.  For each candidate R, `resolve_altindex`
   now evaluates the RMSD using
   `t_opt = mean(B_frac) − mean(R·A_frac)` (LSQ-optimal continuous t)
   rather than the per-atom-wrap stand-in.  The reported `is_zero_t`
   check is on the UNWRAPPED t so a lattice-vector shift (F-space
   no-op, real-space ~|a| Å translation) actually gets applied to
   the PDB.

3. **Loose metric tolerance.**  `_get_altindex_ops` and
   `_enum_alt_rot_origin_candidates` accept a `metric_tol_rel`
   parameter (default `1e-6` keeps existing same-cell callers
   identical).  `resolve_altindex` sets it to `0.05` when
   stretched=True, which surfaces ops that are metric-preserving in a
   NEARBY higher-symmetry holohedry the deposited cell only
   approximates.  Lipox example: monoclinic P2₁ with a/b within 1.6%
   and β within 1.3° of orthorhombic admits the alt-cell 180°-about-z
   (`R = diag(−1,−1,+1)`) at drot=0.00° from the Kabsch LSQ rotation;
   strict 1e-6 tolerance would reject this op.  The existing
   `altindex_refine` branch then reindexes mov's experimental Fobs by
   R and re-refines; **no Fcalc substitution ever happens — mov's
   experimental data is preserved through the entire pipeline**.

Cross-cell outcome for lipox 9o4s→9o4t:
`action=altindex_refine`, R=diag(−1,−1,+1), drot=0.00°,
rmsd_after=1.43 Å (vs LSQ 1.16 Å), refmac R=0.254 after
reindex+rigid-body.  bend_fit's smooth shift field absorbs the
remaining ~1 Å rigid-body residual + ~5% metric stretch.

Same-cell pairs (lyso, dhfr, etc.) are unaffected by the loose
tolerance default — `_cells_match` returns True and the strict-1e-6
enumeration runs as before.  When cells differ even slightly (lyso
3aw6 vs 3aw7 are ~1.5% off), the new path activates: stretch +
origin_only (R=identity, small continuous t).  Lyso regression
confirmed: fr7 RMSD 0.044 Å matches canonical fr5 0.033 Å within
scan-grid noise; d_opt=9.72 Å matches canonical 9.7 Å; predicted
best Rbent 30.2% matches canonical 30.3%.

## altalign.py and the porin obverse/reverse problem

`altalign.py` is a standalone LSQ-based altindex+origin search (run it
instead of a full `fitreso_scan` while iterating on alignment):

1. Kabsch rigid-body fit of matched CA → continuous `R_lsq`.
2. Enumerate `(altindex × symop × origin)` candidates via
   `_enum_alt_rot_origin_candidates` (the basis-change enumeration).
3. Score each in honest cartesian RMSD (one whole-lattice image shift
   for the whole model, **not** per-atom fractional wrapping — the old
   `resolve_altindex` per-coordinate `diff -= round(diff)` is what made
   it fragile).  `discrete` = exact crystallographic translation;
   `comfit` = COM-optimal continuous translation (the theoretical floor).
4. Rank by honest discrete RMSD; `drot` (deviation from `R_lsq`) is
   carried as a cross-check.

`ccp4-python altalign.py mov.pdb ref.pdb [out.pdb] [mov.mtz] [out.mtz]`

**Porin diagnosis (definitive, May 2026).**  3poq and 3pou are *both*
deposited as standard obverse H 3 2 (every reflection obeys −h+k+l=3n).
altalign finds the op that aligns them — a 2-fold,
`R=[[0,-1,0],[-1,0,0],[0,0,-1]]` — and the moving PDB lands an honest
**2.11 Å** (CA, no fit) from 3pou.  But that op is the obverse↔reverse
2-fold: it is metric-preserving yet does **not** normalize H 3 2 (it
maps centering vector (2/3,1/3,1/3) → (2/3,1/3,2/3), a *reverse*
centering vector).  So the *aligned* moving crystal is unavoidably in
the **reverse H setting**.  This is fundamental — the only clean
(normalizing) altindex op gets porin no closer than 18.86 Å, and the
reverse→obverse converter is a 180° rotation about c that re-scrambles
the fit.  gemmi has no name for reverse-hex H 3 2
(`find_spacegroup_by_ops` → `None`), which is why naively reindexing the
moving MTZ then testing it against obverse H 3 2 reports 30%
completeness (exactly ⅓ — 2-of-3 centering reflections look absent).

**Dual-solution writer (implemented).**  When the chosen op is a
non-normalizer in an R-lattice SG (`sg.ext == 'H'`), altalign writes
*two* output pairs of `<stem>_{H32,R32R}.{pdb,mtz}`.  Neither is
strictly better; which is preferred depends on the downstream workflow:

- **H 3 2 + SYMM**: hexagonal axes preserved; atoms and data reindexed
  by (R, t).  Reference workflow is unchanged.  gemmi cannot name
  reverse-hex H32 (`find_spacegroup_by_ops` → `None`), so the SYMM
  records gemmi writes for the MTZ stay obverse — the data is in the
  conjugate (reverse) setting under obverse labelling.  For refmac on
  this output, swap SYMM records with CCP4 `mtzutils` first.
- **R 3 2 :R**: convert the aligned moving crystal to the rhombohedral
  primitive setting via gemmi's `R 3 2:R` basisop (linear part
  `(-y+z, x+z, -x+y+z)`) composed with the obverse↔reverse fractional
  flip `diag(-1,-1,1)` (so the reverse-hex labels from the reindex
  collapse correctly to integer rhomb HKLs).  Primitive, no centering,
  no obverse/reverse ambiguity; gemmi+refmac native.  Cleanest;
  reference must also be converted to R 3 2 :R for downstream
  comparison (or bendfinder taught to bridge settings) since the rhomb
  cell's gemmi-canonical cartesian frame differs from the hex frame.
  Important: when writing the R32 MTZ, the *per-dataset* cells must
  also be updated (not just `mtz.cell`) — CCP4 reads the dataset cell,
  not the file-level cell, and a mismatch triggers a refmac "Large
  differences between cells from pdb and mtz" abort even though
  `mtz.cell` and `pdb.cell` agree.

Porin R32:R pair (May 2026): refmac rigid completes at R=0.37
Rfree=0.24, completeness 89.9% (same as original obverse-hex data,
since the reflection set is identical, just relabelled).  The H32
solution still needs the CCP4 SYMM swap for refmac.

**Generality (beyond R32).**  The obverse↔reverse centering flip is
unique to R-lattice groups (R3, R-3, R32, R3m, R3c, R-3m, R-3c) in
hexagonal H setting.  The general phenomenon — aligning op is
metric-preserving but NOT in the SG's normalizer, so the reindexed
moving crystal lands in a conjugate setting — is broader: it also
occurs with centered groups (C/I/F/A) where a non-normalizer can flip
the centering type, and with pseudo-symmetric cases (the SG point group
< the lattice holohedry).  A pure-P1 crystal cannot hit this even with
a pseudo-rhombohedral cell because conjugating the trivial group gives
the trivial group.  The dual writer's `_normalizes_sg` detection and
the explicit conjugate-group handling are general; the R 3 2 :R
fallback is R-specific (only kicks in when `sg.ext == 'H'`).

## Cubic b-spline boundary fix (interpolate_map)

`scipy.ndimage.map_coordinates` with `order=3, mode='wrap'` can produce overshoot at cell boundaries if adjacent grid rows have large density gradients. This creates spikes in difference maps at x=0, y=0, z=0.

All test-system maps are ASU maps (not full-cell): lysozyme covers ~½ cell in each axis; DHFR covers full X and Z but ¼ of Y; myoglobin is full-cell. For ASU maps, `mode='wrap'` wraps the boundary (e.g. x≈0.5) to x=0, which is a physically unrelated density — creating a flat artefact in the difference map near the ASU boundary.

**Fix**: `interpolate_map` pads by 5 voxels on each side with `np.pad(data, 5, mode='reflect')` before calling `map_coordinates` with `mode='nearest'`. Reflection gives a smooth continuation at any boundary (ASU edge or unit-cell edge), so the IIR prefilter sees no discontinuity. The prefilter decays as 0.268^k, so 5-voxel padding gives <0.15% boundary error at any interior query point. Do not use `pad=2` (insufficient) or `order=1` (eliminates overshoot at the cost of cubic quality everywhere).

**Additional fix for fitreso_scan with CCP4 ASU map input:** reflect-padding the *unmoved* hkl00 boundary is fine, but once the shift field is active (fr20 onwards), `delta(ref_pt)` can push the sample point past the ASU edge. The reflected value is then physically unrelated to the actual reference density at that boundary, producing a concentrated noise plane at y=0.5 / z=0.5 in `diff_norm.map` (lyso fr20 baseline: 176 voxels > 5σ, max |peak| = 14σ). The fix is to expand the moving CCP4 ASU map to the full unit cell via `read_ccp4_fullcell` and use `pad_mode='wrap'`. `fitreso_scan` now defaults `mov_fullcell=None` which auto-enables this for `.map`/`.ccp4`/`.mrc` inputs (MTZ inputs are full-cell by construction via `mtz_to_map_data`). After fix, lyso fr20 boundary peaks drop from 176 → 3 voxels > 5σ. Do not pass `mov_fullcell=False` unless you specifically want ASU-only voxels (almost always wrong for scanning).

## Map→MTZ via direct numpy FFT (replaces broken gemmi prepare_asu_data)

`_map2mtz`, `_write_scan_mtz`, and `_fspace_scale_and_diff` all FFT the
density grid with `np.fft.fftn` and look up F at the required HKLs:

```python
grid_xyz = _grid_xyz_fullcell(data, hdr)            # (NX, NY, NZ); ASU→full via gemmi setup
F_grid   = np.fft.fftn(grid_xyz.astype(np.float64))
# Crystallographic +2πi convention vs numpy's -2πi:
#   F(h,k,l) = F_grid[(-h)%NX, (-k)%NY, (-l)%NZ]
# Multiply by V/N to match gemmi's amplitude convention.
F = F_grid[(-H)%NX, (-K)%NY, (-L)%NZ] * (cell.volume / (NX*NY*NZ))
```

Phases agree with `gemmi map2sf` to 0.003° (verified by reverse-FFT of a
synthetic MTZ → numpy forward-FFT → lookup → recover injected F to 1e-7).
Amplitudes differ only by a global k_fit that the F-space LS fit absorbs.

In `_fspace_scale_and_diff`, F is looked up at the ref MTZ HKL positions
directly — 100% HKL match by construction, no `bdict` intersect that
used to drop ~half the rows.

### Why this replaced `transform_map_to_f_phi + prepare_asu_data`

`prepare_asu_data` silently dropped ~50% of unique reflections for
trigonal/hexagonal SGs in both gemmi versions we tried:
- **v0.6 (CCP4 8)** — orthorhombic-style h≥0 ∧ k≥0 ∧ l≥0 wedge for R 3:H,
  dropping all l<0.  Insulin: 22,991 rows (l ∈ [0, 30]) instead of ~46k.
- **v0.7 (CCP4 9)** — k=0 plane absent (including the entire (0,0,l)
  axis).  Insulin: 23,527 rows (k ∈ [1, 55]) instead of ~46k.

Diagnostic that exposed it:
```
gemmi map2sf --dmin 1 bent.map bent_back.mtz F PHI
diff.com bent.mtz FDM bent_back.mtz
mtzdmp Fdiff.mtz                    # → FDM 51% complete; (0,0,l) 100% missing
```

The map itself is always correct; only `prepare_asu_data` is broken.
The replacement is in commit 6f05d2c.

## SG-ASU vs Friedel-box MTZ output

All map→MTZ writers (`_map2mtz`, `_write_scan_mtz`,
`_fill_missing_with_fcalc`) emit one row per **SG-unique** HKL, not per
Friedel-unique FFT-box HKL.  Writing the same density in N point-group
copies as a P1-equivalent file amplifies the map by N on inverse-FFT
(commit d80eaa5 fix).

The enumerator `_enumerate_sg_asu_hkls(cell, sg, d_min)`:
1. Calls `generate_hkls` to get Friedel-unique HKLs at d≥d_min.
2. Maps each to its SG-canonical representative via the existing
   `get_proper_symops` + `_canonical_hkl` (lexicographic min of
   {R_k^T h} reduced to Friedel-unique form).
3. Filters out systematic absences via
   `gemmi.GroupOps.is_systematically_absent` (which is reliable,
   unlike `prepare_asu_data` / `ReciprocalAsu` for certain SGs).

Result matches CCP4 `unique` to within rounding (insulin H 3 at 1.23 Å:
25,283 vs uniqueify's 25,291; at 2.0 Å: 5,620 vs 5,620).

## Incomplete input MTZs → opt-in row-fill before refinement

Refmac writes FWT/PHWT only for HKLs present in its input MTZ.  PDB-
deposited datasets are routinely incomplete (4e7u: 93%; 4fg3: 97.6%),
so the FWT inherits the same gaps and downstream bent.mtz shows missing
chunks in Coot.

`run_refinement` now:
1. Measures SG-ASU completeness of its input MTZ
   (`_mtz_completeness` returns n_obs / n_expected — Fobs-finite rows
   only; rows with HKL+FreeR but no F don't count).
2. Refuses to proceed below 99% unless caller passes `fill_fcalc=True`
   (CLI: `--fill-fcalc`).  The name is historical; see below.
3. When set, calls `_fill_missing_with_fcalc(in, pdb, out)` which:
   - Strips any pre-existing systematic absences from input rows
     (they're zero by symmetry and confuse downstream tools).
   - Adds rows for missing SG-ASU HKLs with HKL + work-set
     `FreeR_flag` only — `FP`, `SIGFP`, `IMEAN`, `SIGIMEAN` stay NaN.
   - Original observation rows and free-R flags are untouched.

The empty-Fobs rows pass through refmac the same way it already
handles natively-empty rows in some deposits (e.g. lyso 3aw7 ships
with 683 rows that have a FreeR flag but no measured F): refmac emits
`FC`, `FC_ALL`, `FWT` and `PHWT` for them as model-only.  Output FWT/PHWT
coverage = 100% of input rows = 100% of SG-ASU after fill.  phenix.refine
likewise fills its `2FOFCWT`/`PH2FOFCWT` model coefficients across all
input rows.

**Why no longer Fcalc values?**  An earlier version injected raw
`gemmi.DensityCalculatorX` Fcalc into the FP column.  This is on the
absolute electron scale, which **does not match** the deposit's chosen
Fobs scale in general.  DHFR 1rx1/1rx2 carry FP at ~1/18 of the
absolute scale — the filled rows looked like Fobs an order of magnitude
too large, refmac's single bulk-solvent + scale model couldn't span
both populations, and rigid-body diverged (translated the model ~3.7 Å,
DHFR Rbent at fr5 climbed to 127.9% vs canonical 41.5%).  Letting
refmac/phenix generate model density themselves keeps everything on the
target scale and never introduces a second Fobs population.

`pdb_path` is still accepted by `_fill_missing_with_fcalc` for caller
stability but is unused inside.

Plumbed through `fitreso_scan(..., fill_fcalc=False)` → `run_refinement`.

`run_refinement` prefers `refmac5`; if missing it falls back to
`phenix.refine`.  The phenix call uses
`refinement.main.number_of_macro_cycles=N` — the literal parameter
name in current phenix.  Earlier code passed `refinement.main.cycles=N`,
which is not a valid parameter and **causes phenix.refine to crash**
(unknown-parameter error, not a silent no-op).

Insulin 4fg3→4e7u from raw, fill_fcalc=True:
- 4fg3: 5,487 → 5,620 (100%)
- 4e7u: 23,527 → 25,285 (100%)
- bent.mtz at every scan point: 25,285 rows = full H 3 ASU at 1.23 Å.
- Coot view: artifact-free.

## PSDVF.mtz amplitude units

The `dX`, `dY`, `dZ` columns store the **amplitude of the fitted shift
coefficient along the corresponding fractional axis, in Å** (=
sqrt(A²+B²) × cell_<axis>).  Phase (`PHX`, `PHY`, `PHZ`) is unchanged.

Internally `eval_shift_field`, `write_bent_pdb`, etc. expect AB in
fractional units.  `load_fitparams` reads the `BENDFINDER amplitude_units
angstrom` history tag (added by `save_fitparams`) and divides by
cell_<axis> on load to restore fractional AB — so the units change is
display-only.  Files written before this change have no tag and are
read as fractional unchanged.

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

## Test gamut runner (`run_all_tests.com`)

`claude/run_all_tests.com` (tcsh) runs the full 11-test gamut and prints
a PASS/FAIL table with per-test metrics:

1. `test_symm_all_sgs.py` — 65 Sohncke SGs constraint check
2. `magdoff/test_magdoff.py` — synthetic 0.5° rigid-body deformation
3–10. Eight `fitreso_scan` examples (lyso, dhfr, raddam×3, myoglobin,
   insulin, lipox) — all with `fill_fcalc=True` since the deposited
   reference MTZs are below the 99% SG-ASU completeness gate.  Lipox
   additionally exercises the cross-cell pipeline (cell stretch +
   loose-tol altindex; see "Cross-cell pairs" section above).
11. Porin altalign + refmac on the R 3 2 :R output

Each example writes to `<sys>/scan_test*/` (separate from the canonical
`scan_fitreso_fc/` reference runs).  Per-test logs land in
`test_results_<timestamp>/`.  Script cd's up if invoked from inside
`map_bender/`; exits 2 if the working area can't be found, exits 1 if
any test fails.

Run as `./run_all_tests.com` from the working area.

## Empirical results (fitreso scans)

All systems use default parameters (`outlier_sigma=2.5`, `b_sigma=3.0`, `drop_snr=0`, `batch_hkls=100`).

| System | Space group | CA pairs | fr5 RMSD | fr5 Rbent | best Rbent | d_opt | subtract |
|--------|------------|----------|----------|-----------|------------|-------|----------|
| Lyso 3aw6→3aw7 | P4₃2₁2 | 1008 | 0.033 Å | 33.2% | 30.3% | 9.7 Å | ref |
| DHFR 1rx2→1rx1 | P2₁2₁2₁ | 592 | 0.070 Å | 41.5% | 37.7% | 12.6 Å | ref |
| Raddam 5kxk→5kxl | P4₃2₁2 | 976 | 0.087 Å | 21.9% | 11.5% | 20 Å (clamped) | bent |
| Raddam 5kxk→5kxm | P4₃2₁2 | 984 | 0.047 Å | 19.2% |  9.9% | 20 Å (clamped) | bent |
| Raddam 5kxk→5kxn | P4₃2₁2 | 992 | 0.048 Å | 24.1% | 17.8% | 20 Å (clamped) | bent |
| Myoglobin 1mbo→1a6m | P2₁ | 294 | 0.063 Å | 49.8% | 49.8% | 5.5 Å | ref |
| Insulin 4fg3→4e7u | H3 | 801 | 0.510 Å | 60.3% | 60.5% | 5.8 Å | ref (fill_fcalc=True) |
| Porin 3poq→3pou | H 3 2 | 340 | 0.366 Å | 57.7% | 57.8% | 9.5 Å | ref (fill_fcalc=True) — in-bendfinder altindex resolution; obverse/reverse pair, altalign also emits H32+SYMM and R32:R (R32:R refmac-runnable, R=0.37) |
| Lipox 9o4s→9o4t | P2₁ | 795 | 0.283 Å (fr12)¹ | 53.2% (fr12) | 53.2% | 12 Å (clamped) | ref (fill_fcalc=True) — cross-cell pair (~4% expansion); stretch + loose-tol altindex picks 180°-about-z (drot=0.00° from LSQ) in nearly-orthorhombic monoclinic; refmac R=0.254 after reindex+rigid-body |

¹ Lipox numbers are fr12, not fr5 — finer scan points blow up in
wall-time at this molecule size (~6700 P1 atoms × many HKLs).  fr12
already takes ~8 minutes; fr5 would extrapolate to ~hours.  The full
fr5 row will be filled in once the canonical reference run completes.

The `best` row in each `scan_dir/scan_fitreso.log` is the
parabola-vertex re-fit (see [Best d_opt parabola fit](#best-d_opt-parabola-fit)
above).  Lyso and DHFR show genuine improvement at coarser resolution
(d_opt 9–13 Å); raddam's Rbent rises monotonically with finer
resolution so d_opt clamps to the coarsest scan point (fr20) and the
best row gains 8–10 points over fr5 — the radiation-damage signal is
purely low-frequency.  Myoglobin and insulin gain nothing from the
parabola (myoglobin's curve plateaus at fr8; insulin's T→R shift
exceeds the smooth-PSDVF model so no fitreso choice helps).

All "from-raw" runs in `<system>/scan_fitreso_fc/` (May 2026,
fill_fcalc=True propagated through refmac and altindex re-refinement).
Self-test runs: `test_symm_all_sgs.py` 65/65 PASS (max violation
2.68e-13 Å); `magdoff/test_magdoff.py` Test 2 recovery 91.8%/92.2%
(constrained/unconstrained).

Rbent values are post-F-space (k+B) scaling (compute_riso F-LS).  See per-
system README files in `lyso/`, `dhfr/`, `raddam/` for invocation details
and full peak tables.

Notes:
- **Lyso** Rbent plateaus at ~30% by fr10 and barely changes with higher
  resolution — real structural differences remain in the diff map (waters
  near the protein surface are the persistent ±5σ peaks; cell+protein
  flex are absorbed by the PSDVF).  B-factor drops from +16 (hkl00,
  pre-bending mismatch) → +8 (fr10+, true residual B).
- **DHFR** Rbent barely improves (43% → 42%); these crystal forms are
  more dissimilar.  With mov=1rx2 ref=1rx1 and the default `subtract=ref`,
  **FOL ligand atoms appear as +4–7σ positive (green) density** (O4 +6.78σ,
  C14 +9σ at low res) and the **Ca²⁺ at A/300CA(r) appears as a −10σ red
  peak** (Ca²⁺ is in 1rx1 only).  FOL is **not** in 1rx1; the title
  "complexed with" refers to NAP which is in both.  rigid-body refmac
  keeps both FWT outputs on the same scale so B stays near zero.
- **Raddam** (5kxk undamaged → 5kxl/m/n increasingly damaged) is run with
  `subtract=bent` so **positive peaks = features appearing with dose**.
  Findings across the dose series:
    * Persistent strongly-negative CYS/SG peaks (5kxk has more sulfur
      density than the damaged refs) — disulfide breakage/oxidation.
      5kxl: −10σ at A/30CYS/SG and A/115CYS/SG; 5kxn: −13 to −15σ at
      A/30CYS/SG, A/94CYS/SG.
    * **5kxm Cl⁻ accumulation**: +22.4σ positive peak at A/205CL/CL(m) by
      fr5 — chloride radiolysis product accumulating with mid-dose.
    * Waters reorganize — A/394HOH(r) −15.9σ in 5kxl, A/398HOH(r) +13.2σ
      in 5kxn.
    * fr5 RMSD: 5kxl 0.086, 5kxm 0.047, 5kxn 0.051 Å — 5kxm gives the
      cleanest fit (best data quality of the damaged set).
- **Myoglobin** fr5 Rbent 49.8%; heme iron dominates diff map throughout
  (−38.4σ at A/155HEM/FE(m) by fr5, climbing through the scan).
- **Insulin**: high Rbent (~60%) and residual RMSD (~0.5 Å) reflect
  genuine T→R conformational change — LEU B6 shifts ~8 Å between
  T-state (4fg3) and R-state (4e7u), which is outside the smooth
  shift-field model.  Dominant diff peak: −39σ at D/101ZN/ZN(r) (zinc
  position differs between T and R states).  **Requires
  `fill_fcalc=True`** — both deposited MTZs are < 99% SG-ASU complete
  (4e7u 93%, 4fg3 97.6%); without filling, refmac inherits the gaps
  and bent.mtz shows missing-HKL chunks in Coot.  Reference run:
  `scan_fitreso_fc/` (raw inputs → fill → refmac → altindex → scan).
  The hkl01 → hkl02 jump (Rbend 0.1% → 64%, the canonical (1,1,0)
  HKL picking up `|d| = 2.1 Å` in `hkl02/PSDVF.mtz`) **is not an
  over-fit** — inspecting the bent map at that stage shows the
  displacement matches a genuine low-frequency shift between the two
  crystal forms; later iterations damp the same mode to ~1 Å as
  higher-order HKLs absorb part of the signal.
- **Porin** end-to-end from raw inputs (`scan_test/`, May 2026): the
  in-bendfinder altindex resolution (the basis-change enumeration in
  [`_get_altindex_ops`](#basis-change-altindex-enumeration-_get_altindex_ops))
  picks the obverse↔reverse 2-fold and re-refines the moving model.
  hkl00 baseline is **3.14 Å CA RMSD / Rbent 73.6%** — the unbent
  model after altindex resolution is still ~3 Å from 3pou because the
  alignment also lands the crystal in the conjugate (reverse-hex)
  setting; this is benign for the scan (both moving and reference are
  on the same hexagonal grid post-resolve).  Active HKL count
  saturates at **3149 by fr10** (od_margin hit; 340 ASU CAs × 6
  proper ops = 2040 P1 atoms is the small-cell limit relative to the
  ~120 Å hex cell), so fr8/7/6/5 are identical re-fits.  fr5
  RMSD = 0.366 Å, Rbent = 57.7%; best (d_opt = 9.5 Å) Rbent = 57.8%.
  Largest residual peak **+10.5σ at A/244PHE/CB(m)** (1.73 Å) —
  side chain displacement beyond what the smooth PSDVF can follow,
  similar in spirit to insulin's T→R limit.  Total wall time
  ~72 min, dominated by fr12→fr10 (1.4 ks → 4.3 ks per fr-point) as
  HKL count grows in the H 3 2 cell — see [Test gamut runner](#test-gamut-runner-run_all_testscom)
  for tuning options (`batch_hkls`, `drop_snr`) if this matters.
- **Lipox** (soybean lipoxygenase-1 9o4s→9o4t, XFEL room-temperature
  data at 1.95 Å) is the canonical **cross-cell** test case (~4% cell
  expansion, P2₁ pseudo-orthorhombic).  Pipeline path:
  `resolve_altindex` detects cells differ → stretches mov into ref's
  cell (fractions preserved) → loose-tol (`metric_tol_rel=0.05`)
  altindex enumeration surfaces `R=diag(−1,−1,+1)` at drot=0.00°
  from the Kabsch LSQ rotation → `altindex_refine` reindexes mov's
  experimental Fobs by R, re-refines to refmac R=0.254.  Scan:
  hkl00 baseline RMSD 1.33 Å / Rbent 64.9% → fr20 0.345 Å / 53.7%
  → fr12 0.283 Å / 53.2%; parabola d_opt clamps to fr12 (true
  vertex finer than 12 Å but not yet measured).  Persistent
  −8 to −9σ peak at A/350MET/SD(r) is a sulfur not in mov.
  Mov's experimental Fobs is preserved end-to-end — no Fcalc
  substitution.  See [Cross-cell pairs](#cross-cell-pairs-resolve_altindex-cell-stretch--loose-tolerance-altindex)
  section above for the mechanism.
