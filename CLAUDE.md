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
  diff.com                      Riso calculation: scaleit-based isomorphous R-factor
  map_bender/                   git repository
    bendfinder.py
    CLAUDE.md                   this file
    README.md
    origins.com
    LICENSE
    examples/3aw6_3aw7/         canonical lysozyme example data (checked into git)
  lyso/                         lysozyme 3aw6→3aw7 (P4₃2₁2, 1008 CA pairs)
    lyso_fitreso_scan.py
    3aw6.pdb, 3aw7.pdb
    3aw6_2fofc.map, 3aw7_2fofc.map
    fitreso_scan/               hkl00..hkl10, fr5..fr20 output subdirs
  dhfr/                         DHFR 1rx1→1rx2 (P2₁2₁2₁, 636 CA pairs)
    dhfr_fitreso_scan.py
    fitreso_scan/
  myoglobin/                    myoglobin 1mbo→1a6m (P2₁, 294 CA pairs)
    myoglobin_fitreso_scan.py
    fitreso_scan/
  raddam/                       radiation damage 5kxk→5kxl/m/n (P4₃2₁2, ~980 CA pairs)
    raddam_fitreso_scan.py      takes target (5kxl|5kxm|5kxn) as argv[1]
    fitreso_scan_5kxl/, fitreso_scan_5kxm/, fitreso_scan_5kxn/
  lyso_test_031419/             gold-standard reference run (old prototype, RMSD=0.209 Å)
```

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

## Fitreso scan scripts

Each system has a `*_fitreso_scan.py` script that runs three sections and reports per-point: label, RMSD(CA), active HKLs, Rfac (Riso vs reference map), strongest positive and negative difference map peaks with nearest atom.

**Section 1 — hkl00**: zero shift; just resample the moving map onto the reference grid via `interpolate_map`. Establishes the unregistered baseline.

**Section 2 — hkl01..hkl10**: single `bend_fit_progressive` call with `batch_hkls=1, max_canon=11, fitreso_start=100`. The `iter_callback` fires once per canonical HKL added (n_non_dc = n_canon − 1). Efficient: only one initialization (atom matching, outlier rejection) for all 10 checkpoints.

**Section 3 — fr20..fr5**: separate `bend_fit_progressive` calls with `fitreso_end` in [20,15,12,10,8,7,6,5] Å.

**Riso calculation**: `diff.com` (tcsh, uses CCP4 scaleit). Map→MTZ conversion is done with `gemmi.transform_map_to_f_phi` — no sfall required, works for any space group including P2₁.

**Large-map memory**: `eval_shift_field` allocates an (N_voxels × N_hkls) phase matrix. For large maps (raddam: 3.7M voxels × 900+ HKLs ≈ 26 GB) this must be chunked. The raddam scan script uses `eval_shift_field_chunked` in 50k-voxel batches.

## Cubic b-spline boundary fix (interpolate_map)

`scipy.ndimage.map_coordinates` with `order=3, mode='wrap'` can produce overshoot at cell boundaries if adjacent grid rows have large density gradients. This creates spikes in difference maps at x=0, y=0, z=0.

**Fix**: `interpolate_map` pads by 5 voxels on each side with periodic copies (`np.pad(data, 5, mode='wrap')`) before calling `map_coordinates` with `mode='nearest'`. The IIR b-spline prefilter decays as 0.268^k, so 5-voxel padding gives <0.15% boundary error at any interior query point. Do not use `pad=2` (insufficient) or `order=1` (eliminates overshoot at the cost of cubic quality everywhere).

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

Notes:
- Lyso Rfac plateaus at ~33% by fr20 and barely changes with higher resolution — real structural differences remain in the diff map (A/74ASN/O is the persistent −10σ peak).
- DHFR Rfac barely improves (44% → 42%) — these crystal forms are more dissimilar; the FOL ligand and Mn ion dominate the diff map at all resolutions.
- Raddam Rfac *starts* low (13–11%) because the fc maps are nearly identical; huge water peaks (±30–65σ) reflect water molecules appearing/disappearing with radiation dose.
- Myoglobin Rfac pending (gemmi map2mtz re-run in progress); heme iron dominates diff map throughout (±35σ at A/154HEM/FE).
