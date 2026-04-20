# map_bender / bendfinder.com — Claude context

## Project overview

`bendfinder.com` is a tcsh script that fits a smooth, periodic 3D coordinate-shift field between two non-isomorphous crystal forms of the same protein. The shift field is a Fourier series of real sine functions indexed by Miller indices (h,k,l).

The git repo lives at `./map_bender/` (relative to this directory). The working development copy is `./bendfinder.com`. Keep both in sync when making changes — both files are identical and should stay that way.

## Directory structure

```
./                          working area, not in git
  bendfinder.com            current development version (keep in sync with map_bender/)
  map_bender/               git repository
    bendfinder.com          version in git
    README.md
    origins.com
    LICENSE
  lyso_test_031419/         gold-standard reference run
    bendfinder.log          reference results: order 0=0.771, 1=0.413, 2=0.327, 3=0.246, 4=0.215, 5=0.209
    bendfinder.com          script version used for that gold-standard run
    3aw7_refine_001.pdb     moving PDB (lysozyme form 2)
    3aw6_refine_001.pdb     reference PDB (lysozyme form 1)
    3aw7_2fofc.map          2Fo-Fc map for re-sampling
  test1/                    oldest reference run (positive-only HKLs, complex exponential basis)
    bendfinder.com          original script (2017)
    3aw6.pdb, 3aw7.pdb
  test_new/                 sandbox for new runs
  test_goldstd/             benchmark run comparing new vs gold-standard
```

## Algorithm summary

1. Both PDBs expanded from their space group to P1.
2. Atom pairs matched by residue+atom name; fractional Δr = shift at each atom.
3. CCP4 `unique` generates all (h,k,l) to `reso` Å resolution, sorted by d-spacing (low-res first).
4. **Big loop** over nhkls (starthkls to maxhkls, stepping batchhkls):
   - **Slow-FT**: compute DFT over atom positions as starting values for new HKL.
   - **Pre-fit drop** (`drop_frac`): drop DFT amplitudes < drop_frac × max amplitude.
   - **gnuplot5 fit**: jointly refine all active (a, φ) pairs; x/y/z run in parallel.
   - **Post-fit SNR drop** (`drop_snr`): drop fitted params with |a|/σ_a < drop_snr. These get a second chance next iteration. This is the "dirty-beam deconvolution" step.
   - Rebuild fitparams.gnuplot, apply shift field, report RMSD(CA).

## Key parameters

| Param | Default | Role |
|-------|---------|------|
| `nhkls` | 1000 | Max Fourier terms |
| `starthkls` | 5 | Start here |
| `batchhkls` | 1 | Add one per iteration |
| `reso` | 3 | Å resolution cutoff for HKL list |
| `drop_frac` | 0.001 | Pre-fit: fraction of max DFT amplitude |
| `drop_snr` | 1 | Post-fit: SNR threshold (|a|/σ_a) |
| `fitscale` | 1000 | Internal amplitude scale to avoid gnuplot numerical issues |
| `dimensions` | x y z | Which coordinate shifts to fit |
| `geotest` | true | Run refmac5 geometry check (slow — set false for testing) |
| `nofit` | unset | Skip gnuplot; use raw DFT (fast, inaccurate due to dirty-beam) |

## Variable name encoding

Negative HKL indices are encoded with `m` and positive with `p` in gnuplot variable names (gnuplot doesn't allow `-` in variable names). The `0` case has no prefix.

Examples:
- h=0, k=1, l=0  →  `a_dx0_1_0`
- h=1, k=-1, l=0 →  `a_dx1_m1_0`
- h=-1, k=0, l=1 →  `a_dxm1_0_1`

## Important files generated during a run

- `allhkl.txt` — full sorted HKL list (resolution-ordered); `0 0 0` always first
- `ehkl.txt` — encoded HKL names for current iteration
- `fitparams.gnuplot` — all fitted (a, φ, a_err, φ_err) for all active HKLs
- `fitparams_{x,y,z}.gnuplot` — per-dimension split (fed to parallel gnuplot jobs)
- `func.gnuplot` — gnuplot function definitions `dx(x,y,z)=...`
- `fitme` — flat table of atom fractional coords + shifts (input to gnuplot)
- `fitrun_{x,y,z}.log` — gnuplot output including correlation matrix
- `param_correlations.txt` — pairs with |CC| > 0.90 (watch for degeneracy)

## fitparams format

gnuplot's `update` command appends fitted values, so each parameter has **two** `_err` lines — the first is the initialization (1e-10), the second is the actual fitted uncertainty. When parsing, take the **last** occurrence:

```
a_dx0_1_0 = 7.03
a_dx0_1_0_err = 1e-10        ← initialization, ignore
a_dx0_1_0_err = 0.058        ← fitted σ, use this
phi_dx0_1_0 = 0.0023
phi_dx0_1_0_err = 1e-10
phi_dx0_1_0_err = 0.00145
```

## Post-fit SNR pruning (drop_snr)

Added at lines ~796–822 (after gnuplot fit, before correlation analysis). Uses awk to:
1. Store all lines in memory.
2. For each `a_d*` param: record value; last `a_d*_err` line is the fitted σ.
3. Compute SNR = |a| / σ. If SNR < drop_snr, mark as bad.
4. Drop bad params and their associated `phi_d*` lines from fitparams.gnuplot.
5. Dropped params re-enter as "new" on next iteration → recomputed from DFT → likely dropped again by drop_frac if genuinely zero. This prevents dirty-beam ghost terms from accumulating.

## Why the old approach was slower

The old (2017) version used:
- Positive-only h,k,l (missed antisymmetric deformations)
- Complex exponential form `real(a·exp(2πi(hx+ky+lz+φ)))` (required gnuplot complex support)
- HKL enumeration by `max(h,k,l) ≤ order` (not resolution-sorted)

The new approach gets equivalent RMSD with ~3× fewer HKLs and ~2× faster, primarily because resolution sorting prioritises the physically significant large-scale components.

## Gold-standard benchmark (lyso 3aw6/3aw7)

Run from `test_goldstd/` using `map_bender/bendfinder.com` with `nhkls=30 geotest=false`:

| nhkls | RMSD(CA) | Time |
|-------|----------|------|
| 5     | 0.327 Å  | 9 s  |
| 9     | 0.277 Å  | 34 s |
| 20    | 0.245 Å  | 580 s|
| 26    | 0.215 Å  | 1070 s |
| 30    | 0.211 Å  | 1503 s |

Gold-standard (old approach, order 5, 91 HKLs): RMSD=0.209 in 2938 s.

## Helper programs (self-deployed)

`bendfinder.com` is self-deploying. The programs `rmsd`, `map2pdb.com`, `floatgen.c`, and `origins.com` are embedded as shell-here-doc payloads at the end of the script and extracted into `$PATH` on first run (triggered by `deploy_scripts` label).

## Common pitfalls

- `gnuplot5` must be version ≥ 5 with `set fit errorvariables` support. The script tries `gnuplot` first and aliases if it's v5.
- `CCP4_SCR` must be set (CCP4 environment). The script aborts if not.
- `drop_snr=0` disables post-fit SNR pruning entirely.
- The `nofit` option produces poor RMSD (dirty-beam aliasing) and is only useful for rapid testing of the pipeline or starting-value inspection.
- Fitting time scales roughly O(N²) in the number of active parameters per dimension. SNR pruning keeps this from blowing up at large nhkls.
