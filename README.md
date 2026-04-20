# map_bender / bendfinder.com

Finds a smooth, periodic coordinate-shift field that maps one crystal form of a protein onto another. The shift field is expressed as a truncated Fourier series using real Miller-index sine functions — physically motivated basis functions that respect the periodicity of the crystal lattice.

Primary author: James Holton

## What it does

Given two PDB files of the same protein in different (non-isomorphous) crystal forms, `bendfinder.com` computes a smooth 3D vector field **Δr(x,y,z)** such that applying that field to the coordinates of `bendme.pdb` minimises the all-atom RMSD to `reference.pdb`. Optionally, a 2Fo-Fc electron density map can be re-sampled into the reference frame.

The shift field is parameterised as:

```
Δd(x,y,z) = Σ_{hkl}  a_{hkl} / N  ·  sin(2π(hx + ky + lz) + φ_{hkl})
```

where the sum runs over (h,k,l) triplets sorted by resolution (low resolution first), so the most physically meaningful large-scale deformations are captured first.

## Quick start

```tcsh
bendfinder.com  bendme.pdb  reference.pdb  [nhkls=30]
```

Both PDB files must be provided even when applying a pre-fitted function, because they define the two unit cells. The script is self-deploying: the helper programs `rmsd`, `map2pdb.com`, `floatgen`, and `origins.com` are embedded at the end and extracted on first run.

## Requirements

- **tcsh**
- **CCP4 suite** (coordconv, pdbset, mapmask, mapdump, unique, sftools, mtzdump, f2mtz, fft, refmac5)
- **gnuplot ≥ 5** (accessible as `gnuplot` or `gnuplot5`)
- **mapman** from the Uppsala Software Factory RAVE suite (only needed when a map file is provided)

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `nhkls=N` | 1000 | Maximum number of Fourier coefficients to fit |
| `starthkls=N` | 5 | Number of coefficients in the first iteration |
| `batchhkls=N` | 1 | Coefficients added per iteration |
| `reso=R` | 3 | Resolution cutoff (Å) for HKL generation |
| `frac=F` | 1 | How far along the bend path to place the output (0–1) |
| `reject=M` | 0 | Reject atom pairs > M×MAD from median shift |
| `dimensions=xyz` | xyz | Which PDB fields to fit (any subset of `x y z o B`) |
| `drop_frac=F` | 0.001 | Pre-fit: drop DFT coefficients below this fraction of the maximum |
| `drop_snr=S` | 1 | Post-fit: drop coefficients with fitted \|a\|/σ_a < S |
| `geotest=false` | true | Run refmac5 geometry check at each iteration (slow) |
| `nofit` | off | Skip gnuplot fitting; use raw DFT values directly (fast but inaccurate) |
| `deltamaps` | off | Write electron density maps of shift magnitudes Δx, Δy, Δz, Δr |
| `fitparams_N.gnuplot` | — | Apply a previously-fitted function without re-fitting |
| `refit` | — | Force re-fitting even when a fitparams file is provided |

## How it works

1. Both PDB files are expanded from their crystallographic space group to P1.
2. Atom pairs are matched by residue and atom name; fractional coordinate differences give the shift at each atom position.
3. The CCP4 `unique` program generates all (h,k,l) indices out to `reso` Å, sorted by resolution.
4. For each iteration (adding one HKL at a time), a direct Fourier transform (slow-FT) over atom positions gives starting values for the amplitude *a* and phase *φ* of each sine term.
5. **Pre-fit filtering** (`drop_frac`): DFT coefficients below `drop_frac × max(|a|)` are dropped before the gnuplot fit.
6. gnuplot 5's nonlinear Levenberg-Marquardt fitter jointly refines all active amplitudes and phases, minimising the sum of squared residuals over all atom pairs. Jobs for x, y, z run in parallel.
7. **Post-fit SNR filtering** (`drop_snr`): after fitting, coefficients with `|a_fitted| / σ_a < drop_snr` are dropped. This removes statistically insignificant terms whose DFT values were inflated by dirty-beam aliasing (non-uniform atom sampling). Dropped terms re-enter as candidates on the next iteration.
8. The fitted shift field is applied to the moving PDB and RMSD reported. The loop continues until `nhkls` coefficients have been tried.

### Why ±HKL and resolution sorting?

Earlier versions used only positive (h,k,l) and a complex-exponential basis, which requires gnuplot to handle complex arithmetic. The current approach:

- Uses **both positive and negative** indices (via CCP4 `unique`), doubling the effective basis and capturing antisymmetric deformations.
- Sorts by **resolution** (d-spacing), so the largest-scale (most physically meaningful) deformations are fitted first — analogous to fitting low-order terms before high-order ones.
- Uses **real sine functions**, which avoids complex arithmetic in gnuplot and makes the fitted parameters directly interpretable as amplitude and phase.

Benchmark on lysozyme 3aw6/3aw7 (same crystal form, two relative humidity conditions — 71.9% vs 84.2% — causing ~2.5% cell contraction):

| HKLs fitted | RMSD(CA) | Time  | Old approach equivalent        |
|-------------|----------|-------|-------------------------------|
| 5           | 0.327 Å  | 9 s   | Order 2 (19 HKLs)             |
| 20          | 0.245 Å  | 580 s | Order 3 (37 HKLs)             |
| 26          | 0.215 Å  | 1070 s| Order 4 (61 HKLs)             |
| 30          | 0.211 Å  | 1503 s| Order 5 (91 HKLs, 2938 s)    |

Roughly 3× fewer coefficients and 2× faster to reach equivalent accuracy.

## Output files

| File | Description |
|------|-------------|
| `bentN_<pdb1>.pdb` | Input PDB with shift field applied after N coefficients |
| `fitparams_N.gnuplot` | Fitted Fourier coefficients after N iterations |
| `func.gnuplot` | gnuplot function definitions for the current shift field |
| `fitparams_[xyz].gnuplot` | Per-dimension fitted parameters |
| `fitrun_[xyz].log` | gnuplot fitting output including correlation matrix |
| `param_correlations.txt` | Highly correlated parameter pairs (potential degeneracies) |
| `allhkl.txt` | Full sorted HKL list used as basis |
| `bent[N].map` | Re-sampled map at iteration N (when a map file is provided) |

## Files in this repository

| File | Description |
|------|-------------|
| `bendfinder.com` | Main script (self-deploying; embeds helper programs) |
| `origins.com` | Helper: test symmetry origin choices |
| `mapman64_notes.csh` | Notes on building mapman for 64-bit Linux |
| `mapman_regression_test.csh` | Regression test for mapman interpolation |
| `LICENSE` | License |
