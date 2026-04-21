# map_bender / bendfinder.com — Claude context

## Project overview

`bendfinder.com` is a prototype script that fits a smooth, periodic 3D coordinate-shift field between two non-isomorphous crystal forms of the same protein. The shift field is a Fourier series indexed by Miller indices (h,k,l).

The git repo lives at `./map_bender/` (relative to `../`). The working development copy is `../bendfinder.com`. Keep both in sync when making changes — they must stay identical.

## Directory structure

```
../                             working area, not in git
  bendfinder.com                current development version (keep in sync with map_bender/)
  map_bender/                   git repository
    bendfinder.com
    CLAUDE.md                   this file
    README.md
    origins.com
    LICENSE
    examples/3aw6_3aw7/         canonical lysozyme example (checked into git)
  lyso_test_031419/             gold-standard reference run (old prototype version)
    bendfinder.log              order 0=0.771, 1=0.413, 2=0.327, 3=0.246, 4=0.215, 5=0.209
    bendfinder.com              original script used for gold-standard
    3aw7_refine_001.pdb         moving PDB (lysozyme form 2)
    3aw6_refine_001.pdb         reference PDB (lysozyme form 1)
    3aw7_2fofc.map              2Fo-Fc map for re-sampling
  test_temperature/             4kjk / 4kjj — same crystal, two temperatures, P212121
    run4/                       nhkls=100 batchhkls=5: RMSD≈0.28 Å
    run5/                       starthkls=100 single-shot: RMSD=0.316 Å (worse — see nofit)
    nofit_N/                    single-shot nofit at N HKLs (diverges above N~30)
  myoglobin/                    1mbo → 1a6m, P2₁, 151 CA pairs
    run1/                       nhkls=100: RMSD=0.229 Å in 237 s
  raddam/                       radiation damage series vs 5kxk (P4₃2₁2, 169 CA)
    run_5kxl/                   RMSD=0.177 Å in 643 s
    run_5kxm/                   RMSD=0.165 Å in 506 s  (smallest shift)
    run_5kxn/                   RMSD=0.220 Å in 558 s
  scan_3aw6_3aw7/               parameter scan results (geotest=false, nhkls=30)
    batch1…batch5/              batchhkls scan; batch5 is fastest at same RMSD
    snr2/ snr3/                 drop_snr scan
    nhkls20/                    nhkls=20 reference
    nofit_lyso_N/               single-shot nofit divergence tests
  test1/                        oldest reference (positive-only HKLs, complex exponential)
  test_new/                     sandbox
  test_goldstd/                 benchmark run (new vs gold-standard)
```

## Algorithm summary (prototype version)

1. Both PDBs expanded from their space group to P1 (via CCP4 `pdbset`).
2. Atom pairs matched by residue+atom name; fractional Δr = shift at each atom.
3. CCP4 `unique` generates all (h,k,l) to `reso` Å, sorted by d-spacing (low-res first).
4. **Big loop** over nhkls (starthkls → maxhkls, stepping batchhkls):
   - **Slow-FT**: compute DFT over atom positions **from current residuals** as starting values for new HKLs.
   - **Pre-fit drop** (`drop_frac`): drop DFT amplitudes < drop_frac × max amplitude.
   - **gnuplot5 fit**: jointly refine all active (a, φ) pairs; x/y/z run in parallel.
   - **Post-fit SNR drop** (`drop_snr`): drop fitted params with |a|/σ_a < drop_snr. Dropped params re-enter next iteration (dirty-beam deconvolution).
   - Rebuild fitparams.gnuplot, apply shift field, report RMSD(CA).

**Critical**: the slow-FT operates on the *residual* (observed minus current model). Single-shot nofit (all HKLs at once without a prior model) diverges — see below.

## Key parameters

| Param | Default | Role |
|-------|---------|------|
| `nhkls` | 1000 | Max Fourier terms |
| `starthkls` | 5 | Start here |
| `batchhkls` | 1 | Add this many per iteration |
| `reso` | 3 | Å resolution cutoff for HKL list |
| `drop_frac` | 0.001 | Pre-fit: fraction of max DFT amplitude |
| `drop_snr` | 1 | Post-fit: SNR threshold (|a|/σ_a) |
| `fitscale` | 1000 | Internal amplitude scale (avoids gnuplot underflow) |
| `dimensions` | x y z | Which coordinate shifts to fit |
| `geotest` | true | Run refmac5 geometry check (slow — set false for testing) |
| `nofit` | unset | Skip gnuplot; use raw DFT values (see nofit pitfall below) |

## Parameter scan results (3aw6/3aw7, nhkls=30, geotest=false)

| Run | Change | RMSD(CA) | Time |
|-----|--------|----------|------|
| batch1 (baseline) | batchhkls=1 | 0.211 Å | 1503 s |
| **batch5** | **batchhkls=5** | **0.210 Å** | **301 s** |
| batch3 | batchhkls=3 | 0.217 Å | 331 s |
| batch2 | batchhkls=2 | 0.217 Å | 495 s |
| snr3 | drop_snr=3 | 0.211 Å | 710 s |
| snr2 | drop_snr=2 | 0.211 Å | 719 s |
| nhkls20 | nhkls=20 | 0.245 Å | 247 s |

**Best: `batchhkls=5`** — same RMSD as default, 5× faster. Higher drop_snr doesn't help at nhkls=30. nhkls=20 saves ~50 s at cost of 0.035 Å.

Gold-standard (old approach, order 5, 91 HKLs): RMSD=0.209 Å in 2938 s.

## nofit single-shot divergence

`nofit` with `starthkls=nhkls` (all HKLs in one shot) **diverges** because each HKL computes its DFT amplitude from the full unmodeled shift, not the residual. Summing N such terms produces a divergent dirty-beam superposition.

| nhkls | lyso nofit | temp nofit |
|-------|-----------|-----------|
| 5     | 0.497 Å   | 0.607 Å   |
| 30    | **31.2 Å** | 0.665 Å  |
| 300   | 31.2 Å    | **21.9 Å** |
| 1000  | 31.2 Å    | 18.9 Å    |

The incremental loop with nofit is a greedy matching pursuit (CLEAN algorithm) that works by projecting each new HKL onto the current residual. It works but accumulates dirty-beam errors that the gnuplot fit step corrects.

A linear lstsq solver (numpy) eliminates this entirely — fit all HKLs at once, globally optimal, no loop needed.

## Variable name encoding (gnuplot fitparams)

Negative HKL indices encoded with `m`, positive with `p`, zero with nothing:

- h=0, k=1, l=0  →  `a_dx0_1_0`
- h=1, k=-1, l=0 →  `a_dx1_m1_0`
- h=-1, k=0, l=1 →  `a_dxm1_0_1`

## Important runtime files

- `allhkl.txt` — full sorted HKL list; `0 0 0` always first
- `ehkl.txt` — encoded HKL names for current iteration
- `fitparams.gnuplot` — all fitted (a, φ, a_err, φ_err) for active HKLs
- `fitparams_{x,y,z}.gnuplot` — per-dimension split (parallel gnuplot jobs)
- `func.gnuplot` — gnuplot function definitions `dx(x,y,z)=...`
- `fitme` — atom fractional coords + shifts (input to gnuplot)
- `fitrun_{x,y,z}.log` — gnuplot output including correlation matrix
- `param_correlations.txt` — pairs with |CC| > 0.90 (watch for degeneracy)

## fitparams format

gnuplot `update` appends fitted values, so each parameter has two `_err` lines — first is initialization (1e-10), second is actual σ. Always take the **last** occurrence:

```
a_dx0_1_0 = 7.03
a_dx0_1_0_err = 1e-10        ← initialization, ignore
a_dx0_1_0_err = 0.058        ← fitted σ, use this
```

## Post-fit SNR pruning (drop_snr)

Lines ~796–822 (after gnuplot fit, before correlation analysis). awk reads fitparams.gnuplot, computes SNR = |a| / σ_a, drops params below threshold. Output routing: informational `SNR-pruned N params` line is mixed in tempfile; `egrep ^SNR-pruned` displays it, `egrep -v ^SNR-pruned` strips it before writing fitparams.gnuplot (prevents gnuplot parse errors).

## tcsh features to keep in mind

1. **Multi-line awk in single quotes**: tcsh requires `\` at the end of every internal line of a `'...'` string. Bare newlines cause `Unmatched '.` at runtime.
2. **`!` in awk patterns**: tcsh expands `!` for history substitution before quote parsing. `!/pattern/` triggers `Event not found`. Use `index($field,"str")==0` instead.
3. **`next` in awk END block**: invalid in awk. Use `continue` inside a for loop in END.
4. **PDB args are positional**: the script matches `*.pdb` glob in the arg loop. Pass PDB files as bare paths (`file1.pdb file2.pdb`), not `pdb1=file.pdb` — the literal string `pdb1=file.pdb` fails the `-e` existence check.

## Helper programs (self-deployed)

`bendfinder.com` is self-deploying. `rmsd`, `map2pdb.com`, `floatgen.c`, and `origins.com` are embedded as here-doc payloads, extracted into `$PATH` on first run (`deploy_scripts` label).

## Common pitfalls

- `gnuplot5` must be version ≥ 5 with `set fit errorvariables`. Script auto-detects.
- `CCP4_SCR` must be set. Script aborts if not.
- `drop_snr=0` disables SNR pruning entirely.
- Fitting time scales O(N²) in active parameters/dimension. batchhkls=5 keeps gnuplot comfortable; above ~100 simultaneous params, L-M fails to converge.
- On SLURM login nodes: `srun --cpus-per-task=3` (3 for parallel x/y/z gnuplot jobs).

---

## Python / C port plan

### Motivation

The gnuplot (a, φ) nonlinear parameterization forces the incremental batchhkls loop because L-M fails above ~100 simultaneous parameters. Switching to a **linear (A, B) parameterization** eliminates this entirely:

```
dx = Σ_hkl [ A_hkl · sin(2π(hx+ky+lz)) + B_hkl · cos(2π(hx+ky+lz)) ]
```

Build design matrix X (atoms × 2·nhkls), solve `X @ [A,B] = shifts` with lstsq. Recover `a = sqrt(A²+B²)`, `φ = atan2(B,A)`. Uncertainties from covariance diagonal `σ²·(XᵀX)⁻¹`. All HKLs fit simultaneously, no incremental loop, guaranteed convergence.

### Architecture

```
bendfinder.py               main script (replaces bendfinder.com)
  pdb_io.py                 read PDB, match atoms, expand to P1 (via gemmi)
  hkl.py                    enumerate HKLs to resolution cutoff (replaces CCP4 unique)
  fit.py                    design matrix, lstsq, SNR pruning
  apply.py                  evaluate shift field, write output PDB

design_matrix.c (optional)  fast inner loop if numpy is too slow
```

### Phase 1 — Python core

- **pdb_io.py**: read PDB, match atom pairs by resname+resnum+atomname, compute fractional Δr. Use `gemmi` for P1 expansion (replaces `pdbset`/`pdbcur`).
- **hkl.py**: enumerate all (h,k,l) with d ≥ reso Å, sort by d descending, deduplicate Friedel pairs. Pure numpy — replaces CCP4 `unique`.
- **fit.py**: build X; call `scipy.linalg.lstsq`; recover (a, φ, σ_a); SNR-prune (drop HKLs where |a|/σ_a < drop_snr); re-solve; repeat until stable (1–2 rounds typically).
- **apply.py**: evaluate shift field at atom positions; write bent PDB.
- Keep `refmac5` geotest and `rmsd` helper as external calls.

### Phase 2 — C extension (if needed)

Profile first. Likely bottleneck is building X: O(N_atoms × N_hkls) trig calls. If slow:

```c
// design_matrix.c — called via ctypes
void build_design_matrix(
    const double *atoms,  // Nx3 fractional coords
    const int    *hkls,   // Mx3 Miller indices
    int N, int M,
    double *X);           // Nx2M output
```

LAPACK is already used internally by scipy for lstsq; no need to call it from C.

### Phase 3 — Reduce CCP4 dependency

- HKL generation: done in Phase 1 (pure Python/numpy).
- P1 expansion: done in Phase 1 (gemmi).
- Map resampling (deltamaps): keep CCP4 FFT for now; replace with `gemmi` later.
- Geometry check (`geotest`): keep `refmac5` — no suitable free alternative.
- `rmsd` helper: keep C binary or rewrite in Python with gemmi.

### Comparison

| Feature | prototype + gnuplot | Python + lstsq |
|---------|---------------|----------------|
| Parameterization | (a, φ) nonlinear | (A, B) linear |
| Incremental loop | Required (load-bearing) | Not needed |
| Convergence | L-M, can fail >~100 params | Guaranteed |
| Uncertainties | gnuplot errorvariables | Covariance matrix |
| nhkls limit | ~30–50/dimension | 1000+ trivially |
| gnuplot dependency | Required | Eliminated |
| CCP4 dependency | pdbset + unique + refmac5 | refmac5 only (Phase 3: none) |
