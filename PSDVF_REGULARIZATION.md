# PSDVF Regularization — design notes + experiment log

Origin: 2026-07-14/15 diagnostic session on the carbonic-anhydrase pair
(mov = 8sf1 at 1.70 Å → ref = 6klz at 0.90 Å) after user (James) called
out Rbent = 73.8% at hkl00 that `diff.com` said should be ~26%.  What
started as an R-factor bug hunt turned into a broader look at how the
Fourier PSDVF basis handles the finest-fitreso regime, and how much of
what we compute is signal vs bookkeeping vs interpolation artefact.
Working data: `claude/carbonic/`, PNGs there: `wilson_ca.png` (pre-fix),
`wilson_ca_gridfix.png` (post-grid-fix), `wilson_ca_hkl00_vs_fr10.png`
(multi-scan-point overlay), `wilson_psdvf_kink.png` (PSDVF Wilson +
signal-vs-noise envelopes).

---

## 1. The grid-unification bug (implemented; the actual root cause)

### Symptom

`compute_riso` at hkl00 (zero shift, no bending applied) reported
`Rbent = 73.8%` on the CA pair while `diff.com` on raw Fobs pair
reported `Riso ≈ 26% (B = −22.9 Å²)`.  Difference maps showed atomic-
scale "holes inside peaks".

### Diagnosis

At hkl00 the bent map should be identical to the mov map, so bent FDM
should equal mov FWT plus exact zeros in the resolution gap band
(0.9–1.7 Å here).  It didn't.  Shell means:

| d (Å) | ⟨\|FDM_bent\|⟩ (pre-fix) | ⟨\|FWT_ref\|⟩ |
|---|---|---|
| 1.7–2.0 | 27 | 98 |
| 1.3–1.7 | 0.94 | 56 |
| 1.0–1.3 | 5.4 | 33 |
| 0.9–1.0 | 4.8 | 16 |

The "flat ⟨|F|⟩ ≈ 5 background across the entire 0.9–1.7 Å gap band"
is NOT PSDVF-induced (PSDVF is identically zero at hkl00 → no
convolution smearing possible) — it comes from the pipeline's
**cross-grid cubic-spline interpolation** step.

Pipeline before fix:
1. `gemmi.transform_f_phi_to_map(mov_mtz, sample_rate=1.5)` → mov
   density on **mov's own coarse grid** (spacing ≈ 0.57 Å).
2. `interpolate_map(mov_map, ref_grid_shape, ...)` — cubic b-spline
   upsample onto **ref's finer grid** (spacing ≈ 0.30 Å).
3. Forward-FFT → FDM at HKLs to ref's d_min.

Cubic b-spline's spectral response has small ripples in the passband
and a slowly-decaying tail past the source Nyquist — hence the ~5 σ
leakage in the gap band.

### Fix (already in the working tree)

Force both mov and ref MTZ loads onto the **finer** of the two
natural grids, via `gemmi.transform_f_phi_to_map(exact_size=(nu,nv,nw))`
which is a zero-padded FFT.  On the shared grid:

- Missing HKLs past the coarser MTZ's d_min are implicit zeros.
- The subsequent `interpolate_map` for the shift application is
  identity at hkl00 (each ref-grid voxel samples mov at its own
  coordinate — cubic-spline evaluation at knot points is exact).
- Forward-FFT of the bent map recovers exact HKLs where mov had
  data + exact zeros where it didn't.

Files touched:
- `bendfinder.py:mtz_to_map_data(..., exact_size=None)` — passes
  through to gemmi's zero-padded FFT.
- `bendfinder.py:fitreso_scan._load(..., exact_size=None)` — also
  handles CCP4 map inputs via FFT-then-inverse-FFT onto the target
  grid (preserves band-limit; no cubic-interp).
- `fitreso_scan` captures `ref_grid_shape = (nx, ny, nz)` after
  loading ref, then forwards it to every subsequent mov load (pre_mov
  and post-altindex mov).

### Empirical verification (post-fix)

At hkl00:

- **⟨|FDM|⟩ in the 0.9–1.7 Å gap band: 1.4×10⁻³ Å** (was 5.4 pre-fix).
- **max|FDM| in the gap: 0.42** (was ~46 pre-fix).
- **ln⟨F²⟩ in the gap band: −10 to −20** (⟨|F|⟩ ≈ 1e-5 to 1e-4).

Cross-grid leakage is gone; residual is numerical (interp kernel round-
off at fractional Δr), harmless.

### Rbent after grid fix

`hkl00 Rbent = 72.0%` (was 73.8%).  Marginal improvement because the R
factor SUM still includes 151k gap-band HKLs where FDM ≈ 0 and
FWT_ref ≈ 30 → per-HKL R contribution ≈ 1.0 → dominates the sum.

The number is a bookkeeping artefact of comparing "0" (from mov side —
correctly) with "model amplitude" (from ref side).  Not physically
meaningful.  Fixing this needs the compute_riso d_min cap (section 2
below).

---

## 2. The compute_riso + `_fspace_scale_and_diff` d_min cap (pending)

Now that FDM is honest (=0 past mov d_min), the R factor and diff-map
math should ignore that band.  Two edits:

### `compute_riso` (bendfinder.py:5482)

Add `d_min=None` kwarg.  When set, apply `d ≥ d_min` filter after
HKL-tuple match, before the positivity mask.  Backward-compatible.

### `_fspace_scale_and_diff` (bendfinder.py:4229)

Add `d_min_cap=None` kwarg.  When set, restrict k+B non-linear LS to
`d ≥ d_min_cap` HKLs; force DELFWT to zero for `d < d_min_cap` (so
its inverse-FFT is band-limited and the "holes inside peaks" go away).

### Plumbing

In `fitreso_scan`, capture `mov_dmin_orig` and `ref_dmin_orig` from
the RAW input MTZ headers (pattern already used at lines
3397/3464/3608/3843).  Compute `d_cap = max(mov_dmin_orig,
ref_dmin_orig)`.  Pass through to every `compute_riso` and
`_fspace_scale_and_diff` call site inside fitreso_scan (lines 4750,
4760, 4772, 4857, 4865).

### Expected effect (from standalone reproduction)

CA hkl00: `Rbent 72.0% → ~38%` (compute_riso restricted to d ≥ 1.7 Å).
Diff-map holes-in-peaks: gone (DELFWT truncated at mov d_min).

Same-d_min systems (lyso, dhfr, raddam, lipox): `d_cap ≈ mov ≈ ref`
so cap includes every HKL → bit-identical to today.

---

## 3. The Wilson-plot diagnostic on the PSDVF itself

### Setup

For each `fr<N>` scan point's `PSDVF.mtz` (dX/dY/dZ = fitted Fourier
coefficients per canonical HKL, in Å, POST-softpnna_kth-weighted;
SNR = raw pre-Pnn), plot `ln⟨|Δr|²⟩ = ln(dX²+dY²+dZ²)` vs `s² = 1/d²`.

Key data (CA, from `wilson_psdvf_kink.png`):

**Shell median (brown, best row):**
- Below s² ≈ 0.005 (d > 14 Å): decays roughly linearly in log-space
  (signal-dominated).
- Above s² ≈ 0.005: flattens to `ln⟨|Δr|²⟩ ≈ −14.5 to −15` (noise-
  floor-dominated).

**Signal-only fit** (SNR > 3 HKLs, n=34, no floor pollution):
`ln⟨|Δr|²⟩ = −7.49 − 212·s²` — an effective Luzzati-style envelope
with rms-coord-error analog `√(212/(2π²)) ≈ 3.3 Å`.

**Noise floor** (median of SNR < 1.5 HKLs, n=421):
`ln⟨|Δr|²⟩ ≈ −15.52 → |Δr| ≈ 0.43 mÅ`.

**Signal/noise crossover** (extrapolated): `d ≈ 5.1 Å` (past the
current scan's finest fitreso of ~6.5 Å, so we never reach it inside
the scan).

**SNR panel** (raw pre-Pnn per-HKL):
- Median SNR falls below 3 at s² ≈ 0.003 (d ≈ 18 Å).
- Below 2 at s² ≈ 0.006 (d ≈ 13 Å).
- Around 1.5 across everything beyond s² ≈ 0.01 (d < 10 Å).

### What this says

1. **The Fourier basis's low-frequency behavior is textbook-clean.**
   Multiple fitresos agree; envelope is smoothly Gaussian in s².
2. **The high-s² tail is noise-floor-dominated.**  Individual HKLs
   are scattered orders of magnitude around a tiny mean.  softPnna
   is not enough to fully suppress: hundreds of unit-SNR HKLs
   admitted at fr7/best.
3. **The kink is a Wilson-plot pathology, not a basis pathology.**
   Same shape you see on real Wilson plots when signal drops below
   noise.

### Are the high-freq outliers legit?

Three discriminators, none conclusive alone:

- **Per-HKL SNR.** Legit signal → SNR ≥ 3.  Noise excursion → SNR
  near 1.  Discriminating but scatter is huge.
- **Cross-fitreso consistency.** Real signal at (h,k,l) should show
  similar |Δr| at fr10 and fr7 (same physics, more HKLs).  Noise
  jitters randomly across fitresos.
- **Cross-validation holdout.** Real signal reduces held-out atom
  residual; over-fit noise only reduces training-atom residual.
  Infrastructure exists (`holdout_frac`, `cv_callback`) but
  fitreso_scan doesn't use it (section 4 below).

---

## 4. Cross-validation holdout status

Infrastructure already in `bend_fit_progressive` (lines 2650, 2873,
kwargs `holdout_frac`, `holdout_seed`, `cv_callback`) but NONE of
the three `bend_fit_progressive` call sites inside fitreso_scan
(lines 5095, 5232, 5365) pass it.

Turning it on is a small change: add kwargs, forward them, add a
`cv_callback` that logs `train_rmsd / hold_rmsd / hold_prefit` per
scan-row into `scan_fitreso.log`.  Default `holdout_frac=0.0`
preserves current behavior.  Section 6.

---

## 5. Basis choice — is Fourier the right basis?

Considered but tabled for later.  Alternatives:

| Basis | Locality | Reg built-in | Fit cost | Fit for our problem |
|---|---|---|---|---|
| Fourier (current) | global | fitreso cutoff only | linear | good for periodicity |
| Divergence-free Fourier (h · F(Δr) = 0) | global | volume-preserving | 2/3 params | proteins ≈ incompressible on small scales |
| Radial basis functions on atoms | local | smoothness kernel | linear + kernel | "atoms slide, don't collide" — very natural |
| Wavelets | scale-localized | multi-resolution | linear+sparse | good for localized deformations |
| Elastic normal modes | quasi-global, physics-based | mode truncation | expensive setup | in-principle correct basis |
| Voxel-wise Δr + smoothness prior | fully local | Sobolev penalty | very large | max flexibility, needs strong reg |

**Assessment from the PSDVF Wilson diagnostic:** the current basis's
low-frequency behavior is fine.  The over-fit lives in the high-s²
tail where each HKL is a unit-SNR fit.  A per-HKL regularizer on top
of the current basis is likely enough — no need for a full basis
change as first move.  Local basis (RBFs) is worth prototyping if
per-HKL regularization is insufficient.

---

## 6. Regularization strategies — the menu we're trying

Ordered by cheapness / target-precision:

### 6a. Cross-validation holdout as a validation metric (not a regularizer per se)

Turn on `holdout_frac` in fitreso_scan.  Logs `train / hold` RMSD
per iteration.  Doesn't change the fit; provides ground truth on
whether admitting each new HKL helps unseen atoms or just training
atoms.  Prerequisite for evaluating regularizers below.

### 6b. Sobolev smoothness penalty (Fourier-space)

Modify `fit_lstsq_symm`'s Tikhonov ridge from constant λ to
`λ(h) = λ₀ · (1 + β · |h|²)` — high-|h| coefficients incur more
penalty.  Composes with existing `bound_by_obs` ridge.  Kwargs:
`sobolev_lambda`, `sobolev_beta`.  Cheap: one-line change inside
the SVD filter.  Targets the specific over-fit mode (high-|h|
tail).

### 6c. Wilson-envelope post-hoc cap

After SVD, fit envelope on SNR>3 HKLs; establish noise floor from
SNR<1.5 HKLs.  For each HKL: downweight if |Δr|² > α·envelope(s²)
AND SNR < threshold (so we don't kill legit high-SNR outliers).
Applied to the (A, B) block in `save_fitparams` and to the
returned PSDVF coefficients used to build the bent map.  Kwargs:
`wilson_cap_alpha`, `wilson_cap_snr_min`.  Moderately cheap.

### 6d. Divergence-free constraint

Project each canonical HKL's (A_x, A_y, A_z, B_x, B_y, B_z) block
onto the subspace perpendicular to h (both real and imaginary
parts).  Reduces free parameters by 1/3, forces volume-preserving
shift field.  Cheap: per-HKL projection in `fit_lstsq_symm`.
Kwarg: `divfree=True`.  Physics-motivated.

### 6e. Local basis (RBFs on CA positions)

Bigger design change.  Deferred until 6a–6d have been evaluated.

---

## 7. Experiment plan

Test on two systems:

- **CA (8sf1 → 6klz)**: cross-d_min pair, has the over-fit tail
  we're trying to suppress.  Diagnostic: gap-band FDM after
  regularization + Rbent (with cap) + held-out atom RMSD.
- **lyso (3aw6 → 3aw7)**: matched d_min, already clean baseline.
  Diagnostic: regularization must NOT hurt lyso's numbers.
  Regressions here mean the regularizer is too aggressive.

Baseline (post-grid-fix, post-cap): scan both systems.  Then run
under each configuration:

- `+sobolev` — sobolev_lambda tuned on lyso to be a no-op there,
  then apply to CA.
- `+wilson` — wilson_cap_alpha=5, snr_min=2.
- `+divfree` — enforce div-free.
- `+cv` — holdout enabled; check hold vs train ratio.

Compare per system: RMSD, Rbent, bondZ, gap-band FDM ⟨|F|⟩, hold
vs train ratio.  Winning combination is one that leaves lyso
unchanged AND reduces CA's over-fit tail.

---

## 8. Status

- [x] Section 1 (grid unification): implemented, verified on CA.
- [ ] Section 2 (d_min cap): pending.
- [x] Section 3 (Wilson diagnostic): observed, documented above.
- [ ] Section 4 (cv holdout wiring): pending.
- [ ] Section 6 (regularizers): pending.
- [ ] Section 7 (experiments): pending.
