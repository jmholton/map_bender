# sigmaA-weighted Fc fill / Fo-vs-Fc blending for map scaling

Origin: worked out with James on 2026-07-10 during the carbonic-anhydrase
water-occupancy project. Prototype + data live in a *sibling* project:
`/home/jamesh/projects/map_bender/carbonicanydrase/_wocc_batch/sigmaa_diag.py`
(and `sigmaa_diag.png`). This note ports the method here because it is the
principled version of what `_fill_missing_with_fcalc` / `--fill-fcalc` already
does, and it fixes the scale problem documented in `CLAUDE.md` (the raw
`DensityCalculatorX` Fcalc-into-FP bug, ~1/18 scale on DHFR 1rx1/1rx2).

## Problem it solves
When building a common-reference 2Fo-Fc map (esp. sharpening low-res data or
filling missing/beyond-d_min reflections), you need a per-reflection blend of
the experiment (Fo ± SIGF) and the model (Fc) that leans on Fc exactly where Fo
is noisy or absent — without introducing a second, mis-scaled Fobs population.

## The blend (phase = phi_c)
Treat true F as unknown, Fo and Fc as two estimates:

    F_best = ( Fo/so^2 + D*Fc/sD^2 ) / ( 1/so^2 + 1/sD^2 )
    Fo absent  =>  F_best = D*Fc
    w_c (weight on model) = so^2 / (so^2 + sD^2)      # ->1 when Fo noisy/absent

  so  = SIGF (experimental)
  D   = sigmaA(s) * sqrt(SigN(s)/SigP(s))     # the shrunk fill factor, D<1
  sD  = sqrt( (1 - sigmaA^2) * Sig(s) )        # residual width

James' shortcut (good zeroth order): Rfree is ~ the relative error of Fc, so
`sc ~= 1.25 * Rfree * |F|` (the 1.25 = sqrt(pi/2) converts the L1 R-value to an
L2 sigma), giving `w_c = SIGF^2 / (SIGF^2 + (1.25*Rfree*F)^2)`. The per-shell
sigmaA below is that idea done properly, resolution-resolved.

## Extrapolating into empty / beyond-d_min shells — TWO straight lines in s^2
sigmaA is not measurable where there is no Fo, so it becomes a *predicted* model
reliability. Both of these are ~linear in s^2 (Gaussian falloffs); fit on the
shells that have data, extend into the empty ones:

    ln sigmaA(s) ~= ln p - k*s^2      # Luzzati/Read; slope k = coord error, intercept p = completeness
    ln SigN(s)   (Wilson plot)        # SigN = <|F_true|^2>; SigP = <|Fc|^2> is known everywhere (you have Fc)
    D = sigmaA * sqrt(SigN/SigP)

As s grows past d_min, extrapolated sigmaA -> 0, so D*Fc -> 0: the fill FADES
OUT on its own. Self-limiting sharpening, no runaway model features. Caveat:
`ln sigmaA` linear in s^2 is OPTIMISTIC at high res (unmodeled H/solvent/disorder
steepen the real falloff) -> treat extrapolated sigmaA as an upper bound; cross-
check against sigmaA predicted from the refined coordinate error (refmac ML
SU/DPI) and take the smaller.

## TWO production findings (both cost debugging time; bake them in)
1. **Recompute Fc — do NOT use deposited FC_ALL, and do NOT use unscaled
   DensityCalculatorX.** PDB-REDO's `FC_ALL` is truncated to ~0 near d_min
   (8sf1: <Fc>=0.08 vs <Fo>=4.9 in the outer shell) — it destroys exactly the
   high-res shells you need. Recompute with bulk solvent + anisotropic scaling
   to Fo:
   ```
   # relabel FP/SIGFP -> F/SIGF first (that's what --scale-to expects), then:
   gemmi sfcalc --dmin=<d> --scale-to=<data>.mtz:F --to-mtz=fc.mtz model.pdb
   ```
   `--scale-to` optimises k_solv/B_solv + the anisotropic scale and puts Fc on
   the Fo scale — this is the fix for the FP-scale bug in CLAUDE.md. sigmaA is
   scale-free, but D/SigN/SigP need this scaling.
2. **Estimate sigmaA on a HELD-OUT free set.** Working-set sigmaA is overfit-
   optimistic (8sf1: 0.958 free vs 0.971 working). PDB-REDO MTZs use complete
   cross-validation (FREE = 0..N-1, ~5% each); FREE==0 is genuinely held out.

## Estimators (per resolution shell, on FREE==0)
    sigmaA = sum(Fo*Fc) / sqrt( sum(Fo^2) * sum(Fc^2) )   # uncentered amp corr; scale-free
    D      = sum(Fo*Fc) / sum(Fc^2)                       # least-squares slope = fill factor
    SigN/SigP = sum(Fo^2) / sum(Fc^2)
Skip d>~5 A shells when fitting the Luzzati line (bulk-solvent-dominated).

## Prototype results (2 CA cases, the plot)
- **6klz (0.90 A, high-res):** sigmaA ~= 0.98-0.99 FLAT to the edge (slope ~0).
  Model reliable everywhere -> D ~= 1, nothing to shrink/fill.
- **8sf1 (1.70 A, low-res):** clean Luzzati line, ln sigmaA = -0.003 - 0.176*s^2
  (rms coord error ~0.08 A). sigmaA 0.98 -> 0.92 across the measured range;
  extrapolated sigmaA=0.3 only at ~0.38 A.
  KEY: sigmaA (phase reliability) stays high (0.92 at d_min) but the fill factor
  D crashes to 0.29 because SigN/SigP collapses to 0.10 — the *data* is weak at
  d_min, not the phases. So D self-shrinks through the Wilson/SigN term, exactly
  as intended. Report sigmaA AND D separately; they diverge for low-res cases.

## Suggested map_bender integration
Emit `F_best, phi_c` (or a DELFWT-style column) from: recompute Fc (sfcalc
--scale-to) -> per-shell sigmaA/D on FREE==0 -> fit ln sigmaA and Wilson in s^2
-> blend with the formula above, DFc-filling missing/weak rows. This slots into
the existing `--fill-fcalc` path and is scale-safe by construction.
