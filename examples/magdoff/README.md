# Magdoff-Crick synthetic deformation validation

| field | value |
|---|---|
| structure | **7rsa** (ribonuclease A) |
| space group | P 2₁ |
| 248 CA in P1 | a=30.2 Å, b=38.4 Å, c=53.3 Å, β=105.9° |

Two synthetic deformations of the same model — a ground-truth sanity
check that the bend fit recovers a known displacement.  Unlike the
other examples, **there is no second crystal form**; the moving PDB is
just 7rsa.pdb with a known transformation imposed on it.

## Background — the Magdoff & Crick papers

The synthetic test is named for two classic papers, both using
ribonuclease as the case study:

> **Magdoff, B.S. & Crick, F.H.C. (1955)** *Ribonuclease. II. Accuracy
> of measurement and shrinkage.* Acta Cryst. **8**, 461.
>
> **Crick, F.H.C. & Magdoff, B.S. (1956)** *The theory of the method of
> isomorphous replacement for protein crystals. I.* Acta Cryst. **9**,
> 901.

The 1955 paper documents a continuous, reversible sub-Ångström
cell-dimension fluctuation in ribonuclease II (P 2₁) crystals, driven
by small temperature gradients across the mounting capillary and the
associated humidity changes.  The point relevant here: such cell
fluctuations produce measurable intensity changes — enough to confound
isomorphous-replacement measurements unless they're either suppressed
(mount wetter, equalise temperature) or modeled out.

The 1956 paper then develops the first-order theory.  Working out the
expected intensity change in centric and acentric reflections for both
heavy-atom addition and for small displacements of the protein
molecules, it derives explicit linear-in-displacement formulae for
translations, rotations, cell-parameter changes, and "breathing"
motions.  Two key qualitative results: (i) ΔI grows linearly with 1/d,
so the same imposed displacement matters more at high resolution, and
(ii) very small molecular shifts already interfere with isomorphous
replacement at high 1/d while remaining negligible at low 1/d.

Both results are exactly the regime a low-order Fourier shift field
operates in: ΔF is approximately linear in the imposed displacement
field, and the bandlimited cutoff is a feature, not a bug — it
restricts the fit to the resolution range where first-order theory is
valid.  bendfinder solves the inverse problem: given the difference in
atom positions between two crystal forms, recover the shift field.
These tests use **7rsa (ribonuclease A, P 2₁)** in honour of the
original analysis, and Test 1 (uniform cell scaling 0.5 %) is a direct
synthetic analogue of the Magdoff-Crick shrinkage.

## What `make_test_pdbs.py` produces

- **`7rsa_scaled.pdb`** — `CRYST1` `a`, `b`, `c` multiplied by 1.005
  (0.5 % expansion), atom Cartesian coordinates unchanged.  Same atoms,
  slightly different cell → fractional coords shift by ≈ −0.005 · r_cart.
  This is the pure Magdoff-Crick scenario: an isomorphous unit-cell
  expansion with no intramolecular rearrangement.  The `SCALE`/`ORIGX`
  cards are stripped so gemmi re-derives the fractional matrix from the
  modified cell instead of the stale matrix recorded on those cards.
- **`7rsa_rotated.pdb`** — atoms rotated 0.5° in Cartesian space about a
  random axis (RNG seed 42).  Not SG-consistent for a general axis
  (the 2₁ screw mixes the rotation), so an SG-symmetric shift field can
  only fit a projection of the true displacement; the unconstrained fit
  recovers a little more but uses fewer effective HKLs.

## Run

```sh
sh runme.sh
```

Single-core wall time: < 1 min total.

## What `runme.sh` does

1. `wget` 7rsa.pdb from RCSB
2. `ccp4-python make_test_pdbs.py` → `7rsa_scaled.pdb`, `7rsa_rotated.pdb`
3. `ccp4-python ../../bendfinder.py 7rsa_scaled.pdb  7rsa.pdb fitreso=7`
4. `ccp4-python ../../bendfinder.py 7rsa_rotated.pdb 7rsa.pdb fitreso=7`

Each `bendfinder.py` invocation runs the constrained (P 2₁) bend fit at
the 7 Å resolution endpoint that Magdoff-Crick's first-order theory
applies to.  Each writes a `bentNN.pdb` (the moving model with the
recovered shift field applied) and a `fitparams.npz`.

## Expected results

bendfinder.py output for each fit (snip):

| Test | Deformation | Imposed RMSD(CA) | RMSD after bend_fit | recovery |
|---|---|---|---|---|
| 1 | cell ×1.005 | 0.197 Å | **0.011 Å** | **94 %** |
| 2 | 0.5° rotation | 0.335 Å | **0.017 Å** | **95 %** |

The reference `runme.log` is checked in alongside the README.

Why the residual is not zero: the shift field is bandlimited
(`fitreso=7` Å), so it can fit the imposed deformation only down to its
spatial-resolution cutoff.  The residual (1–2 % of the imposed RMSD) is
high-frequency content that lies beyond the PSDVF bandwidth — not a fit
failure within its design range.

If you compute Riso from refmac-recomputed Fcalc at 1.5 Å data
resolution (the data resolution is independent of the PSDVF
bandwidth), you should see ~15 % → 3 % for test 1 and ~24 % → 4 % for
test 2.

## See also

- `examples/lyso/` — real (not synthetic) Magdoff-Crick test on lysozyme
  at two humidity points (3aw6 → 3aw7).
- `../../CLAUDE.md` § "Magdoff synthetic deformation tests" — extended
  notes on both tests.
