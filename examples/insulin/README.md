# Insulin — T → R quaternary transition

| field | value |
|---|---|
| moving (`mov`)  | **4fg3** — insulin hexamer, T-state |
| reference (`ref`) | **4e7u** — insulin hexamer, R-state |
| space group | H 3 |
| CA pairs | 99 (ASU; ~800 in P1) |
| sign convention | `subtract=ref` (default) → diff = bent − ref |

Hexameric insulin in the trigonal lattice in two quaternary states.
The two PDBs deposit different ASU chain assignments — when matched
naively by `(chain, residue, atom_name)`, the implied rigid-body
rotation between the structures is essentially one of H 3's own
symmetry operators (a 3⁻ rotation about c, ~120°) combined with a
non-standard origin shift.  `resolve_altindex` discovers and applies
this combined op (re-refining the moving structure under it) before
the bend.  After that registration, the *real* residual is the T→R
conformational change — LEU B6 swings ~8 Å, the B-chain helix
re-packs, and the Zn coordination differs.  The bend fit absorbs most
of this residual; what's left is the largest difference feature in
the suite.

## Run

```sh
sh runme.sh
```

Single-core wall time: ~5 min.

## What `runme.sh` does

1. `wget` PDB and SF cif for 4fg3 and 4e7u from RCSB
2. `gemmi cif2mtz` each SF cif
3. `ccp4-python ../../bendfinder.py 4fg3.pdb 4e7u.pdb 4fg3.mtz 4e7u.mtz run_refinement scan_dir=scan_fitreso`

## Expected results

| scan point | RMSD(CA) | Rbent | Rbend | top peak |
|---|---|---|---|---|
| hkl00 | 2.612 | 98.1% | 0.0% | −39.3σ D/101ZN/ZN(r) |
| fr20  | 1.635 | 78.2% | 132.9% | −34.6σ D/101ZN/ZN(r) |
| fr10  | 0.864 | 72.8% | 138.8% | −33.3σ D/101ZN/ZN(r) |
| **fr5** | **0.907** | **73.0%** | 140.2% | −33.1σ D/101ZN/ZN(r) |

Notes:

- **What `resolve_altindex` picks** — `R_frac = [[-1,1,0],[-1,0,0],[0,0,1]]`
  with `t_frac = (0, 0, 1/3) mod 1`.  This isn't a *genuine* alt-indexing:
  the rotation is the H 3 three-fold (already an SG operator), combined
  with a non-standard origin shift that maps one of 4fg3's H 3-equivalent
  copies onto 4e7u's atomic positions.  Equivalently, the LSQ rigid-body
  rotation between the two ASUs comes out at 119.88° — within 0.12° of a
  perfect 120° — i.e. the moving and reference deposits chose different
  H 3-equivalent ASUs.  bendfinder's enumeration finds the right combined
  op anyway and re-refines.
- **Residual after registration is real**: hkl00 RMSD = 2.6 Å reflects
  the actual T → R conformational change (~8 Å LEU B6 displacement,
  B-chain re-packing).  Bending reduces this to 0.9 Å — a 65 %
  reduction — which is about what a smooth low-order shift field can
  do for a genuinely large local rearrangement.
- **Zn at D/101 ZN**: the strongest persistent peak (−33σ across all
  resolutions) is the zinc cofactor at the hexamer-axis site, which
  occupies a different position between T and R.  The shift field
  can't follow this discrete site change.
- **Rbend > 100 %**: the bending changes the bent map more than the
  bent map differs from ref — a clear sign that the shift field is
  pushing density around to fit features it cannot fully represent
  (B-chain helix re-packing being the biggest contributor).
- This example illustrates the limit of bendfinder's smooth deformation
  model: it handles the per-monomer rigid-body re-registration cleanly
  (via the altindex+origin step), but the per-residue T→R rearrangement
  is too discrete to fully absorb.  Useful for understanding when to
  reach for a different tool (e.g. molecular dynamics or rebuilding).

## See also

- `../lyso/`, `../dhfr/`, `../temperature/` — smaller smooth changes
  that bendfinder absorbs almost completely.
