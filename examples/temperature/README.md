# Temperature — cryo vs RT same crystal

| field | value |
|---|---|
| moving (`mov`)  | **4kjj** — cryogenic temperature (100 K) |
| reference (`ref`) | **4kjk** — room temperature |
| space group | P 2₁ 2₁ 2₁ |
| CA pairs | 588 |
| sign convention | `subtract=ref` (default) → diff = bent − ref |

Same protein, same crystallographic conditions, only the data-collection
temperature differs.  Most of the change is the expected thermal
contraction at cryo (~0.7 % cell change) plus correlated atomic
displacements — a smooth, well-behaved isomorphous change and a good
test of how cleanly bendfinder handles a "real" small deformation.

## Run

```sh
sh runme.sh
```

Single-core wall time: ~5–10 min.

## What `runme.sh` does

1. `wget` PDB and SF cif for 4kjj and 4kjk from RCSB
2. `gemmi cif2mtz` each SF cif (cifs contain intensities only;
   ctruncate converts to F/SIGF)
3. `ctruncate` IMEAN/SIGIMEAN → F/SIGF
4. `ccp4-python ../../bendfinder.py 4kjj.pdb 4kjk.pdb 4kjj.mtz 4kjk.mtz run_refinement scan_dir=scan_fitreso`

## Expected results

| scan point | RMSD(CA) | Rbent | Rbend | top peak |
|---|---|---|---|---|
| hkl00 | 0.361 | 40.0% | 0.0% | −15.3σ A/202CA/CA(r) |
| fr20  | 0.184 | 30.9% | 27.7% | −15.2σ A/202CA/CA(r) |
| fr10  | 0.105 | 33.2% | 32.7% | −14.5σ A/202CA/CA(r) |
| **fr5** | **0.059** | **39.4%** | 40.4% | −15.5σ A/204NAP/PN(r) |

Notes:

- **Smooth thermal contraction**: hkl00 RMSD of 0.361 Å is dominated
  by the cell change between cryo and RT.  fr10 drops to 0.105 Å — the
  shift field absorbs the bulk of the thermal expansion.
- **Persistent Ca²⁺ peak at A/202CA**: −15σ negative density throughout
  the scan — the calcium ion is in a slightly different position
  between cryo and RT, beyond what the shift field can fit.  Useful
  feature for inspecting whether a metal-binding site has subtly
  reorganised.
- **NAP at A/204** dominates at the highest resolutions (fr7–fr5) —
  NADP analog atoms slightly shift in conformation.
