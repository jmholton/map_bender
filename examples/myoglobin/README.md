# Myoglobin — heme iron and packing differences

| field | value |
|---|---|
| moving (`mov`)  | **1mbo** — sperm-whale Mb, oxy form |
| reference (`ref`) | **1a6m** — sperm-whale Mb, met form (re-refined) |
| space group | P 2₁ |
| CA pairs | 294 |
| sign convention | `subtract=ref` (default) → diff = bent − ref |

Two crystal forms of sperm-whale myoglobin.  Both have the heme group
(HEM) and Fe centre, but the Fe oxidation state, sixth-ligand identity
(O₂ vs H₂O), and crystal packing differ enough that even after a
high-resolution shift-field fit the diff map is dominated by a ±38σ
peak at the heme Fe — the strongest single feature in any of these
examples.

## Run

```sh
sh runme.sh
```

Single-core wall time: ~5–10 min.

## What `runme.sh` does

1. `wget` PDB and SF cif for 1mbo and 1a6m from RCSB
2. `gemmi cif2mtz` each SF cif
3. `ccp4-python ../../bendfinder.py 1mbo.pdb 1a6m.pdb 1mbo.mtz 1a6m.mtz run_refinement scan_dir=scan_fitreso`

## Expected results

| scan point | RMSD(CA) | Rbent | Rbend | k | B (Å²) | top peak |
|---|---|---|---|---|---|---|
| hkl00 | 0.308 | 55.7% | 0.0% | ~0.99 | ~−1 | −32.8σ A/155HEM/FE(m) |
| fr10  | 0.091 | 52.6% | 23.2% | ~0.96 | ~−2 | −41.0σ A/155HEM/FE(m) |
| **fr5** | **0.063** | **49.7%** | 39.4% | ~0.93 | ~−2 | −38.5σ A/155HEM/FE(m) |

The Fe peak at A/155HEM/FE is genuinely huge — Fe has Z=26 vs the
surrounding carbons/nitrogens, and its position/occupancy differs
slightly between the two crystal forms.  No amount of shift-field
bending will absorb this because the difference is chemical
(oxidation state + ligand identity), not spatial.  The persistent
~50 % Rbent reflects the same: these are two genuinely different
crystal states of the same protein.
