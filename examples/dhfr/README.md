# DHFR — ligand binding (FOL)

| field | value |
|---|---|
| moving (`mov`)  | **1rx2** — apo + folate (FOL) + Mn²⁺ + BME |
| reference (`ref`) | **1rx1** — apo + Ca²⁺ (no FOL) |
| space group | P 2₁ 2₁ 2₁ |
| CA pairs | 636 |
| sign convention | `subtract=ref` (default) → diff = bent − ref |

Two crystal forms of *E. coli* DHFR.  1rx2 has the folate (FOL) cofactor
bound; 1rx1 does not.  Both have NAP (NADP analog).  This direction was
chosen so the FOL ligand appears as **positive (green)** difference
density, the standard "what's been added" view.

HETATM residue inventory:

| structure | residues |
|---|---|
| 1rx1 | CA, HOH, NAP |
| 1rx2 | BME, FOL, HOH, MN, NAP |

## Run

```sh
sh runme.sh
```

Single-core wall time: ~10 min.

## What `runme.sh` does

1. `wget` PDB and SF cif for 1rx1 and 1rx2 from RCSB
2. `gemmi cif2mtz` each SF cif
3. `ccp4-python ../../bendfinder.py 1rx2.pdb 1rx1.pdb 1rx2.mtz 1rx1.mtz run_refinement scan_dir=scan_fitreso`

## Expected results

| scan point | RMSD(CA) | Rbent | Rbend | k | B (Å²) | top peak |
|---|---|---|---|---|---|---|
| hkl00 | 0.448 | 43.0% | 0.0% | — | — | −8.4σ A/73THR/O(r) |
| hkl01–fr20 | 0.41–0.21 | 47→38% | 13–28% | ~0.95 | ~−1 | **+7 to +9σ A/161FOL/C14(m)** |
| fr10 | 0.121 | 37.9% | 31.4% | 0.94 | −1.9 | −8.1σ A/300CA/CA(r) |
| **fr5** | **0.070** | **41.5%** | 38.3% | 0.92 | −2.1 | −10.5σ A/300CA/CA(r) |

Key features (with the default `subtract=ref`, positive = density in
bent absent from ref):

- **FOL ligand**: every FOL atom shows positive (green) difference
  density in `scan_fitreso/fr5/bent.mtz` DELFWT.  Strongest atoms at
  fr5: O4 +6.78σ, N3 +6.47σ, C4 +5.94σ, C14 +9σ at lower resolutions.
- **Ca²⁺ at A/300CA**: −10.46σ negative (red) peak — Ca²⁺ is only in
  1rx1 (ref), so the bent map (from 1rx2) has no density here.
- B-factor stays near zero (~−2 Å²) throughout the scan: rigid-body
  refmac keeps both FWT outputs on the same scale.

The bending only marginally reduces Rbent (43% → 42%) because these
crystal forms are genuinely different — the diff map is dominated by
ligand/ion peaks, not by spatial mismatch.

The reference output is checked in as `scan_fitreso.log` for comparison.

## Looking at the maps in Coot

```sh
coot scan_fitreso/fr5/ref.pdb scan_fitreso/fr5/bent.pdb \
     --auto scan_fitreso/fr5/ref.mtz \
     --auto scan_fitreso/fr5/bent.mtz
```

Navigate to residue **FOL/161** to see the +6 to +9σ positive density
that the ligand binding site picks up.
