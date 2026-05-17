# Lysozyme — radiation damage dose series

| field | value |
|---|---|
| moving (`mov`)  | **5kxk** — undamaged, lowest dose |
| reference (`ref`) | **5kxl** / **5kxm** / **5kxn** — increasingly damaged |
| space group | P 4₃ 2₁ 2 |
| CA pairs | ~976–992 (per scan) |
| sign convention | `subtract=bent` → diff = ref − bent |

Lysozyme dose series from a single crystal: 5kxk is the lowest-dose
(effectively undamaged) reference; 5kxl/5kxm/5kxn are the same crystal
re-collected at successively higher absorbed dose.  Sign convention is
inverted (`subtract=bent`) so **positive peaks = density appearing with
dose** (radiolysis products, water reorganization), **negative peaks =
density disappearing with dose** (broken disulfides, leaving waters).

Dose ordering (alphabetical = chronological): **5kxk < 5kxl < 5kxm < 5kxn**.

## Run

```sh
sh runme.sh
```

This downloads all four PDBs + SF cifs, ctruncates intensities → F
(RCSB deposits IMEAN/SIGIMEAN for this series, not FP/SIGFP), then
launches three independent scans (one per damaged dataset).

Single-core wall time: ~10 min × 3 = ~30 min serial.  The three scans
are independent and parallelizable.

## What `runme.sh` does

1. `wget` PDB + SF cif for 5kxk, 5kxl, 5kxm, 5kxn
2. `gemmi cif2mtz` each SF cif (intensities only)
3. `ctruncate` IMEAN/SIGIMEAN → F/SIGF
4. Run `bendfinder.py` three times:
   - `5kxk → 5kxl`  →  `scan_fitreso_5kxl/`
   - `5kxk → 5kxm`  →  `scan_fitreso_5kxm/`
   - `5kxk → 5kxn`  →  `scan_fitreso_5kxn/`
   each with `subtract=bent` so positive peaks = density appearing with dose.

## Expected results

| dose | fr5 RMSD | fr5 Rbent | fr5 top peak | observation |
|---|---|---|---|---|
| **5kxl** | 0.086 Å | 22.1% | −15.9σ A/394HOH/O(r) | first dose: persistent CYS/SG damage (−10σ at A/30CYS, A/115CYS) |
| **5kxm** | 0.047 Å | 19.2% | **+22.4σ A/205CL/CL(m)** | mid-dose: chloride accumulates as radiolysis product |
| **5kxn** | 0.051 Å | 24.6% | +13.2σ A/398HOH/O(r) | highest dose: progressive CYS/SG damage (−15σ at A/30CYS, A/94CYS), water reorganization |

Cross-cutting findings:

- **Disulfide damage signature**: A/30CYS/SG and A/94CYS/SG (and others)
  consistently appear as strong **negative** peaks across all three doses
  — sulfur density present in undamaged 5kxk is lost in the damaged
  references.  Magnitude grows with dose:
    * 5kxl: −10σ at A/30CYS
    * 5kxn: −13 to −15σ at A/30CYS, A/94CYS
- **Chloride ion** at A/205CL/CL(m): in 5kxm this becomes a striking
  +22σ positive feature by fr5 — a chloride is recruited or strongly
  ordered at this site as dose accumulates.
- **Water reorganization**: water-O atoms dominate the strongest
  surface peaks at the highest resolutions (fr7–fr5).  Both new positive
  peaks (waters arriving) and strong negative peaks (waters leaving)
  appear.
- 5kxm gives the cleanest fr5 RMSD (0.047 Å); 5kxn gives the largest
  Rbent (24.6%) consistent with more dose damage.

Reference logs for all three doses are checked in as
`scan_fitreso_5kx{l,m,n}.log` for comparison.

## Looking at the maps in Coot

```sh
coot scan_fitreso_5kxm/fr5/ref.pdb scan_fitreso_5kxm/fr5/bent.pdb \
     --auto scan_fitreso_5kxm/fr5/ref.mtz \
     --auto scan_fitreso_5kxm/fr5/bent.mtz
```

Navigate to A/205CL/CL to see the +22σ chloride accumulation, or to any
CYS/SG to see the disulfide damage signature.
