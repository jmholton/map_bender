# Lysozyme — humidity-driven cell change

| field | value |
|---|---|
| moving (`mov`)  | **3aw6** — lower humidity |
| reference (`ref`) | **3aw7** — higher humidity |
| space group | P 4₃ 2₁ 2 |
| CA pairs | 1008 |
| sign convention | `subtract=ref` (default) → diff = bent − ref |

Classic Magdoff scenario: same protein, slightly different unit cell
(~2.5% cell change between humidity points).  Almost all the difference
density is absorbed by the bending; what remains is real waters near the
protein surface.

## Run

Single command — `runme.sh` downloads the deposited PDBs and SF cifs from
RCSB, converts the SF cifs to MTZ, then invokes `bendfinder.py` with
`run_refinement` (refmac rigid body generates FWT/PHWT).

```sh
sh runme.sh
```

Single-core wall time: ~10 min.

## What `runme.sh` does

1. `wget` PDB and SF cif for 3aw6 and 3aw7 from RCSB
2. `gemmi cif2mtz` each SF cif (cifs have FP/SIGFP already, no ctruncate needed)
3. `ccp4-python ../../bendfinder.py 3aw6.pdb 3aw7.pdb 3aw6.mtz 3aw7.mtz run_refinement scan_dir=scan_fitreso`

## Outputs (in `scan_fitreso/`)

- `scan_fitreso.log` — table of RMSD / Rbent / Rbend / k / B / top peak per scan point
- `ref.mtz`, `unbent.mtz` — refmac-refined FWT/PHWT for 3aw7 and 3aw6 (symlinks)
- `unbent.pdb` — moving model in ref frame (before bending)
- `hkl00`, `hkl01..hkl10`, `fr20..fr5` — per-resolution scan subdirs, each containing:
    - `bent.map`, `bent.mtz` (FFT of bent map: FDM/PHIDM + DELFWT/PHDELWT)
    - `diff_norm.map` (z-scored diff, bent − ref in F-space after k+B scaling)
    - `bent.pdb` (moving model after bending at this resolution)
    - `PSDVF.mtz` (shift-field Fourier coefficients)
    - `ref.pdb/mtz`, `unbent.pdb/mtz` (symlinks)
- `refine_mov/`, `refine_ref/` — refmac outputs for 3aw6, 3aw7
- `altindex_resolve/` — empty (no altindex needed for lyso)

A **reference log** is checked in as `scan_fitreso.log` (this directory)
so you can compare your output to a known-good run without actually
running the scan.

## Expected results

| scan point | RMSD(CA) | Rbent | Rbend | k | B (Å²) | top peak |
|---|---|---|---|---|---|---|
| hkl00 | 0.759 | 53.1% | 0.0% | 0.992 | +16.34 | +6.9σ A/80CYS/SG(r) |
| fr10  | 0.069 | 27.0% | 33.6% | 0.969 | +8.31 | −4.8σ A/143HOH/O(r) |
| **fr5** | **0.034** | **29.6%** | 47.3% | 0.954 | +8.10 | +6.9σ A/183HOH/O(r) |

- **B drops from +16 → +8 Å²** as the PSDVF absorbs spatial mismatch; the
  residual +8 is the true thermal B difference between the two crystal
  forms.
- All persistent peaks (±5–7σ at fr5) are on waters/surface atoms —
  there's no protein structural change beyond what bending captures.

## Looking at the maps in Coot

```sh
coot scan_fitreso/fr5/ref.pdb scan_fitreso/fr5/bent.pdb \
     --auto scan_fitreso/fr5/ref.mtz \
     --auto scan_fitreso/fr5/bent.mtz
```

`bent.mtz` carries DELFWT/PHDELWT — Coot loads it as a difference map
automatically.  Contour at ±3σ.  With the default `subtract=ref`,
green = density present in bent (post-PSDVF moving structure) but absent
in ref; red = density present in ref but absent in bent.
