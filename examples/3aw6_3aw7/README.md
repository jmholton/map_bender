# Example: lysozyme 3aw7 → 3aw6

Hen egg-white lysozyme measured at two relative humidity conditions, both in the same space group but with slightly different unit cells due to crystal dehydration:

| File | PDB ID | Space group | a=b (Å) | c (Å) | RH |
|------|--------|-------------|---------|-------|----|
| `3aw7_refine_001.pdb` | 3aw7 | P4₃2₁2 | 77.0 | 37.3 | 84.2% |
| `3aw6_refine_001.pdb` | 3aw6 | P4₃2₁2 | 79.0 | 38.2 | 71.9% |

The unit cell contracts ~2.5% isotropically as humidity decreases from 84.2% to 71.9%. The shift field captures both this bulk contraction and the local conformational changes that accompany it.

## Files

| File | Description |
|------|-------------|
| `3aw7_refine_001.pdb` | Moving PDB — the structure to be bent |
| `3aw6_refine_001.pdb` | Reference PDB — the target frame |
| `3aw7_2fofc.map` | 2Fo-Fc electron density map for 3aw7 (CCP4 format; any map can be used) |
| `fitparams_30.gnuplot` | Pre-fitted Fourier coefficients at nhkls=30 |
| `bent30_3aw7_refine_001.pdb` | 3aw7 after applying the nhkls=30 shift field |
| `bendfinder_nhkls30.log` | Full output from the nhkls=30 run |
| `runme.csh` | Script to reproduce the fit from scratch |

## Quick start

**Apply the pre-fitted solution** (seconds, no gnuplot needed):

```tcsh
bendfinder.com 3aw7_refine_001.pdb 3aw6_refine_001.pdb fitparams_30.gnuplot
```

This applies the pre-computed shift field and reports RMSD without re-fitting.

**Apply to a map as well** (any CCP4 map in the frame of the moving PDB):

```tcsh
bendfinder.com 3aw7_refine_001.pdb 3aw6_refine_001.pdb 3aw7_2fofc.map fitparams_30.gnuplot
```

**Reproduce the full fit from scratch** (~25 minutes):

```tcsh
cd /tmp/my_test
cp /path/to/examples/3aw6_3aw7/*.pdb .
/path/to/bendfinder.com 3aw7_refine_001.pdb 3aw6_refine_001.pdb nhkls=30 geotest=false
```

Or simply run the provided script from a scratch directory:

```tcsh
mkdir /tmp/my_test && cd /tmp/my_test
/path/to/examples/3aw6_3aw7/runme.csh
```

## Expected results

RMSD(CA) against `3aw6_refine_001.pdb` at each iteration (129 Cα pairs):

| nhkls | RMSD(CA) | Elapsed |
|-------|----------|---------|
| 0     | 0.771 Å  | —       |
| 5     | 0.327 Å  | 9 s     |
| 9     | 0.277 Å  | 34 s    |
| 15    | 0.268 Å  | 144 s   |
| 20    | 0.245 Å  | 580 s   |
| 26    | 0.215 Å  | 1070 s  |
| 30    | 0.211 Å  | 1503 s  |

The full log is in `bendfinder_nhkls30.log`.

## Notes on this example

- The initial RMSD of 0.771 Å (all-atom) drops to 0.327 Å with just 5 Fourier terms because the dominant shift is a rigid-body cell expansion, captured by the (0,1,0) and equivalent low-resolution terms.
- By nhkls=30 the RMSD(CA) reaches 0.211 Å, comparable to the original result (0.209 Å) obtained with 91 positive-only HKLs in ~49 minutes.
- The `drop_snr=1` default means statistically insignificant Fourier terms (|a|/σ < 1) are pruned after each fit, preventing dirty-beam ghost terms from accumulating.
