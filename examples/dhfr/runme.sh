#!/bin/sh
# DHFR ligand binding: 1rx2 (with FOL) -> 1rx1 (without FOL).
# Downloads inputs from RCSB, converts SF cifs to MTZ, runs bendfinder.
set -e

# 1. Download deposited PDBs and structure-factor cifs
for code in 1rx1 1rx2; do
    [ -f ${code}.pdb    ] || wget -q https://files.rcsb.org/download/${code}.pdb
    [ -f ${code}-sf.cif ] || wget -q https://files.rcsb.org/download/${code}-sf.cif
done

# 2. Convert SF cifs → MTZ (gemmi)
for code in 1rx1 1rx2; do
    [ -f ${code}.mtz ] || ccp4-python -c "
import gemmi
doc = gemmi.cif.read('${code}-sf.cif')
mtz = gemmi.CifToMtz().convert_block_to_mtz(gemmi.as_refln_blocks(doc)[0])
mtz.write_to_file('${code}.mtz')
"
done

# 3. Run bendfinder
#    Default subtract=ref so FOL (in 1rx2 mov but not 1rx1 ref) appears
#    as positive (green) difference density.
ccp4-python ../../bendfinder.py 1rx2.pdb 1rx1.pdb 1rx2.mtz 1rx1.mtz \
    run_refinement scan_dir=scan_fitreso

echo
echo "Done.  Compare scan_fitreso/scan_fitreso.log to the reference"
echo "scan_fitreso.log shipped in this directory."
