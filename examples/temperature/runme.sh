#!/bin/sh
# Temperature: 4kjj (cryo) -> 4kjk (RT).  P 2_1 2_1 2_1.
set -e

for code in 4kjj 4kjk; do
    [ -f ${code}.pdb    ] || wget -q https://files.rcsb.org/download/${code}.pdb
    [ -f ${code}-sf.cif ] || wget -q https://files.rcsb.org/download/${code}-sf.cif
done

# SF cifs for this pair have intensities only (no F).  Convert via gemmi,
# then ctruncate IMEAN/SIGIMEAN -> F/SIGF.
for code in 4kjj 4kjk; do
    [ -f ${code}_I.mtz ] || ccp4-python -c "
import gemmi
doc = gemmi.cif.read('${code}-sf.cif')
mtz = gemmi.CifToMtz().convert_block_to_mtz(gemmi.as_refln_blocks(doc)[0])
mtz.write_to_file('${code}_I.mtz')
"
    [ -f ${code}.mtz ] || ctruncate -mtzin ${code}_I.mtz -mtzout ${code}.mtz \
                                    -colin '/*/*/[IMEAN,SIGIMEAN]' \
                                    > ${code}_ctruncate.log 2>&1
done

ccp4-python ../../bendfinder.py 4kjj.pdb 4kjk.pdb 4kjj.mtz 4kjk.mtz \
    run_refinement scan_dir=scan_fitreso
