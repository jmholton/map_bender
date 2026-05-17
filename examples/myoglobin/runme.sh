#!/bin/sh
# Myoglobin: 1mbo (oxy) -> 1a6m (met).  Both sperm-whale Mb in P 2_1.
set -e

for code in 1mbo 1a6m; do
    [ -f ${code}.pdb    ] || wget -q https://files.rcsb.org/download/${code}.pdb
    [ -f ${code}-sf.cif ] || wget -q https://files.rcsb.org/download/${code}-sf.cif
done

for code in 1mbo 1a6m; do
    [ -f ${code}.mtz ] || ccp4-python -c "
import gemmi
doc = gemmi.cif.read('${code}-sf.cif')
mtz = gemmi.CifToMtz().convert_block_to_mtz(gemmi.as_refln_blocks(doc)[0])
mtz.write_to_file('${code}.mtz')
"
done

ccp4-python ../../bendfinder.py 1mbo.pdb 1a6m.pdb 1mbo.mtz 1a6m.mtz \
    run_refinement scan_dir=scan_fitreso
