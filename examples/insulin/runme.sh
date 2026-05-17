#!/bin/sh
# Insulin T->R transition: 4fg3 (T) -> 4e7u (R), H 3 hexamer.
# Triggers altindex resolution + re-refinement inside bendfinder.
set -e

for code in 4fg3 4e7u; do
    [ -f ${code}.pdb    ] || wget -q https://files.rcsb.org/download/${code}.pdb
    [ -f ${code}-sf.cif ] || wget -q https://files.rcsb.org/download/${code}-sf.cif
done

for code in 4fg3 4e7u; do
    [ -f ${code}.mtz ] || ccp4-python -c "
import gemmi
doc = gemmi.cif.read('${code}-sf.cif')
mtz = gemmi.CifToMtz().convert_block_to_mtz(gemmi.as_refln_blocks(doc)[0])
mtz.write_to_file('${code}.mtz')
"
done

ccp4-python ../../bendfinder.py 4fg3.pdb 4e7u.pdb 4fg3.mtz 4e7u.mtz \
    run_refinement scan_dir=scan_fitreso
