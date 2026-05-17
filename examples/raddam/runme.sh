#!/bin/sh
# Lysozyme radiation damage dose series: 5kxk (undamaged) → 5kxl/m/n.
# Three independent scans with subtract=bent so positive peaks =
# density appearing with dose.
set -e

CODES="5kxk 5kxl 5kxm 5kxn"
DAMAGED="5kxl 5kxm 5kxn"

# 1. Download deposited PDBs and structure-factor cifs
for code in $CODES; do
    [ -f ${code}.pdb    ] || wget -q https://files.rcsb.org/download/${code}.pdb
    [ -f ${code}-sf.cif ] || wget -q https://files.rcsb.org/download/${code}-sf.cif
done

# 2. Convert SF cifs → MTZ (gemmi).  RCSB provides intensities only
#    for this series; we ctruncate below to get F/SIGF for refmac.
for code in $CODES; do
    [ -f ${code}_I.mtz ] || ccp4-python -c "
import gemmi
doc = gemmi.cif.read('${code}-sf.cif')
mtz = gemmi.CifToMtz().convert_block_to_mtz(gemmi.as_refln_blocks(doc)[0])
mtz.write_to_file('${code}_I.mtz')
"
    if [ ! -f ${code}.mtz ]; then
        ctruncate -mtzin ${code}_I.mtz -mtzout ${code}.mtz \
                  -colin '/*/*/[IMEAN,SIGIMEAN]' > ${code}_ctruncate.log 2>&1
    fi
done

# 3. Run bendfinder three times (mov=5kxk, ref=5kxl/m/n)
for ref in $DAMAGED; do
    ccp4-python ../../bendfinder.py 5kxk.pdb ${ref}.pdb 5kxk.mtz ${ref}.mtz \
        run_refinement subtract=bent scan_dir=scan_fitreso_${ref}
done

echo
echo "Done.  Compare scan_fitreso_5kx{l,m,n}/scan_fitreso.log to the"
echo "reference scan_fitreso_5kx{l,m,n}.log files shipped here."
