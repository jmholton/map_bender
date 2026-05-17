#!/bin/sh
# Magdoff-Crick synthetic deformation test on ribonuclease A (7rsa).
# 1. wget 7rsa.pdb from RCSB
# 2. python script creates two synthetic deformations
# 3. bendfinder.py recovers each (constrained P 2_1 fit, fitreso 7 A)
set -e

[ -f 7rsa.pdb ] || wget -q https://files.rcsb.org/download/7rsa.pdb

ccp4-python ./make_test_pdbs.py

# Test 1: cell scaling 0.5%
echo
echo "================  Test 1: cell scaling x 1.005 ================"
ccp4-python ../../bendfinder.py 7rsa_scaled.pdb 7rsa.pdb fitreso=7

# Test 2: rigid rotation 0.5 deg
echo
echo "================  Test 2: rigid rotation 0.5 deg ================"
ccp4-python ../../bendfinder.py 7rsa_rotated.pdb 7rsa.pdb fitreso=7
