#!/usr/bin/env ccp4-python
"""make_test_pdbs.py — produce the two synthetic deformations of 7rsa
that the Magdoff-Crick approximation predicts the bendfinder shift
field should recover.

  Test 1: 7rsa_scaled.pdb   — CRYST1 a,b,c scaled by 1.005 (atoms unchanged).
                              Pure isomorphous cell expansion.
  Test 2: 7rsa_rotated.pdb  — atoms rotated 0.5° about a random axis
                              (RNG seed 42).  Not SG-consistent.

Both are written in the script's own directory (alongside 7rsa.pdb,
which runme.sh wgets from RCSB).  No external deps beyond numpy.
"""
import os, math
import numpy as np

HERE    = os.path.dirname(os.path.abspath(__file__))
PDB_REF = os.path.join(HERE, '7rsa.pdb')
SCALE   = 1.005
ROT_DEG = 0.5
SEED    = 42


def modify_pdb(src, dst, fn):
    """Copy a PDB, applying fn(xyz) -> xyz to each ATOM/HETATM record."""
    with open(src) as f, open(dst, 'w') as g:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                xyz = np.array([float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])])
                xp, yp, zp = fn(xyz)
                line = (line[:30] + f"{xp:8.3f}{yp:8.3f}{zp:8.3f}" + line[54:])
            g.write(line)


def scale_cell_pdb(src, dst, scale):
    """Scale CRYST1 a,b,c by `scale`; strip SCALE/ORIGX so gemmi re-derives
    the fractional matrix from the (modified) cell instead of using the
    stale matrix recorded in the SCALE cards."""
    with open(src) as f, open(dst, 'w') as g:
        for line in f:
            if line.startswith('CRYST1'):
                a = float(line[6:15]) * scale
                b = float(line[15:24]) * scale
                c = float(line[24:33]) * scale
                line = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{line[33:]}"
            elif line.startswith(('SCALE', 'ORIGX')):
                continue
            g.write(line)


def rotation_matrix(axis, angle_deg):
    """Rodrigues rotation matrix (axis is a 3-vector, angle in degrees)."""
    n = np.asarray(axis, float); n /= np.linalg.norm(n)
    th = math.radians(angle_deg)
    c, s = math.cos(th), math.sin(th)
    K = np.array([[0, -n[2], n[1]],
                  [n[2], 0, -n[0]],
                  [-n[1], n[0], 0]])
    return np.eye(3)*c + s*K + (1-c)*np.outer(n, n)


# Test 1 — isomorphous cell scaling
scale_cell_pdb(PDB_REF, os.path.join(HERE, '7rsa_scaled.pdb'), SCALE)
print(f'wrote 7rsa_scaled.pdb   (CRYST1 a,b,c x {SCALE}; atoms unchanged)')

# Test 2 — rigid-body rotation
rng     = np.random.default_rng(SEED)
axis    = rng.standard_normal(3); axis /= np.linalg.norm(axis)
R       = rotation_matrix(axis, ROT_DEG)
modify_pdb(PDB_REF, os.path.join(HERE, '7rsa_rotated.pdb'), lambda r: R @ r)
print(f'wrote 7rsa_rotated.pdb  ({ROT_DEG} deg about axis '
      f'[{axis[0]:+.3f},{axis[1]:+.3f},{axis[2]:+.3f}], seed={SEED})')
