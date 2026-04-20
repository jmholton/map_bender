#! /bin/tcsh -f
#
# Example: map the lysozyme crystal form 3aw7 onto 3aw6
#
# 3aw7_refine_001.pdb  -- moving PDB  (P43212, a=b=77.0, c=37.3 A)
# 3aw6_refine_001.pdb  -- reference PDB (P43212, a=b=79.0, c=38.2 A)
#
# Both structures are refined models of hen egg-white lysozyme in different
# crystal forms with slightly different unit cells.  The coordinate shift
# field is dominated by a rigid-body expansion of the cell.
#
# Expected result (nhkls=30, ~25 minutes, RMSD converges to ~0.21 A CA):
#   nhkls=5  RMSD(CA)=0.327 A   (9 s)
#   nhkls=20 RMSD(CA)=0.245 A  (580 s)
#   nhkls=30 RMSD(CA)=0.211 A (1503 s)
#
# To just apply the pre-fitted solution (< 1 minute):
#   bendfinder.com 3aw7_refine_001.pdb 3aw6_refine_001.pdb fitparams_30.gnuplot
#
# To also re-sample the 2Fo-Fc map into the reference frame:
#   bendfinder.com 3aw7_refine_001.pdb 3aw6_refine_001.pdb 3aw7_2fofc.map fitparams_30.gnuplot

set script = `dirname $0`/../../bendfinder.com

$script \
    3aw7_refine_001.pdb \
    3aw6_refine_001.pdb \
    nhkls=30 \
    geotest=false
