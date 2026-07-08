#!/bin/tcsh -f
#
# run_all_tests.com — run the full gamut of bendfinder tests and examples.
#
# Step 0 (always runs first, idempotent):
#   Download deposited PDB + sf-cif files from RCSB and cif2mtz to MTZ
#   into per-system dirs.  Skips any file that already exists, so re-runs
#   are fast.  This is the "fresh clone" entry point — running this script
#   in a directory with just bendfinder.py / altalign.py / test_symm_all_sgs.py
#   is sufficient.
#
# Tests:
#   1. test_symm_all_sgs.py     — 65 Sohncke SGs symmetry constraint check
#   2. magdoff/test_magdoff.py  — synthetic 0.5deg deformation recovery
#
# Examples (fitreso_scan from raw, scan_fitreso/ subdir per system):
#   3. lyso     3aw6 -> 3aw7   P4(3)2(1)2
#   4. dhfr     1rx2 -> 1rx1   P2(1)2(1)2(1)
#   5. raddam   5kxk -> 5kxl   P4(3)2(1)2   subtract=bent
#   6. raddam   5kxk -> 5kxm   P4(3)2(1)2   subtract=bent
#   7. raddam   5kxk -> 5kxn   P4(3)2(1)2   subtract=bent
#   8. myoglobin 1mbo -> 1a6m  P2(1)
#   9. insulin  4fg3 -> 4e7u   H3           fill_asu=True
#  10. lipox    9o4s -> 9o4t   P2(1)        same crystal, just distorted
#                                            (~4% metric drift in a/b/c);
#                                            exercises stretch + loose-tol
#                                            altindex reindex of mov Fobs
#  11. porin    altalign 3poq->3pou + refmac on the R32:R output
#
# Each test writes to test_results_<timestamp>/<name>.log and the script
# prints a PASS/FAIL summary at the end.  Run from the working area
# (the directory holding bendfinder.py, altalign.py, lyso/, dhfr/, ...).
#
# Usage:  ./run_all_tests.com

if (! $?CCP4) source /programs/ccp4-9/bin/ccp4.setup-csh

# CCP4-9 ships single-threaded blas/cblas/lapack — SVDs in bend_fit are 5-10x
# slower than they need to be.  LD_PRELOAD the pthreaded openblas bundled
# with phenix (libopenblasp 0.3.25).  Pair with `srun -c $SLURM_CPUS` below
# so OPENBLAS_NUM_THREADS actually sees the cores (default srun gives 1).
set NCPUS = 64
if (-e /programs/phenix-2.1rc2-6037/lib/libopenblas.so.0) then
    setenv LD_PRELOAD /programs/phenix-2.1rc2-6037/lib/libopenblas.so.0
    setenv OPENBLAS_NUM_THREADS $NCPUS
endif

# CWD sanity check — this script just needs bendfinder.py + altalign.py +
# test_symm_all_sgs.py in the CWD.  The per-system input dirs (lyso/, dhfr/,
# ...) are created and populated by the step-0 setup below.  If invoked from
# inside map_bender/ (the git subdir), cd up; otherwise fail loudly so we
# don't fake PASS off "no such file" grep output.
if (! -e bendfinder.py || ! -e altalign.py || ! -e test_symm_all_sgs.py) then
    if (-e ../bendfinder.py && -e ../altalign.py && -e ../test_symm_all_sgs.py) then
        echo "(cd .. — running from working area, not map_bender/)"
        cd ..
    else
        echo "ERROR: run from a directory containing bendfinder.py, altalign.py, test_symm_all_sgs.py."
        echo "       CWD=`pwd`"
        exit 2
    endif
endif

set ts     = `date +%Y%m%d_%H%M%S`
set logdir = test_results_$ts
mkdir -p $logdir
set summary = $logdir/SUMMARY.txt
echo "logs and per-test details: $logdir/"

@ npass = 0
@ nfail = 0
printf "%-30s %-7s %s\n" "test" "result" "metric" > $summary
printf "%-30s %-7s %s\n" "------------------------------" "------" "------" >> $summary

# ── 0. ensure inputs (download PDBs + sf-cifs from RCSB, cif2mtz) ────────
# Idempotent: skips per-file work whenever the target already exists, so a
# re-run is fast.  Removes a legacy lipox symlink so we land real downloaded
# inputs in lipox/.  magdoff has no deposited SF — only 7rsa.pdb is fetched.
echo ""
echo "==> [ 0/11] ensure inputs (RCSB download + cif2mtz)"
set setup_log = $logdir/00_setup.log
echo "setup log: $setup_log" >& $setup_log
foreach line ( \
    "lyso      3aw6 3aw7" \
    "dhfr      1rx2 1rx1" \
    "raddam    5kxk 5kxl 5kxm 5kxn" \
    "myoglobin 1mbo 1a6m" \
    "insulin   4fg3 4e7u" \
    "porin     3poq 3pou" \
    "lipox     9o4s 9o4t" \
    "magdoff   7rsa" \
)
    set parts = ( $line )
    set sysd  = $parts[1]
    if (-l $sysd) rm $sysd                      # drop legacy symlink (lipox)
    mkdir -p $sysd
    foreach pdbid ( $parts[2-] )
        set pdb = $sysd/$pdbid.pdb
        if (! -e $pdb) then
            echo "    [$sysd] fetching $pdbid.pdb"
            curl -s -f -o $pdb https://files.rcsb.org/download/$pdbid.pdb
            if ($status != 0) then
                echo "    ERROR: curl failed for $pdbid.pdb" | tee -a $setup_log
                rm -f $pdb
            endif
        endif
        if ($sysd == magdoff) continue          # magdoff: no SF needed
        set cif = $sysd/$pdbid-sf.cif
        set mtz = $sysd/$pdbid.mtz
        if (! -e $cif) then
            echo "    [$sysd] fetching $pdbid-sf.cif"
            curl -s -f -o $cif https://files.rcsb.org/download/$pdbid-sf.cif
            if ($status != 0) then
                echo "    ERROR: curl failed for $pdbid-sf.cif" | tee -a $setup_log
                rm -f $cif
            endif
        endif
        if (! -e $mtz && -e $cif) then
            echo "    [$sysd] cif2mtz $pdbid-sf.cif → $pdbid.mtz"
            echo END | cif2mtz hklin $cif hklout $mtz >>& $setup_log
        endif
    end
end
# magdoff/test_magdoff.py is not currently shipped in the git repo; until
# it is, copy from old/magdoff (where the previous gamut lived) if present.
if (! -e magdoff/test_magdoff.py && -e old/magdoff/test_magdoff.py) then
    cp old/magdoff/test_magdoff.py magdoff/test_magdoff.py
    echo "    [magdoff] copied test_magdoff.py from old/magdoff/" | tee -a $setup_log
endif
echo "    setup done"

# ── 1. test_symm_all_sgs ─────────────────────────────────────────────────
echo ""
echo "==> [ 1/11] test_symm_all_sgs.py"
srun -c $NCPUS --job-name=test_symm ccp4-python test_symm_all_sgs.py >& $logdir/01_test_symm.log
set r = `grep "Constrained fit:" $logdir/01_test_symm.log | head -1`
if ("$r" =~ *"65/65"*) then
    echo "    PASS  $r"
    @ npass++
    printf "%-30s %-7s %s\n" "test_symm_all_sgs" PASS "$r" >> $summary
else
    echo "    FAIL  expected '65/65' line"
    @ nfail++
    printf "%-30s %-7s %s\n" "test_symm_all_sgs" FAIL "no_65/65_line" >> $summary
endif

# ── 2. magdoff Tests 1 & 2 ──────────────────────────────────────────────
# Two tests run; each prints "Constrained RMSD(CA) after fit: 0.0XXX Å ...".
# Report Test 2 (rigid-body rotation) only — it's the canonical entry in
# CLAUDE.md and the harder of the two.  tail -1 selects Test 2's lines.
echo ""
echo "==> [ 2/11] magdoff/test_magdoff.py"
if (! -e magdoff/test_magdoff.py) then
    echo "    FAIL  magdoff/test_magdoff.py missing"
    @ nfail++
    printf "%-30s %-7s no_script\n" "magdoff" FAIL >> $summary
else
    ( cd magdoff && srun -c $NCPUS --job-name=magdoff ccp4-python test_magdoff.py ) >& $logdir/02_magdoff.log
    set rcon = `grep "RMSD(CA) after fit:" $logdir/02_magdoff.log | grep "Constrained" | tail -1 | awk '{print $5}'`
    set runc = `grep "RMSD(CA) after fit:" $logdir/02_magdoff.log | grep "Unconstrained" | tail -1 | awk '{print $5}'`
    if ("$rcon" != "" && "$runc" != "") then
        echo "    PASS  Test 2  constrained=$rcon  unconstrained=$runc"
        @ npass++
        printf "%-30s %-7s Test2  constrained=%s  unconstrained=%s\n" "magdoff" PASS "$rcon" "$runc" >> $summary
    else
        echo "    FAIL  could not parse RMSDs from log"
        @ nfail++
        printf "%-30s %-7s parse_error\n" "magdoff" FAIL >> $summary
    endif
endif

# Reusable: run one fitreso_scan example.  Args via env vars (tcsh has no
# functions; using a foreach loop with positional vars per iteration).
#
# For each example: spec = "label sys mov_pdb ref_pdb mov_mtz ref_mtz subtract fill_asu tag"
# tag goes in the scan_dir name (so raddam variants don't collide).

set i = 2
# fill_asu=True for every example because the deposited reference MTZs
# sit below the 99% SG-ASU completeness gate (lyso 3aw7 is 90.2%, etc.).
foreach spec ( \
    "lyso        lyso        3aw6.pdb 3aw7.pdb 3aw6.mtz 3aw7.mtz ref  True ''" \
    "dhfr        dhfr        1rx2.pdb 1rx1.pdb 1rx2.mtz 1rx1.mtz ref  True ''" \
    "raddam_5kxl raddam      5kxk.pdb 5kxl.pdb 5kxk.mtz 5kxl.mtz bent True _5kxl" \
    "raddam_5kxm raddam      5kxk.pdb 5kxm.pdb 5kxk.mtz 5kxm.mtz bent True _5kxm" \
    "raddam_5kxn raddam      5kxk.pdb 5kxn.pdb 5kxk.mtz 5kxn.mtz bent True _5kxn" \
    "myoglobin   myoglobin   1mbo.pdb 1a6m.pdb 1mbo.mtz 1a6m.mtz ref  True ''" \
    "insulin     insulin     4fg3.pdb 4e7u.pdb 4fg3.mtz 4e7u.mtz ref  True ''" \
    "lipox       lipox       9o4s.pdb 9o4t.pdb 9o4s.mtz 9o4t.mtz ref  True ''" \
)
    set parts = ( $spec )
    set label = $parts[1]
    set sys   = $parts[2]
    set movp  = $parts[3]
    set refp  = $parts[4]
    set movm  = $parts[5]
    set refm  = $parts[6]
    set sub   = $parts[7]
    set fill  = $parts[8]
    set tag   = $parts[9]
    if ("$tag" == "''") set tag = ""
    set scan = scan_fitreso$tag
    @ i++
    set log = $logdir/`printf "%02d" $i`_$label.log

    echo ""
    echo "==> [`printf '%2d' $i`/11] $label  ($sys $movp -> $refp, fill_asu=$fill, subtract=$sub)"
    ( cd $sys && rm -rf $scan && srun -c $NCPUS --job-name=$label ccp4-python -c "import sys; sys.path.insert(0,'..'); from bendfinder import fitreso_scan; fitreso_scan(mov_pdb='$movp', ref_pdb='$refp', mov_mtz='$movm', ref_mtz='$refm', scan_dir='$scan', run_refinement_flag=True, refine_cycles=5, fill_asu=$fill, subtract='$sub', scan_all_fr=True)" ) >& $log

    if (! -e $sys/$scan/scan_fitreso.log) then
        echo "    FAIL  $sys/$scan/scan_fitreso.log not written (see $log)"
        @ nfail++
        printf "%-30s %-7s no_scan_log\n" "$label" FAIL >> $summary
        continue
    endif
    set fr5 = `grep '^    fr5' $sys/$scan/scan_fitreso.log | head -1`
    set bst = `grep '^   best' $sys/$scan/scan_fitreso.log | head -1`
    set dopt = `grep '# best row:' $sys/$scan/scan_fitreso.log | head -1 | sed 's/.*d_opt = //' | awk '{print $1}'`
    if ("$fr5" != "") then
        set rmsd  = `echo "$fr5" | awk '{print $2}'`
        set rbent = `echo "$fr5" | awk '{print $4}'`
        set best_str = ""
        if ("$bst" != "") then
            set b_rmsd  = `echo "$bst" | awk '{print $2}'`
            set b_rbent = `echo "$bst" | awk '{print $4}'`
            set best_str = "  best  RMSD=$b_rmsd  Rbent=$b_rbent  d_opt=${dopt}A"
        endif
        echo "    PASS  fr5  RMSD=$rmsd  Rbent=$rbent$best_str"
        @ npass++
        printf "%-30s %-7s fr5  RMSD=%s  Rbent=%s%s\n" "$label" PASS "$rmsd" "$rbent" "$best_str" >> $summary
    else
        echo "    FAIL  no fr5 row in scan_fitreso.log (see $sys/$scan/scan_fitreso.log)"
        @ nfail++
        printf "%-30s %-7s no_fr5_row\n" "$label" FAIL >> $summary
    endif
end

# ── 11. porin altalign + refmac on R32:R ─────────────────────────────────
@ i++
echo ""
echo "==> [`printf '%2d' $i`/11] porin  altalign 3poq->3pou + refmac on R32:R"
set plog = $logdir/`printf "%02d" $i`_porin_altalign.log
( cd porin && rm -rf altalign && mkdir altalign && \
    srun -c $NCPUS --job-name=porin_aa ccp4-python ../altalign.py 3poq.pdb 3pou.pdb \
        altalign/aa.pdb 3poq.mtz altalign/aa.mtz ) >& $plog

set residual = ""
if (-e $plog) set residual = `grep 'honest discrete CA RMSD' $plog | head -1 | awk '{print $7}'`
if (-e porin/altalign/aa_R32R.pdb && -e porin/altalign/aa_R32R.mtz) then
    set rlog = $logdir/`printf "%02d" $i`_porin_refmac.log
    srun -c $NCPUS --job-name=porin_refmac ccp4-python -c "import sys; sys.path.insert(0,'.'); from bendfinder import run_refinement; run_refinement('porin/altalign/aa_R32R.pdb', 'porin/altalign/aa_R32R.mtz', outdir='porin/altalign/refine', n_cycles=5, fill_asu=True)" >& $rlog

    set rf = ""
    if (-d porin/altalign/refine) then
        # refmac prints BOTH "Overall weighted R factor" and "Overall weighted R2
        # factor"; specifying "R factor" excludes the R2 line.
        set rf = `grep 'Overall weighted R factor' porin/altalign/refine/*refmac.log | tail -1 | awk '{print $6}'`
    endif
    if ("$rf" != "") then
        echo "    PASS  altalign residual=${residual} A   refmac R=$rf"
        @ npass++
        printf "%-30s %-7s altalign=%sA  refmac=R%s\n" "porin (altalign+R32:R)" PASS "$residual" "$rf" >> $summary
    else
        echo "    FAIL  refmac did not produce R-factor (see $rlog)"
        @ nfail++
        printf "%-30s %-7s no_refmac_R\n" "porin (altalign+R32:R)" FAIL >> $summary
    endif
else
    echo "    FAIL  altalign did not write R32R outputs (see $plog)"
    @ nfail++
    printf "%-30s %-7s no_R32R_outputs\n" "porin (altalign+R32:R)" FAIL >> $summary
endif

# ── summary ──────────────────────────────────────────────────────────────
echo ""
echo "===================================================================="
cat $summary
echo "===================================================================="
echo "Total: $npass passed, $nfail failed."
echo "Per-test logs: $logdir/"
if ($nfail > 0) exit 1
exit 0
