#!/bin/tcsh -f
#
# run_all_tests.com — run the full gamut of bendfinder tests and examples.
#
# Tests:
#   1. test_symm_all_sgs.py     — 65 Sohncke SGs symmetry constraint check
#   2. magdoff/test_magdoff.py  — synthetic 0.5deg deformation recovery
#
# Examples (fitreso_scan from raw, scan_test/ subdir per system):
#   3. lyso     3aw6 -> 3aw7   P4(3)2(1)2
#   4. dhfr     1rx2 -> 1rx1   P2(1)2(1)2(1)
#   5. raddam   5kxk -> 5kxl   P4(3)2(1)2   subtract=bent
#   6. raddam   5kxk -> 5kxm   P4(3)2(1)2   subtract=bent
#   7. raddam   5kxk -> 5kxn   P4(3)2(1)2   subtract=bent
#   8. myoglobin 1mbo -> 1a6m  P2(1)
#   9. insulin  4fg3 -> 4e7u   H3           fill_fcalc=True
#  10. porin    altalign 3poq->3pou + refmac on the R32:R output
#
# Each test writes to test_results_<timestamp>/<name>.log and the script
# prints a PASS/FAIL summary at the end.  Run from the working area
# (the directory holding bendfinder.py, altalign.py, lyso/, dhfr/, ...).
#
# Usage:  ./run_all_tests.com

if (! $?CCP4) source /programs/ccp4-9/bin/ccp4.setup-csh

set ts     = `date +%Y%m%d_%H%M%S`
set logdir = test_results_$ts
mkdir -p $logdir
set summary = $logdir/SUMMARY.txt
echo "logs and per-test details: $logdir/"

@ npass = 0
@ nfail = 0
printf "%-30s %-7s %s\n" "test" "result" "metric" > $summary
printf "%-30s %-7s %s\n" "------------------------------" "------" "------" >> $summary

# ── 1. test_symm_all_sgs ─────────────────────────────────────────────────
echo ""
echo "==> [ 1/10] test_symm_all_sgs.py"
srun --job-name=test_symm ccp4-python test_symm_all_sgs.py >& $logdir/01_test_symm.log
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

# ── 2. magdoff Test 2 ────────────────────────────────────────────────────
echo ""
echo "==> [ 2/10] magdoff/test_magdoff.py"
( cd magdoff && srun --job-name=magdoff ccp4-python test_magdoff.py ) >& $logdir/02_magdoff.log
# The RMSD line is "Constrained   RMSD(CA) after fit: 0.0276 Å (91.8% ...";
# awk position 5 is the number.  Grep specifically for "after fit:" to avoid
# matching the later "Constrained (symm) : Riso=..." line.
set rcon = `grep "RMSD(CA) after fit:" $logdir/02_magdoff.log | grep "Constrained" | awk '{print $5}'`
set runc = `grep "RMSD(CA) after fit:" $logdir/02_magdoff.log | grep "Unconstrained" | awk '{print $5}'`
if ("$rcon" != "" && "$runc" != "") then
    echo "    PASS  constrained=$rcon  unconstrained=$runc"
    @ npass++
    printf "%-30s %-7s constrained=%s  unconstrained=%s\n" "magdoff" PASS "$rcon" "$runc" >> $summary
else
    echo "    FAIL  could not parse RMSDs from log"
    @ nfail++
    printf "%-30s %-7s parse_error\n" "magdoff" FAIL >> $summary
endif

# Reusable: run one fitreso_scan example.  Args via env vars (tcsh has no
# functions; using a foreach loop with positional vars per iteration).
#
# For each example: spec = "label sys mov_pdb ref_pdb mov_mtz ref_mtz subtract fill_fcalc tag"
# tag goes in the scan_dir name (so raddam variants don't collide).

set i = 2
foreach spec ( \
    "lyso        lyso        3aw6.pdb 3aw7.pdb 3aw6.mtz 3aw7.mtz ref  False ''" \
    "dhfr        dhfr        1rx2.pdb 1rx1.pdb 1rx2.mtz 1rx1.mtz ref  False ''" \
    "raddam_5kxl raddam      5kxk.pdb 5kxl.pdb 5kxk.mtz 5kxl.mtz bent False _5kxl" \
    "raddam_5kxm raddam      5kxk.pdb 5kxm.pdb 5kxk.mtz 5kxm.mtz bent False _5kxm" \
    "raddam_5kxn raddam      5kxk.pdb 5kxn.pdb 5kxk.mtz 5kxn.mtz bent False _5kxn" \
    "myoglobin   myoglobin   1mbo.pdb 1a6m.pdb 1mbo.mtz 1a6m.mtz ref  False ''" \
    "insulin     insulin     4fg3.pdb 4e7u.pdb 4fg3.mtz 4e7u.mtz ref  True  ''" \
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
    set scan = scan_test$tag
    @ i++
    set log = $logdir/`printf "%02d" $i`_$label.log

    echo ""
    echo "==> [`printf '%2d' $i`/10] $label  ($sys $movp -> $refp, fill_fcalc=$fill, subtract=$sub)"
    ( cd $sys && rm -rf $scan && srun --job-name=$label ccp4-python -c "import sys; sys.path.insert(0,'..'); from bendfinder import fitreso_scan; fitreso_scan(mov_pdb='$movp', ref_pdb='$refp', mov_mtz='$movm', ref_mtz='$refm', scan_dir='$scan', run_refinement_flag=True, refine_cycles=5, fill_fcalc=$fill, subtract='$sub')" ) >& $log

    set fr5 = `grep '^    fr5' $sys/$scan/scan_fitreso.log |& head -1`
    if ("$fr5" != "") then
        set rmsd  = `echo "$fr5" | awk '{print $2}'`
        set rbent = `echo "$fr5" | awk '{print $4}'`
        echo "    PASS  fr5  RMSD=$rmsd  Rbent=$rbent"
        @ npass++
        printf "%-30s %-7s fr5  RMSD=%s  Rbent=%s\n" "$label" PASS "$rmsd" "$rbent" >> $summary
    else
        echo "    FAIL  no fr5 row in scan_fitreso.log"
        @ nfail++
        printf "%-30s %-7s no_fr5_row\n" "$label" FAIL >> $summary
    endif
end

# ── 10. porin altalign + refmac on R32:R ─────────────────────────────────
@ i++
echo ""
echo "==> [`printf '%2d' $i`/10] porin  altalign 3poq->3pou + refmac on R32:R"
set plog = $logdir/`printf "%02d" $i`_porin_altalign.log
( cd porin && rm -rf altalign_test && mkdir altalign_test && \
    srun --job-name=porin_aa ccp4-python ../altalign.py 3poq.pdb 3pou.pdb \
        altalign_test/aa.pdb 3poq.mtz altalign_test/aa.mtz ) >& $plog

set residual = `grep 'honest discrete CA RMSD' $plog | head -1 | awk '{print $7}'`
if (-e porin/altalign_test/aa_R32R.pdb && -e porin/altalign_test/aa_R32R.mtz) then
    set rlog = $logdir/`printf "%02d" $i`_porin_refmac.log
    srun --job-name=porin_refmac ccp4-python -c "import sys; sys.path.insert(0,'.'); from bendfinder import run_refinement; run_refinement('porin/altalign_test/aa_R32R.pdb', 'porin/altalign_test/aa_R32R.mtz', outdir='porin/altalign_test/refine', n_cycles=5, fill_fcalc=True)" >& $rlog

    # refmac prints BOTH "Overall weighted R factor" and "Overall weighted R2
    # factor"; the trailing-space match keeps only the R-not-R2 line.
    set rf = `grep 'Overall weighted R factor' porin/altalign_test/refine/*refmac.log |& tail -1 | awk '{print $6}'`
    if ("$rf" != "") then
        echo "    PASS  altalign residual=$residual A   refmac R=$rf"
        @ npass++
        printf "%-30s %-7s altalign=%sA  refmac=R%s\n" "porin (altalign+R32:R)" PASS "$residual" "$rf" >> $summary
    else
        echo "    FAIL  refmac did not produce R-factor"
        @ nfail++
        printf "%-30s %-7s no_refmac_R\n" "porin (altalign+R32:R)" FAIL >> $summary
    endif
else
    echo "    FAIL  altalign did not write R32R outputs"
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
