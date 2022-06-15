#! /bin/tcsh -f
#
#    origins.com                    - James Holton 9-18-19
#
#    script for translating one PDB to a number of alternative origins
#       and indexing conventions
#    and checking if the symmetry-expanded atoms are close to those of
#    another PDB.
#
#    Atoms with the same name and resiude number are compared
#
set awk = awk
$awk 'BEGIN{print}' >& /dev/null
if($status) set awk = gawk
alias awk $awk


set SG        = ""
set right_pdb = ""
set wrong_pdb = ""
set outfile = neworigin.pdb

set max_opeaks = 30
set min_opeak = 2.5
set skip_noid = 1

mkdir -p ${CCP4_SCR} >&! /dev/null
set tempfile = ${CCP4_SCR}/origins_temp$$
if($?debug) set tempfile = ./origins_temp
################################################################################
goto Setup
Help:
cat << EOF

usage: $0 right_origin.pdb wrong_origin.pdb $SG [correlate] [nochains]

where: right_origin.pdb    is relative to your "desired" origin
       wrong_origin.pdb    is relative to another origin
       $SG        is your space group

wrong_origin.pdb will be moved to every possible origin in ${SG}
and then checked to see if the moved, and symmetry-expanded atoms
line up with the ones in right_origin.pdb.

The results of the origin choice that give the best agreement
to the reference pdb will be copied to $outfile

Using the word "correlate" on the command line signals the program
to ignore atom names in the comparision and instead use the correlation
coefficient of calculated electron density maps as the "similarity score"

By default, the program breaks up the input PDB files into any "chains"
specified therein.  You can turn this off by using the word "nochains" on 
the command line.

Using the word "otherhand" will also check for hand-inversion.

Using the word "altindex" will check for alternative indexing conventions.

Using the word "noorigins" will turn off the origin shift and just search
for symmetry-allowed alignments.

Using the word "fast" will skip atom-by-atom label re-assignments and stop
if a match is better than rmsd=1 or CC=0.8.

EOF
exit 9
Return_from_Setup:
################################################################################

set SCORE = "rmsd"
if($?CORRELATE) set SCORE = " CC "

# get unit cell from either pdb
set CELL = `awk '/^CRYST/{print $2,$3,$4,$5,$6,$7}' $right_pdb $wrong_pdb | head -1`
if("$CELL" == "") then
    echo "ERROR: no CRYST card in $right_pdb"
    goto Help
endif
# just in case these are different?
set right_cell = `awk '/^CRYST/{print $2,$3,$4,$5,$6,$7}' $right_pdb | head -1`
set wrong_cell = `awk '/^CRYST/{print $2,$3,$4,$5,$6,$7}' $wrong_pdb | head -1`
if("$right_cell" == "") set right_cell = "$CELL"
if("$wrong_cell" == "") set wrong_cell = "$right_cell"
# make the wrong_cell the default, so we are moving in allowed space
set CELL = ( $wrong_cell )

# try to get space group?
if("$SG" == "") then
    set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $right_pdb $wrong_pdb | head -1`
    set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
endif
if("$SG" == "") then
    # hmm, throw an error or go on? ...
    set SG = P1
endif
set right_latt = `echo $SG | awk '{print substr($0,1,1)}'`

# check for H3/R3 weirdness
set test = `echo $CELL | awk '{print ( $4+0==90 && $5+0==90 && $6+0==120 )}'`
if("$right_latt" == "R" && "$test" == "1" ) then
    # probably should use hexagonal system
    set SG = `echo $SG | awk '{gsub("R","H");print}'`
endif
if("$right_latt" == "H" && "$test" == "0" ) then
    # probably should use rhombohedral system
    set SG = `echo $SG | awk '{gsub("H","R");print}'`
endif

# check if wrong pdb has different info
set wrong_pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $right_pdb $wrong_pdb | head -1`
set wrong_SG = `awk -v pdbSG="$wrong_pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
if("$wrong_SG" == "") set wrong_SG = $SG
set wrong_latt = `echo $wrong_SG | awk '{print substr($0,1,1)}'`
set test = `echo $wrong_cell | awk '{print ( $4+0==90 && $5+0==90 && $6+0==120 )}'`

if("$wrong_latt" == "R" && "$test" == "1" ) then
    # probably should use hexagonal system
    set wrong_SG = `echo $wrong_SG | awk '{gsub("R","H");print}'`
endif
if("$wrong_latt" == "H" && "$test" == "0" ) then
    # probably should use rhombohedral system
    set wrong_SG = `echo $wrong_SG | awk '{gsub("H","R");print}'`
endif
set wrong_latt = `echo $wrong_SG | awk '{print substr($0,1,1)}'`



# space group for correlation calculations
set CCSG = $SG

# preemtive cleanup
rm -f ${tempfile}_right[0-9][0-9][0-9].pdb >& /dev/null
rm -f ${tempfile}_wrong[0-9][0-9][0-9].pdb >& /dev/null


set reindexings = +X,+Y,+Z

# need to use the "lattice cell" for re-indexing pdb file
set lattCELL = ( $wrong_cell )
echo "CELL $wrong_cell" | pdbset xyzin $wrong_pdb xyzout ${tempfile}.pdb > /dev/null
cat ${tempfile}.pdb |\
awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00 80.00"} \
    /^ANIS/{next} {print}' |\
cat >! ${tempfile}reindexme.pdb

if(! $?ALTINDEX) goto otherhand

echo "finding all possible re-indexing operators..."
echo "within $wrong_cell on latt= $wrong_latt"
# find all possible reindexing operations
othercell << EOF | tee ${tempfile}othercell.log > /dev/null
$wrong_cell $wrong_latt
EOF
#set reciprocal_symops = `awk '/close to target:/{getline;while(/]$/){gsub("[][]","");print $NF;getline}}' ${tempfile}othercell.log | sort -u`

set lattCELL = `awk '/Lattice unit cell after reindexing/{print $6,$7,$8,$9,$10,$11}' ${tempfile}othercell.log`
echo "lattice cell for re-indexing coordinates: $lattCELL"

# need to use the "lattice cell" for re-indexing pdb file
# dont forget to change it back
echo "CELL $lattCELL" | pdbset xyzin $wrong_pdb xyzout ${tempfile}lattcell.pdb > /dev/null
cat ${tempfile}lattcell.pdb |\
awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00 80.00"} \
    /^ANIS/{next} {print}' |\
cat >! ${tempfile}reindexme.pdb

# get unit cell shifts associated with each re-indexing?

# gather all symmetry operators from all possible space groups for this cell
cat ${tempfile}othercell.log |\
awk '$1~/^</ && $NF~/>$/' | sort -u |\
awk '{split($0,w,"[<>]");print "SPACEGROUP",w[2]}' |\
cat - ${CLIBD}/symop.lib |\
awk '/^SPACEGROUP/{pdbSG=substr($0,12);\
            #print "SG:",pdbSG;\
            ++seen[pdbSG];next}\
        /[XYZ]/ && ! /[PCIFRH]/ && p==1 {print $1 $2 $3 $4}\
    $5 ~ /^PG/{p=0} {split($0,w,"\047");pdbSG=w[2];if(seen[pdbSG])p=1}' |\
sort -u >! ${tempfile}all_symops

# and eliminate ones that are already in the space group to be used
cat ${CLIBD}/symop.lib |\
awk -v SG=$CCSG '/[XYZ]/ && ! /[PCIFRH]/ && p==1 {print $1 $2 $3 $4}\
    $5 ~ /^PG/{p=0} $4 == SG{p=1}' |\
awk '{print "GOT",$0}' >! ${tempfile}SG_symops

echo "+X,+Y,+Z" >! ${tempfile}reindexings.txt
cat ${tempfile}SG_symops ${tempfile}all_symops |\
awk '/^GOT/{++got[substr($0,4)];next}\
    ! got[$0]{print}' |\
awk '{gsub(" ","");print}' |\
awk -F "," '{for(i=1;i<=NF;++i){print $i};print ""}' |\
awk '/^[^-]/{$0="+"$0} {print}' |\
awk 'NF==0{print substr(op,2);op=""} NF>0{op=op"," $0}' |\
cat >> ${tempfile}reindexings.txt

# get list of all alternative indexings to try
set reindexings = `cat ${tempfile}reindexings.txt`
echo "$#reindexings alternative indexing operators."

if(! $?debug) then
    rm -f ${tempfile}symops
    rm -f ${tempfile}all_symops
    rm -f ${tempfile}SG_symops
    rm -f ${tempfile}othercell.log
endif

# filter out equivalent operators?
echo "filtering out redundant re-indexing operators"

# do a coarser map grid sampling than normal (~ 3A resolution)
if(! $?reso) set reso = 3
set fakeCELL = `echo $lattCELL $reso | awk '{print $1/$NF/2,$2/$NF/2,$3/$NF/2,$4,$5,$6}'`
set BADD     = `echo $reso | awk '{printf "%d", 79*($1/3)^2}'`
echo "CELL $fakeCELL" | pdbset xyzin $wrong_pdb xyzout ${tempfile}.pdb > /dev/null

# make an all-carbon version of this model
cat ${tempfile}.pdb |\
awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00  5.00"}\
     /^ANIS/{next} {print}' |\
cat >! ${tempfile}mapme.pdb
sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}right.map << EOF-sfall > /dev/null
MODE ATMMAP
CELL $fakeCELL
SYMM $CCSG
BADD $BADD
EOF-sfall
# recover coarse grid spacing 
set coarseGRID = `echo "GO" | mapdump MAPIN ${tempfile}right.map | awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}'`
rm -f ${tempfile}.pdb >& /dev/null

# now that we have a grid, make a safe version of whole file

set n = 0
foreach reindexing ( $reindexings )
    @ n = ( $n + 1 )
    #echo "applying $reindexing to $wrong_pdb"
    pdbset xyzin ${tempfile}reindexme.pdb xyzout ${tempfile}reindexed.pdb << EOF > /dev/null
    symgen $reindexing
    BFACTOR 80
EOF

    sfall xyzin ${tempfile}reindexed.pdb mapout ${tempfile}test_${n}.map << EOF > /dev/null
    mode atmmap
    resolution 3
    GRID $coarseGRID
    symm $SG
EOF
end
echo "maps done"

# now look for maps that are prefectly correlated, these indicate one of the ops is redundant
set redundant = ""
rm -f CC_vs_ij.log
foreach i ( `seq 1 $n` )
  @ k = ( $i + 1 )
  foreach j ( `seq $k $n` )
    set test = `echo $i $j $redundant | awk '{for(k=3;k<=NF;++k){if($1==$k || $2==$k){print 1;exit}}}'`
    if("$test" == "1") continue
    echo "correlate section" |\
    overlapmap mapin1 ${tempfile}test_${i}.map \
            mapin2 ${tempfile}test_${j}.map mapout /dev/null |\
    awk '/Total correlation/{print $NF}' >! ${tempfile}correlation
    set CC = `awk '{print $1}' ${tempfile}correlation`
    #echo "$CC $i $j" | tee -a CC_vs_ij.log
    set test = `echo $CC | awk '{print ( $1 > 0.95 )}'`
    if($test) then
       # the jth operator is redundant with teh ith operator
       set redundant = ( $redundant $j )
    endif
  end
end
foreach dupe ( $redundant )
    set reindexings[$dupe] = ""
end
set reindexings = ( $reindexings )
echo "$#reindexings re-indexing operators remain:"
if(! $?debug) rm -f ${tempfile}correlation >& /dev/null
if(! $?debug) rm -f ${tempfile}test_*.map >& /dev/null
if(! $?debug) rm -f CC_vs_ij.log >& /dev/null


otherhand:
if($?OTHERHAND) then
    echo "adding other hand..."
    echo " $reindexings " |\
    awk '{for(i=1;i<=NF;++i)print $i}' |\
    awk '{print;\
          gsub("+","q");gsub("-","+");gsub("q","-");\
          print}' |\
    cat >! ${tempfile}newops.txt
    set reindexings = `sort -u ${tempfile}newops.txt`
    echo "$#reindexings re-indexing operators remain:"
    rm -f ${tempfile}newops.txt
endif


# now re-discover direct and reciprocal space operations
# first, make a reference mtz file for wrong_pdb expanded to P1
echo CELL $wrong_cell |\
  pdbset xyzin ${tempfile}reindexme.pdb xyzout ${tempfile}ref.pdb > /dev/null
sfall xyzin ${tempfile}ref.pdb hklout ${tempfile}sfalled.mtz << EOF > /dev/null
mode sfcalc xyzin
resolution 3
symm $SG
EOF
cad hklin1 ${tempfile}sfalled.mtz hklout ${tempfile}P1.mtz << EOF > /dev/null
labin file 1 all
outlim P1
EOF
cad hklin1 ${tempfile}P1.mtz hklout ${tempfile}reindexme.mtz << EOF > /dev/null
labin file 1 all
symm 1
EOF
cp -p ${tempfile}sfalled.mtz ${tempfile}reindexme.mtz

# now loop over all real-space reindexing operators and see what they do in reciprocal space
rm -f ${tempfile}xyz_hkl_cell.log
foreach reindexing ( $reindexings )
    pdbset xyzin ${tempfile}reindexme.pdb xyzout ${tempfile}.pdb << EOF > /dev/null
    symgen $reindexing
    BFACTOR 80
    SPACE $SG
EOF
    echo CELL $wrong_cell |\
    pdbset xyzin ${tempfile}.pdb xyzout ${tempfile}reindexed.pdb > /dev/null
    echo "tolerance 10" |\
    pointless xyzin ${tempfile}reindexed.pdb \
              hklin ${tempfile}reindexme.mtz \
              hklout ${tempfile}reindexed.mtz |\
    tee ${tempfile}pointless.log |\
    awk '/Alternative /,NF==0' |\
    awk '{print $4+0,$2}' | sort -gr |\
    awk '{gsub(/[\]\[]/,"");print $2;exit}' |\
    sort -u >! ${tempfile}reindexing.txt
    set hkl = `cat ${tempfile}reindexing.txt`
    set cell = `echo head | mtzdump hklin ${tempfile}reindexed.mtz | awk '/Cell Dimensions/{getline;getline;print}'`
    echo "$reindexing   $hkl   $cell" |\
     awk '{printf("%-20s %-10s %s %s %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' |\
     tee -a ${tempfile}xyz_hkl_cell.log
end


if (! $?BYFILE) then
    # break up "right" file into its respective chains
    cat $right_pdb |\
    awk '/^ATOM|^HETAT/{resnum=substr($0, 23, 4)+0;segid = substr($0, 22, 1); \
         if(last_segid=="")last_segid=segid;\
             if(resnum<last_resnum || (segid!=last_segid)) print "BREAK";\
                last_resnum=resnum;last_segid=segid} {print}' |\
    awk -v tempfile=${tempfile}_right 'BEGIN{chain="001"} /^ATOM|^HETAT/{segid = substr($0, 22, 1)}\
     /^BREAK/{print "REMARK chain", segid > tempfile chain ".pdb";\
              chain=sprintf("%03d", chain+1)}\
      /^ATOM|^HETATM/{print "ATOM  "substr($0,7,66) > tempfile chain ".pdb"}\
          END{print "REMARK chain", segid > tempfile chain ".pdb"}'
    # chains are now separated into files
    # ${tempfile}_right###.pdb

    # break up "wrong" file into its respective chains
    cat $wrong_pdb |\
    awk '/^ATOM|^HETAT/{resnum=substr($0, 23, 4)+0;segid = substr($0, 22, 1); \
         if(last_segid=="")last_segid=segid;\
             if(resnum<last_resnum || (segid!=last_segid)) print "BREAK";\
                last_resnum=resnum;last_segid=segid} {print}' |\
    awk -v tempfile=${tempfile}_wrong 'BEGIN{chain="001"} /^ATOM|^HETAT/{segid = substr($0, 22, 1)}\
     /^BREAK/{print "REMARK chain", segid > tempfile chain ".pdb";\
              chain=sprintf("%03d", chain+1)}\
      /^ATOM|^HETAT/{print "ATOM  "substr($0,7,66) > tempfile chain ".pdb"}\
          END{print "REMARK chain", segid > tempfile chain ".pdb"}'
    # chains are now separated into files
    # ${tempfile}_wrong###.pdb
endif

# automatically switch to "nochains" mode if there is only one chain per PDB?
set test = `ls -1 ${tempfile}_right???.pdb ${tempfile}_wrong???.pdb | wc -l`
if("$test" == "2") set BYFILE

if($?BYFILE) then
    echo "REMARK CHAIN _" >! ${tempfile}_right001.pdb
    cat $right_pdb |\
    awk '/^ATOM|^HETAT/{print "ATOM  "substr($0,7,66)}' |\
    cat >> ${tempfile}_right001.pdb
    echo "REMARK CHAIN _" >> ${tempfile}_right001.pdb

    echo "REMARK CHAIN _" >! ${tempfile}_wrong001.pdb
    cat $wrong_pdb |\
    awk '/^ATOM|^HETAT/{print "ATOM  "substr($0,7,66)}' |\
    cat >> ${tempfile}_wrong001.pdb
    echo "REMARK CHAIN _" >> ${tempfile}_wrong001.pdb
    
endif

# add the proper headers
foreach file ( ${tempfile}_right[0-9][0-9][0-9].pdb ${tempfile}_wrong[0-9][0-9][0-9].pdb )

    echo "END" >> $file

    # make sure the cell gets encoded
    set CELL = "$right_cell"
    if("$file" =~ "*_wrong*") set CELL = "$wrong_cell"

    # calculate the center of mass
    pdbset xyzin $file xyzout ${tempfile}.pdb << EOF >! ${tempfile}.log
    CELL $CELL
    CHAIN " "
    COM
EOF
    egrep "^REMARK" $file >! ${tempfile}new.pdb
    egrep -v "REMARK" ${tempfile}.pdb >> ${tempfile}new.pdb
    mv ${tempfile}new.pdb ${tempfile}.pdb
    
    # make a pdb file of an atom at the COM of this file
    awk '/^CRYST/ || /^SCALE/' ${tempfile}.pdb >! ${tempfile}COM.pdb
    cat ${tempfile}.log |\
    awk '$1=="Center"{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\n", $4, $5, $6;exit}' |\
    cat >> ${tempfile}COM.pdb
    echo "END" >> ${tempfile}COM.pdb
    
    # get fractional version of the COM too
    coordconv XYZIN ${tempfile}COM.pdb XYZOUT ${tempfile}.xyz << EOF >& /dev/null
    CELL $CELL
    INPUT PDB
    OUTPUT FRAC
EOF

    # put this center as a remark in the PDB file
    cat ${tempfile}.log |\
    awk '$1=="Center"{print "COM", $4, $5, $6}' |\
    cat - ${tempfile}.xyz |\
    awk '{printf "REMARK " $0; getline; print "  ", $2, $3, $4}' |\
    cat >! $file
    cat ${tempfile}.pdb |\
    awk '/^REMARK/ && $NF=="chain"{print "REMARK chain _";next}\
          /^REMARK/{print}' >> $file
    cat ${tempfile}.pdb |\
    awk '! /REMARK/{print substr($0,1,66)}' >> $file
    rm -f ${tempfile}.pdb ${tempfile}COM.pdb ${tempfile}.xyz ${tempfile}.log >& /dev/null

    # now "$file" has its center of mass in its header
end

# display stats
echo -n "reference chains: "
foreach file ( ${tempfile}_right[0-9][0-9][0-9].pdb )
    awk '/^REMARK chain/{printf "%s ", $3;exit}' $file
end
echo "    ( "`basename $right_pdb`" )"
echo -n "  subject chains: "
foreach file ( ${tempfile}_wrong[0-9][0-9][0-9].pdb )
    awk '/^REMARK chain/{printf "%s ", $3;exit}' $file
end
echo "    ( "`basename $wrong_pdb`" )"
echo "CELL $CELL"
echo "SG   $SG"
echo ""







################################
# get tranformations from a space group

# get symmetry operations from space group
cat ${CLIBD}/symop.lib |\
awk -v SG=$SG '/[XYZ]/ && ! /[PCIFRH]/ && p==1 {print $1 $2 $3 $4}\
$5 ~ /^PG/{p=0} $4 == SG{p=1}' |\
cat >! ${tempfile}symops
set symops = `cat ${tempfile}symops`

# get origin-shift operators from a space group (at the bottom of this script)
set table = $0
cat $table |\
awk '/TABLE OF ALLOWED ORIGIN SHIFTS/,/END OF THE TABLE/' |\
awk -v SG=$SG '/^[PCIFRH]/ && $0 ~ SG" "{getline; while(NF>0){\
    print $1 "," $2 "," $3;getline};exit}' |\
cat >! ${tempfile}origins
# format: dX,dY,dZ
set origins = `cat ${tempfile}origins`
#set origins = `awk 'BEGIN{for(y=0;y<1;y+=0.01){print "0,"y",0";print "1/2,"y",0";print "0,"y",1/2";print "1/2,"y",1/2"}}'`

if($?NO_ORIGINS || "$origins" == "") then
    set origins = "0,0,0"
    set best_reindexing_origin = "+X,+Y,+Z_0,0,0"
    set CCSG = 1
endif

rm -f ${tempfile}origins
rm -f ${tempfile}symops

if( $?NO_DECONV ) then
    echo "skipping deconvolution step, using origins: $origins"
    set reindexings = +X,+Y,+Z
    set best_reindexing_origin = "+X,+Y,+Z_0,0,0"
    set reindexing_origins = `echo $origins | awk '{for(i=1;i<=NF;++i){print "+X,+Y,+Z_"$i}}'`
    echo "$reindexings  h,k,l $CELL" |\
     awk '{printf("%-20s %-10s %s %s %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' |\
    cat >! ${tempfile}xyz_hkl_cell.log
    goto skip_deconv
endif


deconvolute:
# preemtive:
rm -f ${tempfile}all_reindexing_origins.txt

# use deconvolution to find optimal shift
foreach reindexing ( $reindexings )
    echo "applying $reindexing to $wrong_pdb"
    set CELL = `awk -v key="$reindexing" '$1==key{print $3,$4,$5,$6,$7,$8;exit}' ${tempfile}xyz_hkl_cell.log`
    set reindexing_hkl = `awk -v key="$reindexing" '$1==key{print $2;exit}' ${tempfile}xyz_hkl_cell.log`
    echo "new cell: $CELL"

    # apply the re-indexing operation, which is done in lattice cell
    pdbset xyzin ${wrong_pdb} xyzout ${tempfile}lattcell.pdb << EOF > /dev/null
    CELL $lattCELL
EOF
    pdbset xyzin ${tempfile}lattcell.pdb xyzout ${tempfile}symmed.pdb << EOF > /dev/null
    symgen $reindexing
EOF
    # and then put back the relevant cell and SG in the header
    pdbset xyzin ${tempfile}symmed.pdb xyzout ${tempfile}reindexed.pdb << EOF > /dev/null
    CELL $CELL
    space $SG
EOF


    echo "deconvoluting maps..."

    # do a coarser map grid sampling than normal, ~ 2A resolution
    set reso = 2

    # make a "mask" of all possible origin shifts
    echo -n "" >! ${tempfile}xyz.txt
    if("$SG" == "P1") then
        echo "0 0 0" >> ${tempfile}xyz.txt
    else
        foreach origin ( $origins )
            # break up the origin string
            set Xf = `echo "$origin" | awk -F "," '{print $1}'`
            set Yf = `echo "$origin" | awk -F "," '{print $2}'`
            set Zf = `echo "$origin" | awk -F "," '{print $3}'`
            set X = `echo "$Xf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`
            set Y = `echo "$Yf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`
            set Z = `echo "$Zf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`
        
            echo "$X $Y $Z $Xf $Yf $Zf" |\
            awk '{print $1+0,$2+0,$3+0}\
             $4=="x"{for(x=0;x<1;x+=0.001)print x,$2+0,$3+0}\
             $5=="y"{for(y=0;y<1;y+=0.001)print $1+0,y,$3+0}\
             $6=="z"{for(z=0;z<1;z+=0.001)print $1+0,$2+0,z}\
             $4=="x=y=z"{for(x=0;x<2;x+=0.001)print x,x,x}' |\
            cat >> ${tempfile}xyz.txt
        end
    endif
    # make a frac coord file from all possible origin positions
    cat ${tempfile}xyz.txt |\
    awk 'NF>2{++n; printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%2s%3s%3s %1s\n", \
               n, $1, $2, $3, 80, 1, "38", n, "H", "", "IUM", " "}' |\
    cat >! ${tempfile}.frac
    coordconv XYZIN ${tempfile}.frac \
          XYZOUT ${tempfile}originmask.pdb << EOF-conv >& /dev/null
CELL $CELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv
    sfall xyzin ${tempfile}originmask.pdb hklout ${tempfile}originmask.mtz << EOF-sfall > /dev/null
MODE sfcalc xyzin
CELL $CELL
SYMM 1
RESO $reso
EOF-sfall
    fft hklin ${tempfile}originmask.mtz mapout ${tempfile}originmask.map << EOF >! ${tempfile}.log
    labin F1=FC PHI=PHIC
    scale F1 1 80
    reso $reso
    symm P1
EOF
    # turn it into a binary mask
    set scale = `awk '/Maximum density/ && $NF>0{print 1/$NF}' ${tempfile}.log`
    rm -f ${tempfile}.log
    if("$scale" == "") set scale = 1
    echo "scale factor $scale 0" |\
      mapmask mapin ${tempfile}originmask.map mapout ${tempfile}new.map > /dev/null
    mv ${tempfile}new.map ${tempfile}originmask.map
    echo "scale factor 0 1" |\
      mapmask mapin ${tempfile}originmask.map mapout ${tempfile}one.map > /dev/null

#    echo "mask cut 0" |\
#    mapmask mapin ${tempfile}originmask.map mskout ${tempfile}originmask.msk > /dev/null
#    echo "maps mult" |\
#    mapmask mapin ${tempfile}one.map mskin ${tempfile}originmask.msk \
#     mapout ${tempfile}originmask.map > /dev/null
    if("$SG" == "P1") then
        # P1 origin is good everywhere, so no mask
        cp ${tempfile}one.map ${tempfile}originmask.map > /dev/null
    endif
    rm -f ${tempfile}one.map ${tempfile}originmask.mtz
    rm -f ${tempfile}originmask.pdb ${tempfile}.frac ${tempfile}xyz.txt


    # make an all-carbon version of right model
    echo "CELL $right_cell" |\
       pdbset xyzin $right_pdb xyzout ${tempfile}.pdb >! ${tempfile}.log
    if($status) then
        set BAD = "cannot set unit cell of right_pdb"
        goto exit
    endif
    egrep "^CRYST1|^ATOM|^HETAT" ${tempfile}.pdb |\
    awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00 80.00"} {print}' |\
    cat >! ${tempfile}mapme.pdb
    sfall xyzin ${tempfile}mapme.pdb hklout ${tempfile}right.mtz << EOF-sfall >! ${tempfile}.log
    MODE sfcalc xyzin
    CELL $right_cell
    SYMM $SG
    RESO $reso
EOF-sfall
    if($status) then
        set BAD = "cannot sfall on right_pdb"
        goto exit
    endif
    # make an all-carbon version of wrong, re-indexed model
    echo "CELL $CELL" |\
      pdbset xyzin ${tempfile}reindexed.pdb xyzout ${tempfile}.pdb > /dev/null
    egrep "^CRYST1|^ATOM|^HETAT" ${tempfile}.pdb |\
    awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00 80.00"} {print}' |\
    cat >! ${tempfile}mapme.pdb
    sfall xyzin ${tempfile}mapme.pdb hklout ${tempfile}wrong.mtz << EOF-sfall >! ${tempfile}.log
    MODE sfcalc xyzin
    CELL $CELL
    SYMM $SG
    RESO $reso
EOF-sfall
    if($status) then
        set BAD = "cannot sfall on wrong_pdb"
        goto exit
    endif

    rm -f ${tempfile}del.mtz
    sftools << EOF > /dev/null
read ${tempfile}right.mtz
read ${tempfile}wrong.mtz
set labels
Fright
PHIright
Fwrong
PHIwrong
calc ( COL Fq PHIdel ) = ( COL Fright PHIright ) ( COL Fwrong PHIwrong ) /
calc COL W = COL Fq
select COL Fq > 1
calc COL W = 1 COL Fq /
select all
calc F COL Fdel = COL W 0.5 **
write ${tempfile}del.mtz col Fdel PHIdel
y
stop
EOF
    fft hklin ${tempfile}del.mtz mapout ${tempfile}del.map << EOF >! ${tempfile}.log
    labin F1=Fdel PHI=PHIdel
    reso $reso
    symm P1
EOF
    # make sure that we define "sigma" for the unmasked map
    echo "scale sigma 1 0" |\
    mapmask mapin ${tempfile}del.map mapout ${tempfile}zscore.map  > /dev/null
    mapmask mapin1 ${tempfile}zscore.map mapin2 ${tempfile}originmask.map \
     mapout ${tempfile}pickme.map << EOF >! ${tempfile}.log
mode mapin1 mapin2
maps mult
EOF
    peakmax mapin ${tempfile}pickme.map xyzfrc ${tempfile}peak.txt << EOF >! ${tempfile}.log
    output frac
#    threshold 3
#    numpeaks $max_opeaks
EOF
    # this scale compensates for the mask, restoring "sigma" units
    set scale = `awk '/peaks higher than the threshold/{print $(NF-1)/$9}' ${tempfile}.log`
#    set scale = `awk '/Rms deviation from mean density/{print $NF}' ${tempfile}.log`
#    echo | mapdump mapin ${tempfile}zscore.map | grep density
#    echo | mapdump mapin ${tempfile}pickme.map | grep density
#    head ${tempfile}peak.txt
#    set scale = 1
    echo "fractional origin shift   Z-score"
    cat ${tempfile}peak.txt |\
    awk -v scale=$scale '/ATOM|^HETAT/{printf "%7.4f %7.4f %7.4f   %s\n", $3,$4,$5,$6/scale}' |\
    awk -v min_opeak=$min_opeak '$4>min_opeak' |\
    awk 'NR==1{max=$4;print;next} \
       max>=9 && $4>max/3 || max<9' |\
     head -n $max_opeaks |\
    tee ${tempfile}likely_origins.txt
    wc -l ${tempfile}likely_origins.txt | awk '{print $1,"origins added."}'

    # gather stats in case we need them
    if(! $?maxpeak) set maxpeak = 0
    set thismax = `awk -v scale=$scale '{print $6/scale}' ${tempfile}peak.txt |& sort -gr |& head -n 1`
    set maxpeak = `echo $thismax $maxpeak | awk '$2>$1{$1=$2} {print $1}'`

    # clean up a bit
    rm -f ${tempfile}right.mtz ${tempfile}wrong.mtz
    rm -f ${tempfile}del.mtz ${tempfile}del.map
    rm -f ${tempfile}originmask.map ${tempfile}pickme.map ${tempfile}zscore.map
    rm -f ${tempfile}peak.txt ${tempfile}.log

    # now integrate these shifts with the "official" origin list
    echo "ORIGINS $origins" |\
    cat - ${tempfile}likely_origins.txt |\
    awk '/^ORIGINS/{\
      for(i=2;i<=NF;++i){\
        ++n;\
        split($i,xyz,",");\
        for(k in xyz){\
          if(xyz[k]~/\//){\
            split(xyz[k],w,"/");\
            xyz[k]=w[1]/w[2];\
          }\
        }\
        x0[n]=xyz[1];y0[n]=xyz[2];z0[n]=xyz[3];\
      }\
    }\
    NF==4{x=$1;y=$2;z=$3;score=$4;\
      minfd=999;\
      for(i in x0){\
        xp=x0[i];yp=y0[i];zp=z0[i];\
        if(xp=="x")xp=x;\
        if(yp=="y")yp=y;\
        if(zp=="z")zp=z;\
        if(xp=="x=y=z")xp=yp=zp=x;\
        origin[i]= xp","yp","zp;\
        fdx=sqrt((x-xp)^2);if(fdx>0.9)fdx-=1;\
        fdy=sqrt((y-yp)^2);if(fdy>0.9)fdy-=1;\
        fdz=sqrt((z-zp)^2);if(fdz>0.9)fdz-=1;\
        fd=sqrt(fdx^2+fdy^2+fdz^2);\
        if(fd<minfd){minfd=fd;best=i};\
      };\
      print score,1-minfd,origin[best], x,y,z;\
    }' |\
    awk -v op="$reindexing" '{print $0,op}' |\
    tee -a ${tempfile}all_reindexing_origins.txt > /dev/null
    # format: score rel_height x,y,z  x y z  reindexing 

    rm -f ${tempfile}likely_origins.txt

    set goodenough = `awk '$1>50' ${tempfile}all_reindexing_origins.txt | wc -l`
    if(("$goodenough" == "1") && ($?SPEEDUP)) then
        echo "thats good enough..."
        break
    endif
end

set test = `cat ${tempfile}all_reindexing_origins.txt | wc -l`
echo "$test possible origins to explore."
if("$test" == "" || "$test" == "0") then
    set test = `echo $min_opeak $maxpeak | awk '{print ($1>$2)}'`
    if($test) then
        set min_opeak = `echo $maxpeak | awk '{print $1*0.9}'`
        echo "re-setting min_opeak to $min_opeak "
        goto deconvolute
    endif
    set test = `echo $min_opeak 1 | awk '{print ($1>$2)}'`
    if($test) then
        echo "re-setting min_opeak to 1 "
        set min_opeak = 1
        goto deconvolute
    endif
    set BAD = "no significant origin peaks"
    goto exit
endif

sort -gr ${tempfile}all_reindexing_origins.txt >! ${tempfile}neworigins.txt
set reindexing_origins = `awk '{ro=$7"_"$3} ! seen[ro]{print ro} {++seen[ro]}' ${tempfile}neworigins.txt`
set best_reindexing = `awk '{print $7;exit}' ${tempfile}neworigins.txt`
set CELL = `awk -v key="$best_reindexing" '$1==key{print $3,$4,$5,$6,$7,$8;exit}' ${tempfile}xyz_hkl_cell.log`
set reindexing_hkl = `awk -v key="$best_reindexing" '$1==key{print $2;exit}' ${tempfile}xyz_hkl_cell.log`

#rm -f ${tempfile}all_reindexing_origins.txt ${tempfile}neworigins.txt

#set reindexing_origins = `seq -1 0.05 2 | awk '{print "+X,+Y,+Z_0,"$1",0"}'`
#echo "GOTHERE: $reindexing_origins"

skip_deconv:
########################################
# now run through all origins, and chain pairings
echo "chain                           origin "
echo "r    s      reindexing          dX  dY  dZ  symop                $SCORE"

again:
echo -n "" >! ${tempfile}scores
foreach right_chain ( ${tempfile}_right[0-9][0-9][0-9].pdb )

# get fractional coordinate limits that cover the "right" chain
coordconv xyzin $right_chain xyzout ${tempfile}.xyz << EOF >> /dev/null
INPUT PDB
OUTPUT FRAC
END
EOF
cat ${tempfile}.xyz |\
awk -v del=0.5 'BEGIN{xmin=ymin=zmin=99999999} \
$2<xmin{xmin=$2} $2>xmax{xmax=$2}\
$3<ymin{ymin=$3} $3>ymax{ymax=$3}\
$4<zmin{zmin=$4} $4>zmax{zmax=$4}\
END{print xmin-del, xmax+del, ymin-del, ymax+del, zmin-del, zmax+del}' |\
cat >! ${tempfile}.xyzlim
set xyzlim = `cat ${tempfile}.xyzlim`
rm -f ${tempfile}.xyzlim
rm -f ${tempfile}.xyz


if($?CORRELATE) then
    # we want to do comparison with maps, not atoms

    if(! $?coarseGRID) then
        # do a coarser map grid sampling than normal (~ 3A resolution)
        if(! $?reso) set reso = 3
        set fakeCELL = `echo $CELL $reso | awk '{print $1/$NF/2,$2/$NF/2,$3/$NF/2,$4,$5,$6}'`
        set BADD     = `echo $reso | awk '{printf "%d", 79*($1/3)^2}'`
        echo "CELL $fakeCELL" | pdbset xyzin $right_chain xyzout ${tempfile}.pdb > /dev/null

        # make an all-carbon version of this model
        cat ${tempfile}.pdb |\
        awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00  5.00"} {print}' |\
        cat >! ${tempfile}mapme.pdb
        sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}right.map << EOF-sfall > /dev/null
        MODE ATMMAP
        CELL $fakeCELL
        SYMM $CCSG
        BADD $BADD
EOF-sfall
        # recover coarse grid spacing 
        set coarseGRID = `echo "GO" | mapdump MAPIN ${tempfile}right.map | awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}'`
        rm -f ${tempfile}.pdb >& /dev/null
    endif
    set GRID = "$coarseGRID"

    # re-calculate the "right" map with reduced grid
    cat $right_chain |\
    awk '/^ATOM|^HETAT/{++n;$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00  5.00"} {print}' |\
    awk '/^ATOM|^HETAT/{++n;}\
     n==1{++n;\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\
         {print}' |\
    cat >! ${tempfile}mapme.pdb
    sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}right.map << EOF-sfall > /dev/null
    MODE ATMMAP
    CELL $right_cell
    SYMM $CCSG
    GRID $GRID
    BADD $BADD
EOF-sfall
    rm -f ${tempfile}mapme.pdb >& /dev/null
    # determine grid spacing to use for "wrong" map
    set GRID = `echo "GO" | mapdump MAPIN ${tempfile}right.map | awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}'`
#    set xyzlim = "0 1 0 1 0 1"
endif

unset SKIP_CHAINPAIR

# now loop over the wrong chains
foreach wrong_chain ( ${tempfile}_wrong[0-9][0-9][0-9].pdb )

# no need to continue if there are no identities
if($?SKIP_CHAINPAIR) then
    set test = `echo $SKIP_CHAINPAIR $right_chain $wrong_chain | awk '{print ( $1==$3 && $2==$4 )}'`
    unset SKIP_CHAINPAIR
    if($test && $skip_noid) continue
endif

# try every possible origin choice
foreach reindexing_origin ( $reindexing_origins )

if($?SKIP_CHAINPAIR) then
    set test = `echo $SKIP_CHAINPAIR $right_chain $wrong_chain | awk '{print ( $1==$3 && $2==$4 )}'`
    if($test && $skip_noid) break
endif

# break up the origin string
set reindexing = `echo $reindexing_origin | awk -F "_" '{print $1}'`
set CELL = `awk -v key="$reindexing" '$1==key{print $3,$4,$5,$6,$7,$8;exit}' ${tempfile}xyz_hkl_cell.log`
set reindexing_hkl = `awk -v key="$reindexing" '$1==key{print $2;exit}' ${tempfile}xyz_hkl_cell.log`

set origin = `echo $reindexing_origin | awk -F "_" '{print $2}'`
set Xf = `echo "$origin" | awk -F "," '{print $1}'`
set Yf = `echo "$origin" | awk -F "," '{print $2}'`
set Zf = `echo "$origin" | awk -F "," '{print $3}'`
set X = `echo "$Xf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`
set Y = `echo "$Yf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`
set Z = `echo "$Zf" | awk -F "/" 'NF==2{$1/=$2} {printf "%.10f",$1}'`

# apply this re-indexing to the "wrong" pdb
pdbset xyzin ${wrong_chain} xyzout ${tempfile}lattcell.pdb << EOF > /dev/null
CELL $lattCELL
EOF
# re-index the PDB in the lattice cell
pdbset xyzin ${tempfile}lattcell.pdb xyzout ${tempfile}symmed.pdb << EOF > /dev/null
symgen $reindexing
EOF
# and then put back the relevant cell and SG in the header
pdbset xyzin ${tempfile}symmed.pdb xyzout ${tempfile}wrong_reindexed.pdb << EOF > /dev/null
CELL $CELL
space $SG
EOF

# get the center of mass in orthogonal coordinates
pdbset xyzin ${tempfile}wrong_reindexed.pdb xyzout ${tempfile}.pdb << EOF >! ${tempfile}.log
COM
EOF
set wrong_COM = `awk '$1=="Center"{print $4,$5,$6}' ${tempfile}.log`

# make a pdb file of an atom at the COM of this file
awk '/^CRYST/ || /^SCALE/' ${tempfile}.pdb >! ${tempfile}COM.pdb
echo $wrong_COM |\
awk '{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\n", $1, $2, $3;exit}' |\
cat >> ${tempfile}COM.pdb
echo "END" >> ${tempfile}COM.pdb

# get fractional version of the COM too
coordconv XYZIN ${tempfile}COM.pdb XYZOUT ${tempfile}.xyz << EOF >& /dev/null
CELL $CELL
INPUT PDB
OUTPUT FRAC
EOF
set wrong_COM_frac = `awk '{print $2,$3,$4}' ${tempfile}.xyz`


# retrieve info from file headers
set right_COM = `awk '/^REMARK COM/{print $3, $4, $5;exit}' $right_chain`
#set wrong_COM = `awk '/^REMARK COM/{print $3, $4, $5;exit}' $wrong_chain`
set shift_COM = `echo "$wrong_COM $right_COM" | awk '{print $4-$1,$5-$2,$6-$3}'`
set right_COM_frac = `awk '/^REMARK COM/{print $6, $7, $8;exit}' $right_chain`
#set wrong_COM_frac = `awk '/^REMARK COM/{print $6, $7, $8;exit}' $wrong_chain`
set shift_COM_frac = `echo "$wrong_COM_frac $right_COM_frac" | awk '{print $4-$1,$5-$2,$6-$3}'`
set right_chain_ID = `awk '/^REMARK chain/{print $3}' $right_chain`
set wrong_chain_ID = `awk '/^REMARK chain/{print $3}' $wrong_chain`
if("$right_chain_ID" == "") set right_chain_ID = "_"
if("$wrong_chain_ID" == "") set wrong_chain_ID = "_"

# add translation along any polar axes
#set X = `echo "$opt_frac_shift $Xf $X" | awk '/x/{$NF=$1} {print $NF}'`
#set Y = `echo "$opt_frac_shift $Yf $Y" | awk '/y/{$NF=$2} {print $NF}'`
#set Z = `echo "$opt_frac_shift $Zf $Z" | awk '/z/{$NF=$3} {print $NF}'`

# calculate how to "level" the centers of the atom constellations
# NOTE: this is not a good idea if the chains are not rigid-body related
#set polar_slip = `echo "$opt_frac_shift $Xf $Yf $Zf" | awk '! /x/{$1=0} ! /y/{$2=0} ! /z/{$3=0} {print $1,$2,$3}'`

# now see which (if any) symops/cell translations will bring the 
# shifted COM of "wrong_chain" anywhere near the COM of the "right" chain

# apply current origin shift to "wrong" COM
set new_COM = `echo "$wrong_COM_frac $X $Y $Z" | awk '{print $1+$4, $2+$5, $3+$6}'`
gensym << EOF >! ${tempfile}.log
CELL $CELL
SYMM $SG
atom X $new_COM
XYZLIM $xyzlim
EOF

# only use symops that showed up in the gensym result
cat ${tempfile}.log |\
awk '/List of sites/,/atoms generated from/' |\
awk 'NF>10{print $NF, $2, $3, $4, $5, $6, $7}' |\
sort -u -k1,1 >! ${tempfile}symop_center
# format: symop xf yf zf X Y Z

if($?CORRELATE && ! $?best_reindexing_origin) then
    echo "1 $right_COM_frac $new_COM" >! ${tempfile}symop_center
endif

# loop over the symmetry operations
foreach line ( `awk '{print NR}' ${tempfile}symop_center` )

    # skip all other origins if a good one has been found
    if(("$good_origin" != "")&&("$good_origin" != "$X $Y $Z") && "$CCSG" != "1") then
        # echo "skipping: $good_origin   $X $Y $Z"
        continue
    endif
    # skip rest of "right_chain" if we have already found a match
    if("$good_right" == "$right_chain") then
        # echo "skipping: $good_right == $right_chain"
        continue
    endif
    if("$good_wrong" == "$wrong_chain") then
        # echo "skipping: $good_wrong == $wrong_chain"
        continue
    endif

    # retrieve symmetry operation
    set symop = `awk -v line=$line 'NR==line{print $1}' ${tempfile}symop_center`
    set symop = $symops[$symop]
    
    # progress meter
    echo "$right_chain_ID $wrong_chain_ID  $reindexing  $Xf $Yf $Zf $symop" |\
     awk '{printf "%s vs %s by %-20s @ %3s %3s %3s %-20s", $1, $2, $3, $4,$5,$6, $7}'
    
    # move the reindexed wrong PDB to the new position
    pdbset xyzin ${tempfile}wrong_reindexed.pdb \
          xyzout ${tempfile}symmed.pdb << EOF >! ${tempfile}.log
    CELL $CELL
    SYMGEN $symop
    SHIFT FRAC $X $Y $Z
    COM
EOF
    set new_COM = `awk '$1=="Center"{print $4, $5, $6}' ${tempfile}.log`
    # get fractional version too
    awk '/^CRYST/ || /^SCALE/' ${tempfile}wrong_reindexed.pdb >! ${tempfile}COM.pdb
    cat ${tempfile}.log |\
    awk '$1=="Center"{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\n", $4, $5, $6;exit}' |\
    cat >> ${tempfile}COM.pdb
    echo "END" >> ${tempfile}COM.pdb
    coordconv XYZIN ${tempfile}COM.pdb XYZOUT ${tempfile}.xyz << EOF > /dev/null
    CELL $CELL
    INPUT PDB
    OUTPUT FRAC
EOF
    set new_COM_frac = `awk '{print "  ", $2, $3, $4}' ${tempfile}.xyz`
    
    # now find the nearest cell translation to bring these COMs together
    set cell_shift = `echo "$right_COM_frac $new_COM_frac" | awk '{printf "%.0f %.0f %.0f", $1-$4+0, $2-$5+0, $3-$6+0}'`
    
    # move the "symmed" "wrong" chain into the right cell
    pdbset xyzin ${tempfile}symmed.pdb \
          xyzout ${tempfile}moved.pdb << EOF >! ${tempfile}.log
    CELL $CELL
    SHIFT FRAC $cell_shift
EOF
 
    # calclulate RMS distance between chains
    set score = `awk -f ${tempfile}rmsd.awk $right_chain ${tempfile}moved.pdb | awk '/RMSD\(all/{print -$2}'`

    if($?CORRELATE) then
        # make sure all atoms and residues are readabel
        cat ${tempfile}moved.pdb |\
        awk '/^ATOM|^HETAT/{$0="ATOM      1  CA  ALA     1    "substr($0,31,25)" 1.00  5.00"} {print}' |\
            awk '/^ATOM|^HETAT/{++n;}\
         n==1{++n;\
          print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\
          print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\
             {print}' |\
        cat >! ${tempfile}mapme.pdb
        # create a map of these atoms
        sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}moved.map << EOF-sfall > /dev/null
        MODE ATMMAP
        CELL $CELL
        SYMM $CCSG
        GRID $GRID
        BADD $BADD
EOF-sfall
        rm -f ${tempfile}mapme.pdb >& /dev/null
        
        # compare to the "right" map
        echo "correlate section" |\
        overlapmap mapin1 ${tempfile}right.map \
               mapin2 ${tempfile}moved.map mapout /dev/null |\
        awk '/Total correlation/{print $NF}' >! ${tempfile}correlation
        set score = `awk '{print $1}' ${tempfile}correlation`
        
        # display CC
        echo "$score" | awk '$1!~/[^0-9]/{print "  -";exit} {printf "%5.2f\n", $1}'
    else
        # display rmsd
        echo "$score" | awk '$1~/[0-9]/{printf "%5.1f\n", -$1; exit} {print "no identities"}'
    endif
    rm -f ${tempfile}symmed.pdb ${tempfile}moved.pdb
    if("$score" == "") set SKIP_CHAINPAIR = ( $right_chain $wrong_chain )
    if($?SKIP_CHAINPAIR && $skip_noid) break
    if("$score" == "") continue

    # see if fit is "good enough"
    set goodenough = `echo $score | awk '{print ($1*$1<1)}'`
    if($?CORRELATE) set goodenough = `echo $score | awk '{print ($1>0.8)}'`
    if(($goodenough == 1)&&($?SPEEDUP)) then
        echo "thats good enough..."
        set good_origin = "$X $Y $Z"
        set good_right = "$right_chain"
        set good_wrong = "$wrong_chain"
    endif
    
    # make a running log
    echo "$right_chain $wrong_chain $reindexing $Xf $Yf $Zf $X $Y $Z $cell_shift $symop $score " >> ${tempfile}scores

end
end
end
end

echo ""

# select out the "best" matches (on the same origin)
sort -nr -k14 ${tempfile}scores |\
awk '{ros=$3" "$4" "$5" "$6;chain_ros=$2" "ros} ! seen[chain_ros]{++seen[chain_ros];\
        score[ros]+=$NF;++count[ros]}\
     END{for(ros in score) print score[ros], count[ros], ros}' |\
sort -nr >! ${tempfile}best_oriscores 
# format: score n X,Y,Z dX dY dZ
set best_reindexing_origin = `head -1 ${tempfile}best_oriscores | awk '{print $3, $4, $5, $6}'`

if($?CORRELATE && "$CCSG" != "1" && ! $?NO_ORIGINS) then
    # need to go back and do this again for symops
    set reindexing_origins = `echo $best_reindexing_origin | awk '{print $1"_"$2","$3","$4}'`
    set CCSG = 1
    set good_right = ""
    set good_wrong = ""
    echo "now checking symops"
    goto again
endif

# extract the best operators consistent with this origin
sort -nr -k14 ${tempfile}scores |\
awk -v best_reindexing_origin="$best_reindexing_origin" '$3" "$4" "$5" "$6==best_reindexing_origin{print}' |\
awk '! seen[$2]{print; seen[$2]=1}' |\
sort >! ${tempfile}best_scores 
#format r_chain w_chain X,Y,Z 1/2 1/2 1/2 0.000 0.000 0.000  1 0 1  X,Y,Z  rmsd


echo -n "" >! ${tempfile}out.pdb
foreach line ( `awk '{print NR}' ${tempfile}best_scores` )
    # retrieve chain info
    set right_chain = `awk -v line=$line 'NR==line{print $1}' ${tempfile}best_scores`
    set wrong_chain = `awk -v line=$line 'NR==line{print $2}' ${tempfile}best_scores`

    set right_chain_ID = `awk '/^REMARK chain/{print $3}' $right_chain`
    set wrong_chain_ID = `awk '/^REMARK chain/{print $3}' $wrong_chain`
    if("$right_chain_ID" == "") set right_chain_ID = "_"
    if("$wrong_chain_ID" == "") set wrong_chain_ID = "_"
    

    # retrieve transofmation info
    set reindexing = `awk -v chain=$wrong_chain '$2==chain{print $3}' ${tempfile}best_scores`
    set X     = `awk -v chain=$wrong_chain '$2==chain{printf "%.2f", $7}' ${tempfile}best_scores`
    set Y     = `awk -v chain=$wrong_chain '$2==chain{printf "%.2f", $8}' ${tempfile}best_scores`
    set Z     = `awk -v chain=$wrong_chain '$2==chain{printf "%.2f", $9}' ${tempfile}best_scores`
    set x     = `awk -v chain=$wrong_chain '$2==chain{print $7+$10}' ${tempfile}best_scores`
    set y     = `awk -v chain=$wrong_chain '$2==chain{print $8+$11}' ${tempfile}best_scores`
    set z     = `awk -v chain=$wrong_chain '$2==chain{print $9+$12}' ${tempfile}best_scores`
    set symop = `awk -v chain=$wrong_chain '$2==chain{print $13}' ${tempfile}best_scores`

    set CELL = `awk -v key="$reindexing" '$1==key{print $3,$4,$5,$6,$7,$8;exit}' ${tempfile}xyz_hkl_cell.log`
    set reindexing_hkl = `awk -v key="$reindexing" '$1==key{print $2;exit}' ${tempfile}xyz_hkl_cell.log`

    # retrieve score
    set score  = `awk -v chain=$wrong_chain '$2==chain{print -$14}' ${tempfile}best_scores`
    if($?CORRELATE) set score = `echo $score | awk '{print -$1}'`
    
    if("$score" == "") then
        echo "unable to match $wrong_pdb chain $wrong_chain_ID to any chain in $right_pdb"
        continue
    endif
        
    # apply it
    echo "$reindexing $x $y $z $symop" |\
    awk '{printf "applying %-15s %5.2f %5.2f %5.2f %-15s\n", $1, $2, $3, $4, $5}'
    echo "to $wrong_pdb chain $wrong_chain_ID --> $outfile chain $right_chain_ID  ($SCORE = $score)"

    if("$right_chain_ID" == "_") set right_chain_ID = " "
    # apply this re-indexing to the "wrong" pdb
    pdbset xyzin ${wrong_chain} xyzout ${tempfile}lattcell.pdb << EOF > /dev/null
    CELL $lattCELL
EOF
    # apply the re-indexing operation, which is done in lattice cell
    pdbset xyzin ${tempfile}lattcell.pdb xyzout ${tempfile}symmed.pdb << EOF > /dev/null
    symgen $reindexing
EOF
    # and then put back the relevant cell and SG in the header
    pdbset xyzin ${tempfile}symmed.pdb xyzout ${tempfile}reindexed.pdb << EOF > /dev/null
    CELL $CELL
    space $SG
EOF
    # now apply the symop and origin shift
    pdbset xyzin ${tempfile}reindexed.pdb xyzout ${tempfile}moved.pdb << EOF >> /dev/null
    CELL $CELL
    CHAIN "$right_chain_ID"
    SYMGEN $symop
    SHIFT FRAC $x $y $z
EOF
    # update the output file
    egrep "^ATOM|^HETAT" ${tempfile}moved.pdb >> ${tempfile}out.pdb
    
    # warn about different origins
    if(! $?last_reindexing_origin) set last_reindexing_origin
    if(("$last_reindexing_origin" != "$reindexing $X $Y $Z")&&("$last_reindexing_origin" != "")) then
        echo "WARNING: $wrong_pdb chain $wrong_chain_ID is on a different origin! "
    endif
    set last_reindexing_origin = "$reindexing $X $Y $Z"
end
echo "END" >> ${tempfile}out.pdb

pdbset xyzin ${tempfile}out.pdb xyzout $outfile << EOF > /dev/null
CELL $CELL
SPACE $SG
EOF

if($?BYFILE) then
    # no point in messing around, just apply the shift to the original file
    pdbset xyzin $wrong_pdb xyzout ${tempfile}lattcell.pdb << EOF > /dev/null
    CELL $lattCELL
EOF
    # apply the re-indexing operation, which is done in lattice cell
    pdbset xyzin ${tempfile}lattcell.pdb xyzout ${tempfile}symmed.pdb << EOF > /dev/null
    symgen $reindexing
EOF
    # and then put back the relevant cell and SG in the header
    pdbset xyzin ${tempfile}symmed.pdb xyzout ${tempfile}out.pdb << EOF > /dev/null
    CELL $CELL
    space $SG
EOF
    # and now apply the symop and origin shift
    pdbset xyzin ${tempfile}out.pdb xyzout $outfile << EOF > /dev/null
CELL $CELL
SPACE $SG
SYMGEN $symop
SHIFT FRAC $x $y $z
EOF

# now recover the operation in other conventions
lsqkab XYZIN1 $outfile XYZIN2 $wrong_pdb  << EOF >! ${tempfile}lsq.log
FIT ATOM 1 to 99999
MATCH    1 to 99999
END
EOF
set lsq_polar = `awk '/OMEGA PHI CHI/ && $1~/^SPHERICAL/{print $(NF-2), $(NF-1), $NF} /ALPHA=BETA=GAMMA=CHI=0.0/{print 0,0,0}' ${tempfile}lsq.log`
set lsq_trans = `awk '/TRANSLATION VECTOR IN AS/{print $(NF-2), $(NF-1), $NF}' ${tempfile}lsq.log`
cat << EOF
general transformation for pdbset, maprot, ncsmask, dm, etc.:
ROTATE POLAR $lsq_polar
TRANSLATE $lsq_trans
EOF
endif

echo "reindexing operator for data: $reindexing_hkl"

if($?SPEEDUP) goto exit

# report final, overall score
if($?CORRELATE) then
    # calculate the overall map correlation
    set logfile = /dev/null

    # create a simplified file for SFALL, cast into the final cell
    egrep "^CRYST1" $outfile | head -n 1 >! ${tempfile}mapme.pdb
    cat $right_pdb |\
    awk '/^END|^REM/{print}\
            /^ATOM|^HETAT/{++n;printf("ATOM      1  CA  ALA %5d    %25s 1.00 80.00\n",n,substr($0,31,25))}' |\
    awk '/^ATOM|^HETAT/{++n;}\
     0 && n==1{++n;\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\
         {print}' |\
    cat >> ${tempfile}mapme.pdb
    sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}right.map << EOF-sfall > $logfile
    MODE ATMMAP
#    CELL $right_cell
    CELL $CELL
    SYMM $SG
EOF-sfall
    # also record the "right" atom codes
    cat $right_pdb |\
    awk '/^ATOM|^HETAT/{++n;printf "code %04d %s\n",n,$0}' |\
    cat >! ${tempfile}atomcodes.txt
    # recover grid spacing 
    set GRID = `echo "GO" | mapdump MAPIN ${tempfile}right.map | awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}'`

    # calculate the "label" map so we can assign atom-specific correlations
    sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}label.map << EOF-sfall >> $logfile
    MODE ATMMAP RESMOD
    CELL $CELL
    SYMM $SG
    GRID $GRID
EOF-sfall

    # calculate map for "moved" model
    cat $outfile |\
    awk '/^CRYST|^END|^REM/{print}\
            /^ATOM|^HETATM/{++n;printf("ATOM      1  CA  ALA %5d    %25s 1.00 80.00\n",n,substr($0,31,25))}' |\
    awk '/^ATOM|^HETAT/{++n;}\
     0 && n==1{++n;\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\
         {print}' |\
    cat >! ${tempfile}mapme.pdb
    # also record the "wrong" atom codes
    cat $outfile |\
    awk '/^ATOM|^HETAT/{++n;printf "code %04d %s\n",n,$0}' |\
    cat >! ${tempfile}wrong_atomcodes.txt
    sfall xyzin ${tempfile}mapme.pdb mapout ${tempfile}moved.map << EOF-sfall > $logfile
    MODE ATMMAP
    CELL $CELL
    SYMM $SG
    GRID $GRID
EOF-sfall
# why do they make me do these things... sfall map can sometimes be too small
#dd if=/dev/zero bs=1024 count=10000 >>& ${tempfile}moved.map
#    sfall xyzin ${tempfile}mapme.pdb HKLOUT ${tempfile}moved.mtz << EOF-sfall > $logfile
#    MODE SFCALC XYZIN
#    CELL $CELL
#    SYMM $CCSG
#EOF-sfall
#    fft hklin ${tempfile}moved.mtz mapout ${tempfile}moved.map << EOF-fft > $logfile
#    LABIN F1=FC PHI=PHIC
#    GRID $GRID
#EOF-fft


    rm -f ${tempfile}mapme.pdb >& /dev/null
    # compare to the "right" map
    echo "correlate residue" |\
    overlapmap mapin1 ${tempfile}right.map \
           mapin2 ${tempfile}moved.map \
           mapin3 ${tempfile}label.map mapout /dev/null |\
    awk '/Main-corr-coef/,/Total/{print}' | tee ${tempfile}correlation |\
    awk '$1+0>0{print}' |\
    awk '$2+0>0{print $1,$2;next} $3+0>0{print $1,$3}' |\
    awk '{print "ATOMCC",$0}' >! ${tempfile}atom_corr.txt
    # format: ATOMCC rightatomnum CC

    # encode the atomic CCs as an occupancy in a new PDB file
    cat ${tempfile}atomcodes.txt ${tempfile}atom_corr.txt $right_pdb |\
    awk '/^code /{name[$2+0]=substr($0,11,26);next} /^ATOMCC /{occ[$2]=$3;next} ! /^ATOM/ {print}\
      /^ATOM|^HETATM/{++n;\
    if(name[n]=="" || occ[n]+0<0.01){name[n]=substr($0,1,26);occ[x]=substr($0,55,6)+0};\
    printf "%26s%28s%6.2f%s\n", name[n],substr($0,27,28),occ[n],substr($0,61)}' |\
    cat >! newlabel_$outfile

    set score = `awk '/Total/{print $NF;print $(NF-1)}' ${tempfile}correlation | sort -nr | head -1`

    rm -f ${tempfile}right.map >& /dev/null
    rm -f ${tempfile}moved.map >& /dev/null
    rm -f ${tempfile}label.map >& /dev/null
    rm -f ${tempfile}correlation >& /dev/null    
    rm -f ${tempfile}atom_corr.txt >& /dev/null    
    rm -f ${tempfile}atomcodes.txt >& /dev/null    
    rm -f ${tempfile}wrong_atomcodes.txt >& /dev/null
else
    # the overall RMSD is the score
    set score = `awk -f ${tempfile}rmsd.awk $right_pdb $outfile | awk '/RMSD\(all/{print $2} /no atom pairs/{print "no identities"}'`
endif
echo "Overall ${SCORE}: $score"



#clean up
exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif
if($?debug) exit
rm -f ${tempfile}rmsd.awk >& /dev/null
rm -f ${tempfile}out.pdb >& /dev/null
rm -f ${tempfile}moved.pdb >& /dev/null
rm -f ${tempfile}best_scores >& /dev/null
rm -f ${tempfile}scores >& /dev/null
rm -f ${tempfile}.log >& /dev/null
rm -f ${tempfile}.xyz >& /dev/null
rm -f ${tempfile}COM.pdb >& /dev/null
rm -f ${tempfile}symop_center >& /dev/null
rm -f ${tempfile}_wrong???.pdb >& /dev/null
rm -f ${tempfile}_right???.pdb >& /dev/null
rm -f ${tempfile}* >& /dev/null

exit

# Notes on hand flipping

cat rh.pdb |\
awk -v SG=$SG '/^CRYST1/{a=$2;c=$4;\
     if(SG=="F4132"){dorx=dory=dorz=a*0.25};\
     if(SG=="I41"){dorx=a*0.5};\
     if(SG=="I4122"){dorx=a*0.5;dorz=c*0.25};\
  } ! /^ATOM|^HETAT/{print} /^ATOM|^HETAT/{;\
  x=substr($0,31,8)+0;y=substr($0,39,8)+0;z=substr($0,47,8)+0;\
  printf("%s%8.3f%8.3f%8.3f%s\n", substr($0,1,30),dorx-x,dory-y,dorz-z,substr($0,55))}' |\
cat >! lh.pdb




################################################################################

Setup:
set good_origin = ""
set good_right = ""
set good_wrong = ""

foreach arg ( $* )

    if(("$arg" =~ *.pdb)||("$arg" =~ *.brk)) then
    # warn about probable mispellings
    if(! -e "$arg") then
        echo "WARNING: $arg does not exist"
        continue
    endif
    # make sure its really a pdb file
    egrep -l "^ATOM |HETATM " "$arg" >& /dev/null
    if($status) then
        echo "WARNING: $arg contains no atoms! "
        continue
    endif
    
    if(-e "$right_pdb") then
        set wrong_pdb = "$arg"
    else
        set right_pdb = "$arg"
    endif
    continue
    endif

    # space group
    if("$arg" =~ [PpCcIiFfRrHh][1-6]*) then
    set temp = `echo $arg | awk '{print toupper($1)}'`
    if($?CLIBD) then
        set temp = `awk -v SG=$temp '$4 == SG {print $4}' $CLIBD/symop.lib | head -1`
    endif
    if("$temp" != "") then
        # add this SG to the space group list
        set SG = "$temp"
    endif
    endif

    # flags and options
    if("$arg" == "correlate") then
        # correlation instead of RMSD
        set CORRELATE
    endif
    if("$arg" == "otherhand") then
        set OTHERHAND
    endif
    if("$arg" == "altindex") then
        set ALTINDEX
    endif
    if("$arg" == "nochains" || "$arg" == "byfile") then
        # whole file is one chain
        set BYFILE
    endif
    if("$arg" == "noorigins" || "$arg" == "noorigin") then
        # just symmetry search
        set NO_ORIGINS
    endif
    if("$arg" =~ nodecon*) then
        # dangerous for polar space groups
        set NO_DECONV
    endif
    if("$arg" == "fast") then
        # skip several steps
        set SPEEDUP
    endif
    if("$arg" =~ reso=*) then
        # fuzziness of correlation
        set test = `echo $arg | awk -F "=" '$2+0>0{print $2+0}'`
        if("$test" != "") then
            set reso = $test
        else
            echo "WARNING: cannot understand $arg"
        endif
    endif
    if("$arg" =~ min_opeak=*) then
        # limit the number of possible origins
        set test = `echo $arg | awk -F "=" '$2+0>0{print $2+0}'`
        if("$test" != "") then
            set min_opeak = $test
        else
            echo "WARNING: cannot understand $arg"
        endif
    endif
    if("$arg" =~ max_opeaks=*) then
        # limit the number of possible origins
        set test = `echo $arg | awk -F "=" '$2+0==int($2+0) && $2+0>0{print int($2+0)}'`
        if("$test" != "") then
            set max_opeaks = $test
        else
            echo "WARNING: cannot understand $arg"
        endif
    endif
end

if((! -e "$right_pdb")||(! -e "$wrong_pdb")) then
    goto Help
endif


# deploy rmsd awk program
cat << EOF-script >! ${tempfile}rmsd.awk
#! $awk -f
#
#   Calculate RMSD of atoms with the same name in two PDB files
#
#    The PDB feild:
#           |<--- here -->|
#ATOM      1  N   ALA A 327      40.574  34.523  43.012  1.00 34.04
#
#    is used to determine if two atoms are a "pair"
#
BEGIN {
if(! atom) atom = "CA"
maxXYZ = maxdB = 0
max_atom_XYZ = max_atom_dB = 0
maxXYZ_ID = maxdB_ID = "-"
max_atom_XYZ_ID = max_atom_dB_ID = "-"
}

/^ATOM/{
    # read in values (watching for duplications)
    ID = substr(\$0,12,15)
    ++count[ID]

    if(count[ID] == 1)
    {
    # not initialized yet
    X[ID] = substr(\$0, 31, 8)+0
    Y[ID] = substr(\$0, 39, 8)+0
    Z[ID] = substr(\$0, 47, 8)+0
    B[ID] = substr(\$0, 61, 6)+0
    }
    
    if(count[ID] == 2)
    {
    ++pairs
    
    # seen this before, subtract values
    dX     = X[ID] - substr(\$0, 31, 8)
    dY     = Y[ID] - substr(\$0, 39, 8)
    dZ     = Z[ID] - substr(\$0, 47, 8)
    dB[ID] = B[ID] - substr(\$0, 61, 6)
    
    # get drift (and add up squares of drifts)
    sqrD   = dX*dX + dY*dY + dZ*dZ
    dXYZ[ID] = sqrt(sqrD)

    # remember maximum shifts
    if(dXYZ[ID] > maxXYZ) {maxXYZ  = dXYZ[ID]; maxXYZ_ID = ID }
    if(dB[ID]*dB[ID] > maxdB*maxdB) {maxdB = dB[ID]; maxdB_ID = ID }

    # maintain mean-square sums    
    sumXYZ += sqrD
    sumB   += dB[ID]*dB[ID]

    # separate stats for special atom type
    if(ID ~ atom)
    {
        ++atom_pairs

        # maintain separate mean-square sums
        sum_atom_XYZ += sqrD
        sum_atom_B   += dB[ID]*dB[ID]
        
        # remember maximum drifts too
        if(dXYZ[ID] > max_atom_XYZ) {max_atom_XYZ  = dXYZ[ID]; max_atom_XYZ_ID = ID }
        if(dB[ID]*dB[ID] > max_atom_dB*max_atom_dB) {max_atom_dB = dB[ID]; max_atom_dB_ID = ID }        
    }
    # debug output
    if(debug)
    {
        printf("%s moved %8.4f (XYZ) %6.2f (B)\\n", ID, dXYZ[ID], dB[ID])
    }    
    }

    if(count[ID] > 2)
    {
    print "WARNING: " ID " appeared more than twice! "
    }
}


END{
    
    if((pairs+0 == 0)&&(! xlog)) 
    {
    print "no atom pairs found"
    exit
    }
    rmsXYZ = sqrt(sumXYZ/pairs)
    rmsB = sqrt(sumB/pairs)
    if(atom_pairs+0 != 0)
    {
    rms_atom_XYZ = sqrt(sum_atom_XYZ/atom_pairs)
    rms_atom_B = sqrt(sum_atom_B/atom_pairs)
    }
    

    if(! xlog) 
    {
    print pairs " atom pairs found"
    print "RMSD("atom" )= " rms_atom_XYZ " ("atom_pairs, atom " pairs)"
    print "RMSD(all)= " rmsXYZ " ("pairs" atom pairs)"
    print "RMSD(Bfac)= " rmsB

    print "MAXD(all)= " maxXYZ "\\tfor " maxXYZ_ID
    print "MAXD(Bfac)= " maxdB "\\tfor " maxdB_ID
    
    # final check for orphan atoms 
    for(id in count)
    {
        if(count[id]<2) print "WARNING: " id " only found once"
    }
    }
    else
    {
    printf "%10.8f %10.8f %10.5f %10.8f %8.2f \\n", rms_atom_XYZ, rmsXYZ, rmsB, maxXYZ, maxdB
    }
}
EOF-script
chmod a+x ${tempfile}rmsd.awk


goto Return_from_Setup


# TABLE OF ALLOWED ORIGIN SHIFTS
These origin shifts were determined emprirically using 100 randomly-placed atoms
that were shifted around with pdbset and checked with SFALL for identical 
amplitudes to the 0 0 0 origin.  They should be correct for the CCP4 convention
of symmetry. (I.E. R3 and R32 have a=b=c, but H3 and H32 have gamma=120)
Note that 0.25 0.25 0.25 is not an allowed origin shift for F4132, despite
that it may be convenient to check this when inverting the hand. 

P1 
   x   y   z

P2 P21 C2 
   0   y   0
   0   y 1/2
 1/2   y   0
 1/2   y 1/2

P222 P2221 P21212 P212121 C2221 C222 I222 I212121 F432 F4132 
   0   0   0
   0   0 1/2
   0 1/2   0
   0 1/2 1/2
 1/2   0   0
 1/2   0 1/2
 1/2 1/2   0
 1/2 1/2 1/2

F222 F23  
   0   0   0
   0   0 1/2
   0 1/2   0
   0 1/2 1/2
 1/2   0   0
 1/2   0 1/2
 1/2 1/2   0
 1/2 1/2 1/2
 1/4 1/4 1/4
 1/4 1/4 3/4
 1/4 3/4 1/4
 1/4 3/4 3/4
 3/4 1/4 1/4
 3/4 1/4 3/4
 3/4 3/4 1/4
 3/4 3/4 3/4

P4 P41 P42 P43 I4 I41 
   0   0   z
 1/2 1/2   z

P422 P4212 P4122 P41212 P4222 P42212 P4322 P43212 I422 I4122 
   0   0   0
   0   0 1/2
 1/2 1/2   0
 1/2 1/2 1/2

P3 P31 P32 
   0   0   z
 1/3 2/3   z
 2/3 1/3   z

R3 
 x=y=z

H3 
   0   0   z
 1/3 2/3   z
 2/3 1/3   z

P312 P3112 P3212 
   0   0   0
   0   0 1/2
 1/3 2/3   0
 2/3 1/3   0
 1/3 2/3 1/2
 2/3 1/3 1/2

P321 P3121 P3221 P622 P6122 P6522 P6222 P6422 P6322 
   0   0   0
   0   0 1/2

R32 
   0   0   0
 1/2 1/2 1/2

H32 
   0   0   0
   0   0 1/2
 1/3 2/3 2/3
 2/3 1/3 1/3
 1/3 2/3 1/6
 2/3 1/3 5/6

P6 P61 P65 P62 P64 P63 
   0   0   z

P23 P213 P432 P4232 I432 P4332 P4132 I4132 
   0   0   0
 1/2 1/2 1/2

END OF THE TABLE
