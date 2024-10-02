#! /bin/tcsh -f
#
##############################################################################
#    Find a smooth, periodic function to describe the shift 
#    between two non-isomorphous forms of the same structure
#
#    Function is based on HKL index equations, but not directly related 
#    to low-resolution structure factors
#
#    -James Holton  5-13-22
#    -comments by David Case, Aug. 2017
#
#    Usage : bendfinder.com bendme.pdb reference.pdb [bendme.map] 
#                          [order=n] [frac=x] [reject=m] [dimensions=list]
#                          [fitparams_3.gnuplot]
#
#    where :
#     order   number of orders of 3D polysine functions to use in the fit. 
#             Default : $maxorder, set in script as 8
#     frac    how far along the bend path to place the output.  
#             Default : $frac, set in script as 1
#     geotest use refmac5 to evaluate geometry at each step
#             Default : $geotest
#     reject  reject outliers more than m MADs from the median
#     dimensions  list of parameters to include
#     fitparams_3.gnuplot  previously fitted bend function to apply to a provided
#             pdb or map file.  Skips the long fitting step requiring gnuplot5
#
#    requires :
#        gnuplot version 5 to be in the PATH as "gnuplot" or "gnuplot5"
#        ccp4 programs to be in the PATH : coordconv, pdbset, mapmask, mapdump
#
#        if bendme.map is provided :
#           mapman from the RAVE suite from the Uppsala Software Factory
##############################################################################

#========================================================================
#    Setting defaults :

# requires gnuplot version 5
set path = ( $path `dirname $0` . /programs/gnuplot5/bin/ )

# location of mapman?
set mapman = mapman

# how many orders of hkl-like functions to use
set maxorder = 8
set startorder = 0

# how far across the transition to go?
set frac = 1

# reject outliers more than m MADs from the median
set reject = 0

# speed things up with higher gnuplot FIT_LIMIT ?
set FIT_LIMIT = ""
set fitscale  = 1000

# maybe decide to ignore some params
set dimensions = ( x y z o B )
set dimensions = ( x y z )

# do geometry test (can be slow)
set geotest = true

# movable PDB
set pdb1 = ""
# reference PDB
set pdb2 = ""
# 2fofc map for pdb1 to be re-sampled back into frame of pdb2
set mapfile = ""
# user-provided fitparams
set user_fitparams = ""

# option to make four "delta maps" for each coordinate shift in space
# that is, delta-x delta-y delta-z, and the rms of all three: delta-r
#set MAKE_DELTAMAPS

set starttime = `date +%s`

#========================================================================
#    help message if run with no arguments

if($#argv == 0) then
    cat << EOF 
usage: $0 bendme.pdb reference.pdb [bendme.map] [order=n] \
      [frac=x] [reject=6] [dimensions=xyz] \
      [fitparams_x.gnuplot] [refit] [deltamaps]

where:
bendme.pdb     defines the frame to be distorted to match the reference
reference.pdb  defines the reference frame and cell

bendme.map     is a map that will be spline interpolated into the reference frame

order          is the number of orders of 3D polysine functions to use in the fit. Default: $maxorder
frac           is how far along the bend path to place the output.  Default: $frac
reject         reject outlier shifts > 6x the median absolute deviation. Default: no rejection
dimensions     aspects of PDB files to fit, default: xyz (not occ or B)
fitparams_x.gnuplot  - use a previously-fitted function, skip the fit step
refit          do the fit step even if a fitparams.gnuplot is provided
deltamaps      generate electron density maps of the shift magnitudes in x,y,z and r over whole cell

both PDB files must always be provided because they contain the two unit cells.
EOF
    exit 9
endif

if(! $?CCP4_SCR) then
    set BAD = "CCP4 is not set up."
    goto exit
endif

# pick temp filename and decide on debug log
set logfile = /dev/null
mkdir -p ${CCP4_SCR} >&! /dev/null
set tempfile = ${CCP4_SCR}/bendfinder$$temp
if(-w /dev/shm/ ) then
    set tempfile = /dev/shm/bendfinder$$temp
endif

#========================================================================
# find right version of gnuplot
set gnuplot_version  = `echo show version | gnuplot  |& awk '$1=="Version"{print $2}'`
set gnuplot5_version = `echo show version | gnuplot5 |& awk '$1=="Version"{print $2}'`
set test = `echo $gnuplot_version | awk '{print ($1>=5)}'`
if("$test" == "1") then
    # gnuplot command runs, and is actually version 5
    alias gnuplot5 gnuplot
    set gnuplot5 = gnuplot
else
    # gnuplot is not version 5, or may not exist
    set test = `echo $gnuplot5_version | awk '{print ($1>=5)}'`
    if("$test" != "1") then
        set BAD = "cannot find gnuplot version 5"
        goto exit
    endif
    # gnuplot5 command runs gnuplot version 5
    if("$gnuplot_version" == "") then
        # gnuplot command doesnt work, but gnuplot5 is installed
        alias gnuplot gnuplot5
    endif
    set gnuplot5 = gnuplot5
endif


#========================================================================
# find or generate other needed items

set test = `echo 0 | floatgen | od -c | awk 'NR==1{print ( / \\0  \\0  \\0  \\0/ )}'`
if("$test" != "1") then
    echo "floatgen doesnt work.  Trying to create it."
    if(! $?DEPLOYED) goto deploy_scripts
endif

set test = `echo -n "" | map2pdb.com | egrep "one atom per map voxel" | wc -l`
if("$test" != "1") then
    echo "map2pdb.com script not found.  Trying to create it."
    if(! $?DEPLOYED) goto deploy_scripts
endif

set test = `echo -n "" | rmsd | egrep "no atom pairs found" | wc -l`
if("$test" != "1") then
    echo "rmsd script not found.  Trying to create it."
    if(! $?DEPLOYED) goto deploy_scripts
endif


return_from_deploy:

set test = `echo 0 | floatgen | od -c | awk 'NR==1{print ( / \\0  \\0  \\0  \\0/ )}'`
if("$test" != "1") then
    set BAD = "floatgen does not work.  Please compile http://bl831.als.lbl.gov/~jamesh/bin_stuff/floatgen.c"
    goto exit
endif

set test = `echo -n "" | map2pdb.com | egrep "one atom per map voxel" | wc -l`
if("$test" != "1") then
    set BAD = "map2pdb.com script not found.  Please download from: http://bl831.als.lbl.gov/~jamesh/scripts/map2pdb.com"
    goto exit
endif

set test = `echo -n "" | rmsd | egrep "no atom pairs found" | wc -l`
if("$test" != "1") then
    set BAD = "rmsd script not found.  Please download from: http://bl831.als.lbl.gov/~jamesh/scripts/rmsd"
    goto exit
endif

#========================================================================
#    Reading command-line arguments:

foreach arg ( $* )
    if("$arg" =~ *.pdb) then
        if(-e "$pdb1" && -e "$pdb2") then
            set outpdb = $arg
            continue
        endif
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        if("$pdb1" == "") then
            set pdb1 = "$arg"
        else
            set pdb2 = "$arg"
        endif
        continue
    endif
    if("$arg" =~ *.map) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set mapfile = "$arg"
    endif
    if("$arg" =~ *.gnuplot) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set user_fitparams = "$arg"
    endif
    if("$arg" == debug) then
        set debug = 1
    endif
    if("$arg" =~ dimensions=*) then
        set dimensions = `echo $arg | awk -F "=" '{print $2}' | awk -F "," '{for(i=1;i<=NF;++i)print $i}'`
    endif
    if("$arg" =~ order=*) set maxorder = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ frac=*) set frac = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ reject=*) set reject = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ startorder=*) set startorder = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ geotest=*) set geotest = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ fit_limit=*) set FIT_LIMIT = `echo $arg | awk -F "=" '{print "FIT_LIMIT=",$2}'`
    if("$arg" =~ fitscale=*) set fitscale = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" == refit) set REFIT
    if("$arg" == deltamaps) set MAKE_DELTAMAPS
    if("$arg" == nodeltamaps) unset MAKE_DELTAMAPS
end

if($?debug) then
    set tempfile = tempfile
    set logfile = /dev/stdout
endif


#========================================================================
#    check if we need mapman

if(-e "$mapfile") then
    set test = `echo quit | $mapman |& egrep "Toodle pip" | wc -l`
    if("$test" != "1") then
        # try another name?
        set mapman = lx_mapman
        set test = `echo quit | $mapman |& egrep "Toodle pip" | wc -l`
    endif
    if("$test" != "1") then
        set BAD = "unable to run mapman"
        goto exit
    endif
endif




#========================================================================
#    first, check if we have two cells
if(! -e "$pdb2") then
    set BAD = "two PDB files are required to provide both unit cells"
    goto exit
endif

# check if PDB files are ridiculously different
set test = `rmsd $pdb1 $pdb2 | awk '/RMSD.CA/ && $3+0<10' | wc -l`
if("$test" != "1") then
    set BAD = "pdb files too different.  Maybe run origins.com first."
    goto exit
endif

# come up with shorthand names
set id1 = `basename $pdb1 .pdb`
set id2 = `basename $pdb2 .pdb`
if("$id1" == "$id2" || "$id1" == "" || "$id2" == "") then
    set id1 = "${id1}1"
    set id2 = "${id2}2"
endif

#========================================================================
#    extract symmetry from input pdb files and convert to fractional coordinates:

set n = 0
set ids = ( $id1 $id2 ) 
foreach pdb ( $pdb1 $pdb2 )
    @ n = ( $n + 1 )

    if(! -e "$pdb") continue

    set id = $ids[$n]

    set CELL = `awk '/CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdb`
    set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdb | head -1`
    if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
    if("$pdbSG" == "P 21") set pdbSG = "P 1 21 1"
    if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
    set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
    if("$SG" == R3 && $CELL[6] == 120.00) set SG = H3
    if("$SG" == "") set SG = P1

    # duplicates dont matter, only ordinal number
    coordconv xyzin $pdb xyzout ${id}.frac << EOF > $logfile
    CELL $CELL 
    INPUT PDB
    OUTPUT FRAC
    END
EOF

    # eliminate all duplicate occurrences of coordconv-equivalent atom IDs
    cat $pdb |\
    awk '/^ANIS/{next}\
      ! /^ATOM|^HETAT/{print;next}\
      {id=substr($0,12,5) substr($0,18,9)}\
     ! seen[id]{print} {++seen[id]}' |\
    cat >! ${tempfile}deduped.pdb

    echo "expanding $pdb from $SG to P1"
    echo "SYMGEN $SG" |\
     pdbset xyzin ${tempfile}deduped.pdb xyzout P1.pdb > $logfile

    coordconv xyzin P1.pdb xyzout ${tempfile}.frac << EOF > $logfile
    CELL $CELL 
    INPUT PDB
    OUTPUT FRAC
    END
EOF
    cat ${tempfile}.frac |\
    awk '{id=substr($0,61,15);print $0,++seen[id]}' |\
    cat >! ${id}_P1.frac

end


################################################################################
# catch user-provided fitparams, reconstruct functions, and skip ahead

if(-e "$user_fitparams") then
    echo "using existing $user_fitparams "
    cat $user_fitparams |\
    awk '/^a/{hkl=substr($1,5);x=substr($1,4,1);a[x hkl]=$3}\
       /^phi/{hkl=substr($1,7);x=substr($1,6,1);p[x hkl]=$3}\
       END{for(xhkl in a) print xhkl,a[xhkl]+0,p[xhkl]+0}' |\
    sort |\
    tee ${tempfile}apsi.txt |\
    awk '{print "a_d"$1" = "$2;\
        print "phi_d"$1" = "$3;}' >! ${tempfile}user_fitparams.gnuplot
    cp -p ${tempfile}user_fitparams.gnuplot fitparams.gnuplot

    cat $user_fitparams |\
    awk '/^a/{hkl=substr($1,length($1)-2);\
       print substr(hkl,1,1),substr(hkl,2,1),substr(hkl,3,1);}' |\
    sort -u | egrep -v '^0 0 0$' |\
    cat >! hkl.txt
    set order = `awk '{print $1;print $2;print $3}' hkl.txt | sort -gr | awk '{print;exit}'`

    # construct bending functions from these coefficients
    rm -f func.gnuplot
    foreach x ( x y z o B )
        set mult = $fitscale
        if( "$x" == "B" ) set mult = 1
        # define the function to fit in gnuplot5
        egrep "^$x" ${tempfile}apsi.txt |\
        awk -v mult=$mult 'NR==1{x=substr($1,1,1);\
                print "d"x"(x,y,z)=real(\\"}\
            {h=substr($1,2,1);k=substr($1,3,1);l=substr($1,4,1);\
                printf(" +a_d%s%d%d%d/%s",x,h,k,l,mult);\
                if($3+0==0){print "\\";next}\
                printf("*exp(2*pi*i*(%d*x+%d*y+%d*z+phi_d%s%d%d%d))\\\n",h,k,l,x,h,k,l);\
            } END{print ")"}' |\
        tee -a func.gnuplot >> $logfile
    end

    set reject = 0
    set maxorder = $order
    set startorder = $order
    if($?REFIT) then
        echo "re-fitting, starting from $user_fitparams"
        goto measure_shifts
    endif
    # no fitting jobs to run
    if(-e "$pdb1") then
        echo "applying $user_fitparams to $pdb1"
        goto bend_pdb
    endif
    if(-e "$mapfile") then
        echo "applying $user_fitparams to $mapfile"
        goto bend_map
    endif
    set BAD = "please provide a pdb or map to bend."
    goto exit
endif


measure_shifts:
#========================================================================
#    Compare the input frame to the reference frame:

# get the shifts in x,y,z,o,B for every atom in common
cat ${id1}_P1.frac ${id2}_P1.frac |\
awk '{ID=substr($0,61,20);\
     x=substr($0,6,10);y=substr($0,16,10);z=substr($0,26,10);\
    o=substr($0,46,5);b=substr($0,36,10);++seen[ID]}\
    seen[ID]==1{X[ID]=x;Y[ID]=y;Z[ID]=z;Occ[ID]=o;Bfac[ID]=b;next}\
    seen[ID]==2{\
       print X[ID],Y[ID],Z[ID],x-X[ID],y-Y[ID],z-Z[ID],o-Occ[ID],b-Bfac[ID],"|" ID;\
       next}\
    {print "# error",ID,"seen",seen[ID],"times"}' |\
grep -v HOH |\
cat >! fitme
# format xf yf zf  dxf dyf dzf do dB |ID

set bigfrac = `awk '/CA/{print sqrt($4*$4+$5*$5+$6*$6),"for",substr($0,index($0,"|")+1)}' fitme | sort -gr | awk '{print;exit}'`
echo "largest fractional CA shift: $bigfrac"
set test  = `echo $bigfrac 0.1 | awk '{print ($1>$NF)}'`
if( $test ) then
    set BAD = "largest fractional shift of $bigfrac cells is too big."
    goto exit
endif

# make sure this does not exist
rm -f ${tempfile}remaining.pdb >& /dev/null

#========================================================================
#    Big loop over orders:

echo -n "" >! fitparams.gnuplot
foreach order ( `seq $startorder 1 $maxorder` ) 

set deltaT = `date +%s | awk -v s=$starttime '{print $1-s}'`
echo ""
echo "=================  Starting order $order fit at $deltaT s"
echo ""

# now generate the basis functions
echo "generating basis functions out to order $order"
echo $order |\
awk '{m=$1;\
    for(h=0;h<=m;++h)for(k=0;k<=m;++k)for(l=0;l<=m;++l){\
      r=sqrt(h*h+k*k+l*l);\
      if(r==0 || r>m) continue;\
      print h,k,l;}}' |\
tee hkl.txt |\
  wc -l |\
  awk -v d="$dimensions" '{print $1+1," complex coefficients for each dimension:",d}'

foreach x ( x y z o B )

    # define new set to be filled in with last iteration
    echo $x |\
    cat - hkl.txt |\
    awk 'NR==1{x=$1;print "a_d"x"000 = new";next}\
       {h=$1;k=$2;l=$3;\
          printf("a_d%s%d%d%d = new\n",x,h,k,l);\
          printf("phi_d%s%d%d%d = new\n",x,h,k,l);\
       }' |\
    tee -a fitparams.gnuplot >> $logfile
end


# filter out duplicates, why would there be any?
cat fitparams.gnuplot |\
awk '! seen[$1]{print;++seen[$1]}' |\
cat >! ${tempfile}
mv ${tempfile} fitparams.gnuplot


if(-e "$user_fitparams") goto after_rejects

#========================================================================
#    do slow-ft for starting values
echo "doing slow-ft to get starting values..."
echo -n "" >! fitparams_init.gnuplot
foreach h_k_l ( 0_0_0 `awk '{print $1"_"$2"_"$3}' hkl.txt` )
    set hkl = `echo $h_k_l | awk -F "_" '{print $1,$2,$3}'`
    set h = $hkl[1]
    set k = $hkl[2]
    set l = $hkl[3]

    # skip if we already have a value
    set test = `egrep "a_dx${h}${k}${l} " fitparams.gnuplot | awk '{print ($NF+0 != 0)}'`
    if("$test" == "1") continue

    echo "$hkl  $fitscale" |\
    cat - fitme |\
    awk 'NR==1{h=$1;k=$2;l=$3;mult=$4;pi=atan2(1,1)*4;h000=(h==0 && k==0 && l==0);next;}\
      {x=$1;y=$2;z=$3;++n;\
      cosine=cos(2*pi*(h*x+k*y+l*z));\
      sine  =sin(2*pi*(h*x+k*y+l*z));\
      sumxc+=$4*cosine;\
      sumxs+=$4*sine;\
      sumyc+=$5*cosine;\
      sumys+=$5*sine;\
      sumzc+=$6*cosine;\
      sumzs+=$6*sine;\
      sumoc+=$7*cosine;\
      sumos+=$7*sine;\
      sumBc+=$8*cosine;\
      sumBs+=$8*sine;}\
    END{axc=sumxc/n;axs=sumxs/n;\
        ayc=sumyc/n;ays=sumys/n;\
        azc=sumzc/n;azs=sumzs/n;\
        aoc=sumoc/n;aos=sumos/n;\
        aBc=sumBc/n;aBs=sumBs/n;\
        rx=2*sqrt(axc*axc+axs*axs);if(rx==0)rx=1e-13;\
        ry=2*sqrt(ayc*ayc+ays*ays);if(ry==0)ry=1e-13;\
        rz=2*sqrt(azc*azc+azs*azs);if(rz==0)rz=1e-13;\
        ro=2*sqrt(aoc*aoc+aos*aos);if(ro==0)ro=1e-13;\
        rB=2*sqrt(aBc*aBc+aBs*aBs);if(rB==0)rB=1e-10;\
        px=atan2(axs/rx,axc/rx);if(px==0)px=1e-10;\
        py=atan2(ays/ry,ayc/ry);if(py==0)py=1e-10;\
        pz=atan2(azs/rz,azc/rz);if(pz==0)pz=1e-10;\
        po=atan2(aos/ro,aoc/ro);if(po==0)po=1e-10;\
        pB=atan2(aBs/rB,aBc/rB);if(pB==0)pB=1e-10;\
        printf("a_dx%d%d%d = %g\n",h,k,l,rx*mult);\
        printf("a_dy%d%d%d = %g\n",h,k,l,ry*mult);\
        printf("a_dz%d%d%d = %g\n",h,k,l,rz*mult);\
        printf("a_do%d%d%d = %g\n",h,k,l,ro*mult);\
        printf("a_dB%d%d%d = %g\n",h,k,l,rB);\
        if(! h000)printf("phi_dx%d%d%d = %g\n",h,k,l,px);\
        if(! h000)printf("phi_dy%d%d%d = %g\n",h,k,l,py);\
        if(! h000)printf("phi_dz%d%d%d = %g\n",h,k,l,pz);\
        if(! h000)printf("phi_do%d%d%d = %g\n",h,k,l,po);\
        if(! h000)printf("phi_dB%d%d%d = %g\n",h,k,l,pB);\
        }' |\
tee -a fitparams_init.gnuplot >> $logfile

end
echo "done with slow-ft"

#========================================================================
#    combine with previous values
cat fitparams_init.gnuplot fitparams.gnuplot |\
awk '$3!="new"{v[$1]=$3} {print $1,"=",v[$1]}' |\
awk '! seen[$1]{print;++seen[$1]}' |\
cat >! ${tempfile}
mv ${tempfile} fitparams.gnuplot


after_rejects:
#============================
#    re-entry point after editing fitparams or shift data

# clean up fitparams so that values don't get too small to fit?
cat fitparams.gnuplot |\
awk '$3+0!=0 && $3!="1e-10" && sqrt($3*$3)<1e-5{$3=1e-5} {print}' |\
cat >! ${tempfile}
set test = `diff fitparams.gnuplot $tempfile | egrep "^<" | wc -l`
if($test != 0) then
    echo "WARNING: augmenting $test very small values:"
    diff fitparams.gnuplot $tempfile | egrep "^<" 
endif
mv ${tempfile} fitparams.gnuplot


if(-e "$user_fitparams") then
    echo "using $user_fitparams"
    cp $user_fitparams fitparams.gnuplot
endif

# construct fit commands for gnuplot to match this param file
cat fitparams.gnuplot |\
awk '/^a/{hkl=substr($1,5);x=substr($1,4,1);a[x hkl]=$3}\
       /^phi/{hkl=substr($1,7);x=substr($1,6,1);p[x hkl]=$3}\
       END{for(xhkl in a) print xhkl,a[xhkl]+0,p[xhkl]+0}' |\
sort >! ${tempfile}apsi.txt 
#  pairwise table of fittable values

# zero value becomes a flag for do-not-fit
set n = 3
foreach x ( x y z o B )
    @ n = ( $n + 1 )

    # generate the fit command
    egrep "^$x" ${tempfile}apsi.txt |\
    awk -v n=$n 'NR==1{x=substr($1,1,1);\
            printf("fit real(d%s(x,y,z)) \"fitme\" using 1:2:3:%d:(1) via ",x,n);first=1}\
        {xhkl=$1}\
        $2!=0{\
            if(! first)printf(",");first=0;\
            printf("a_d%s",xhkl);}\
        $3!=0{\
            if(! first)printf(",");first=0;\
            printf("phi_d%s",xhkl);}\
        END{print ""}' |\
    tee fit${x}.gnuplot >> $logfile
end

# construct bending functions from these coefficients
rm -f func.gnuplot
foreach x ( x y z o B )
    set mult = $fitscale
    if( "$x" == "B" ) set mult = 1

    # define the function to fit in gnuplot5
    egrep "^$x" ${tempfile}apsi.txt |\
    awk -v mult=$mult 'NR==1{x=substr($1,1,1);\
            print "d"x"(x,y,z)=real(\\"}\
        {h=substr($1,2,1);k=substr($1,3,1);l=substr($1,4,1);\
            printf(" +a_d%s%d%d%d/%s",x,h,k,l,mult);\
            if($3+0==0){print "\\";next}\
            printf("*exp(2*pi*i*(%d*x+%d*y+%d*z+phi_d%s%d%d%d))\\\n",h,k,l,x,h,k,l);\
        } END{print ")"}' |\
    tee -a func.gnuplot >> $logfile
end



#========================================================================
#    this is the part that needs gnuplot5

echo "launching fitting jobs..."
onintr killfit
set jobsrunning
foreach x ( x y z o B )

    egrep "_d${x}" fitparams.gnuplot >! fitparams_${x}.gnuplot

    # make sure there is something to do?
    set test = `awk '/^a_/ && $NF>1e-6' fitparams_${x}.gnuplot | wc -l | awk '{print ( $1 <= 0)}'`
    if( "$test" == "1" ) continue
    set test = `echo $x $dimensions | awk '{for(i=2;i<=NF;++i)if($1==$i)print "1"}'`
    if( "$test" != "1" ) continue
    set test = `awk '{print ( $NF == "via" )}' fit${x}.gnuplot`
    if( "$test" == "1" ) continue

    gnuplot5 << EOF >&! fitrun_${x}.log &
    i = sqrt(-1)
    c(x) = column(x)
    set dummy x,y,z
    $FIT_LIMIT
    load 'func.gnuplot'
    load 'fitparams_${x}.gnuplot'
    load 'fit${x}.gnuplot'
    update 'fitparams_${x}.gnuplot'
EOF

end
wait
unset jobsrunning
killfit:
onintr
if($?jobsrunning) then
    killall $gnuplot5
    sleep 1
    killall -9 $gnuplot5
endif

# check that jobs worked?
foreach x ( x y z o B )
    if(! -e fitrun_${x}.log) continue
    set test = `grep -i singular fitrun_${x}.log | wc -l`
    if($test) then
        mv fitrun_${x}.log bad_fitrun_order${order}_${x}.log
        echo "WARNING: fit for $x is unstable. see bad_fitrun_order${order}_${x}.log"
echo "GOTHERE"
exit
    endif
end


#========================================================================
#    zero-fill anything we didnt fit
foreach x ( x y z B o )
    set test = `echo $x $dimensions | awk '{for(i=2;i<=NF;++i)if($1==$i)print "1"}'`
    if("$test" != "1") then
        awk '{$3=0;print}' fitparams_${x}.gnuplot >! ${tempfile}.txt
        mv ${tempfile}.txt fitparams_${x}.gnuplot
    endif
end

#========================================================================
#    combine them back together
echo -n "" >! fitparams.gnuplot
foreach x ( x y z o B )
    cat fitparams_${x}.gnuplot >> fitparams.gnuplot
end

# filter out ridiculously small amplitudes?
echo $CELL    $fitscale |\
cat - fitparams_[xyz].gnuplot |\
awk 'NR==1{a=$1;b=$2;c=$3;mult=$NF}\
   /^a_dx/{print $0,sqrt($NF*$NF)*a/mult}\
   /^a_dy/{print $0,sqrt($NF*$NF)*b/mult}\
   /^a_dz/{print $0,sqrt($NF*$NF)*c/mult}' |\
sort -k5g > /dev/null

# look for parameters that are far too correlated?
echo -n >! param_correlations.txt
foreach x ( x y z o B )
    if(! -e fitrun_${x}.log) continue

    awk '/correlation matrix/,""' fitrun_${x}.log |\
    awk 'NR==1{getline;next}\
       {col[++j]=$1;for(i=1;i<NF;++i) CC[i,j]=$(i+1)}\
       END{for(i in col)for(j in col) if(CC[i,j]!="" && i!=j)\
         print CC[i,j],col[i],col[j]}' |\
    cat >> param_correlations.txt
end
sort -gr param_correlations.txt |\
awk '$1>0.90 && $3~/^a/{print "eliminate",substr($3,3),$1,"CC with",$2}' |\
cat - fitparams.gnuplot |\
awk '/^eliminate/{++elim[$2];next}\
  {split($1,w,"_");xhkl=w[2]}\
  ! elim[xhkl]{print}' |\
cat >! fitparams_edited.gnuplot
diff fitparams.gnuplot fitparams_edited.gnuplot


echo -n "biggest coefficient: "
awk '/^a_/' fitparams_x.gnuplot | sort -k3gr | awk '{print;exit}'


bend_pdb:
#========================================================================
#    2nd round of gnuplot jobs to evaluate bend function for every atom

echo "applying smoothed shifts to $id1 "
ln -sf `pwd`/${id1}.frac ${tempfile}this.frac
if(! -e "$pdb2") then
    set BAD = "need destination cell!  please provide two PDB files."
    goto exit
endif
set CELL2 = `awk '/CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdb2`

set n = 1
foreach x ( x y z B o )
  @ n = ( $n + 1 )

  gnuplot << EOF &
frac = $frac
i = sqrt(-1)
c(x) = column(x)
load 'func.gnuplot'
load 'fitparams.gnuplot'

set table "${tempfile}atom_${x}.txt"
plot '${tempfile}this.frac' using (c(0)+1):(c($n)+frac*real(d${x}(c(2),c(3),c(4))))
EOF
end
wait

#========================================================================
#    stitch results together

echo "combining dimensions together..."
awk '$3=="i"{++n;print n,$2,"x"}' ${tempfile}atom_x.txt >! ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"y"}' ${tempfile}atom_y.txt >> ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"z"}' ${tempfile}atom_z.txt >> ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"o"}' ${tempfile}atom_o.txt >> ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"B"}' ${tempfile}atom_B.txt >> ${tempfile}newpos.txt

cat ${tempfile}newpos.txt ${tempfile}this.frac |\
awk '$3=="x"{x[$1]=$2;next}\
     $3=="y"{y[$1]=$2;next}\
     $3=="z"{z[$1]=$2;next}\
     $3=="o"{o[$1]=$2;next}\
     $3=="B"{B[$1]=$2;next}\
 {++n;an=$1;ID=substr($0,51);\
  printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%1s\n", \
       an, x[n], y[n], z[n], B[n], o[n], ID}' |\
cat >! ${tempfile}bent.frac

echo "converting back to orthogonal coordinates."
coordconv XYZIN ${tempfile}bent.frac \
         XYZOUT ${tempfile}bent.pdb << EOF-conv >> $logfile
CELL $CELL2
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv

# make sure it has the space group in the header
echo "SPACE $SG" | pdbset xyzin ${tempfile}bent.pdb xyzout bent${order}_${id1}.pdb >> $logfile
echo "transformed $id1 to match $id2 is bent${order}_${id1}.pdb"

# since pdbset and coordconv screw up names, transplant these from the original
# keeping only XYZ and B from the bending results
head -n 4 bent${order}_${id1}.pdb | egrep "^CRYST|^SCALE" >! ${tempfile}.pdb
cat $pdb1 |\
awk '/^ATOM|^HETAT/{++n;\
      print "REST",n,"|",substr($0,67);\
      print substr($0,1,30),"OLDNAME",n}' |\
cat bent${order}_${id1}.pdb - |\
awk '/^REST/{rest[$2]=substr($0,index($0,"|")+2)}\
    / OLDNAME /{n=$NF;print substr($0,1,30) xyzOB[n] rest[n];next}\
   /^ATOM|^HETAT/{++n;xyzOB[n]=substr($0,31,3*8+2*6)}' |\
cat >> ${tempfile}.pdb
mv ${tempfile}.pdb bent${order}_${id1}.pdb

echo "residual after shifts:"
if(! -e ${tempfile}remaining.pdb) cp $pdb2 ${tempfile}remaining.pdb
cat ${tempfile}remaining.pdb bent${order}_${id1}.pdb | egrep -v "HETAT" |\
 awk '$NF != "H"' |\
 rmsd | grep -v WARN

if("$geotest" == "true" || "$geotest" == "True") then
    echo "evaluating chemical geometry after distortion"
    refmac5 xyzin bent${order}_${id1}.pdb xyzout ${tempfile}junk.pdb << EOF >&! ${tempfile}geo.log
    refi type ideal
    ncyc 1
EOF
    cat ${tempfile}geo.log |\
    awk '/rmsBOND/{print $7,$8,$9,$10,$11;getline;\
           getline;print $7,$8,$9,$10,$11}' 
endif


#========================================================================
#     now reject outliers, if specified

if("$reject" != 0) then

    cat ${tempfile}remaining.pdb bent${order}_${id1}.pdb | egrep -v "HETAT" |\
     awk '$NF != "H"' |\
     rmsd -v debug=1 |\
     grep -v WARN |\
     tee ${tempfile}moves.txt |\
     awk '/moved/{print substr($0,23,8)}' |\
       sort -g |\
       awk '# compute median old fashioned way\
            {++n;v[n]=$1} \
            END{med=v[int(n/2)];\
              for(i=1;i<=n;++i){print sqrt((v[i]-med)**2),med}}' |\
       sort -g |\
       awk '# median absolute deviation\
            {++n;v[n]=$1;med=$2} \
            END{mad=v[int(n/2)];\
             print med,mad}' |\
     awk -v reject=$reject '{print $1+reject*$2+0.1}' >! ${tempfile}thresh.txt
    set threshold = `cat ${tempfile}thresh.txt`
    echo "deeming residuals > $threshold A as too big."
    rm -f ${tempfile}thresh.txt

    cat ${tempfile}moves.txt |\
    awk -v threshold=$threshold '/moved/ && substr($0,23,8)+0>threshold' |\
    cat >! ${tempfile}rejects.txt
    set nrejects = `cat ${tempfile}rejects.txt | wc -l`
    echo "rejecting $nrejects outliers > $reject * MAD from median"

    cat ${tempfile}rejects.txt ${tempfile}remaining.pdb |\
    awk '/moved/{++rejected[substr($0,1,15)];next} {id=substr($0,12,15)}\
       /^ATOM|^HETAT/ && ! rejected[id]{print}' |\
    cat >! ${tempfile}new.pdb
    mv ${tempfile}new.pdb ${tempfile}remaining.pdb

    egrep "^CRYST|^SCALE" bent${order}_${id1}.pdb >! ${tempfile}rejectme.pdb
    cat ${tempfile}rejects.txt ${tempfile}remaining.pdb |\
    awk '/moved/{++rejected[substr($0,1,15)];next} {id=substr($0,12,15)}\
       /^ATOM|^HETAT/ && rejected[id]{print}' |\
    cat >> ${tempfile}rejectme.pdb

    coordconv xyzin ${tempfile}rejectme.pdb xyzout ${tempfile}rejectme.frac << EOF > $logfile
    CELL $CELL2
    INPUT PDB
    OUTPUT FRAC
    END
EOF
    awk '{print $0,"REJECT"}' ${tempfile}rejectme.frac |\
    cat - ${id2}_P1.frac |\
    awk '{id=substr($0,61,15)} \
        $NF=="REJECT"{++rejected[id]}\
        ! rejected[id]{print}' |\
    cat >! ${tempfile}new.frac
    mv ${tempfile}new.frac ${id2}_P1.frac

    # get the shifts in x,y,z,o,B for every atom in common
    cat ${id1}_P1.frac ${id2}_P1.frac |\
    awk '{ID=substr($0,61,20);\
         x=substr($0,6,10);y=substr($0,16,10);z=substr($0,26,10);\
        o=substr($0,46,5);b=substr($0,36,10);++seen[ID]}\
        seen[ID]==1{X[ID]=x;Y[ID]=y;Z[ID]=z;Occ[ID]=o;Bfac[ID]=b;next}\
        seen[ID]==2{\
           print X[ID],Y[ID],Z[ID],x-X[ID],y-Y[ID],z-Z[ID],o-Occ[ID],b-Bfac[ID],"|" ID;\
           next}\
        {print "# error",ID,"seen",seen[ID],"times"}' |\
    grep -v HOH |\
    cat >! fitme
    # format xf yf zf  dxf dyf dzf do dB |ID

    if($nrejects > 1) then
        echo "fitting again"
        goto after_rejects
    endif
endif

if(! -e "$user_fitparams") then
    cp fitparams.gnuplot fitparams_${order}.gnuplot
    cp func.gnuplot func_${order}.gnuplot
else
    if($?REFIT) cp fitparams.gnuplot fitparams_refit.gnuplot
endif


if(! $?MAKE_DELTAMAPS) goto bend_map

echo ""
echo "making maps of bend function values: delta-[xyzr].map"

# need to make a map basis for this
echo "1 1 1 1 1" >! ${tempfile}dummy.hkl
f2mtz hklin ${tempfile}dummy.hkl hklout ${tempfile}dummy.mtz << EOF > $logfile
CELL $CELL2
SYMM 1
labou H K L F PHI
ctypo H H H F P
EOF
fft hklin ${tempfile}dummy.mtz mapout ${tempfile}dummy.map << EOF > $logfile
labin F1=F PHI=PHI
reso 4
EOF
mapmask mapin ${tempfile}dummy.map mapout ${tempfile}dumpme.map << EOF > $logfile
XYZLIM CELL
AXIS X Y Z
EOF

# convert map grid points into a pdb file
map2pdb.com ${tempfile}dumpme.map | tee ${tempfile}map2pdb.log > $logfile
set header = `awk '/^map header is/{print $4}' ${tempfile}map2pdb.log`

# get all stats
echo "go" | mapdump mapin ${tempfile}dumpme.map >! ${tempfile}mapdump
set GRID    = `awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump`
set LIMITS  = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${tempfile}mapdump`
set AXIS    = `awk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump | awk '! /[^ XYZ]/'`

set voxels = `awk '/Number of columns, rows, sections/{print $7*$8*$9; exit}' ${tempfile}mapdump`
set header = `ls -l ${tempfile}dumpme.map | awk -v voxels=$voxels '{print $5-4*voxels}'`
echo "map header is $header bytes"

#========================================================================
#    round of gnuplot jobs to evaluate bend function at every map grid point

foreach x ( x y z )

  gnuplot << EOF &
i = sqrt(-1)
c(x) = column(x)
load 'func.gnuplot'
load 'fitparams.gnuplot'

set table "${tempfile}delta_${x}.txt"
plot 'mapdump.frac' using (c(0)+1):(real(d${x}(c(2),c(3),c(4))))
EOF
end

gnuplot << EOF &
i = sqrt(-1)
c(x) = column(x)
load 'func.gnuplot'
load 'fitparams.gnuplot'

set table "${tempfile}delta_r.txt"
r(x,y,z) = sqrt(real(dx(c(2),c(3),c(4)))**2 +real(dy(c(2),c(3),c(4)))**2 +real(dz(c(2),c(3),c(4)))**2)
plot 'mapdump.frac' using (c(0)+1):(r(c(2),c(3),c(4)))
EOF

wait


foreach x ( x y z r )

echo "assembing delta-${x} map"
head -c ${header} ${tempfile}dumpme.map >! ${tempfile}.map
awk '$3=="i"{print $2}' ${tempfile}delta_${x}.txt |\
floatgen >> ${tempfile}.map

echo scale factor 1 |\
 mapmask mapin ${tempfile}.map mapout delta-${x}.map > $logfile

echo | mapdump mapin delta-${x}.map | grep dens

end









bend_map:
#========================================================================
#     go on to next order now, if there is no map


# no map and user-defined bend = nothing left to do
if(! -e "$mapfile" && -e "$user_fitparams") goto cleanup
if(! -e "$mapfile") continue



#========================================================================
#    apply bend function to map file

echo "simplifying $mapfile ..."
mapmask mapin $mapfile mapout ${tempfile}dumpme.map << EOF > $logfile
XYZLIM ASU
AXIS X Y Z
EOF

# convert map grid points into a pdb file
map2pdb.com ${tempfile}dumpme.map | tee ${tempfile}map2pdb.log
set header = `awk '/^map header is/{print $4}' ${tempfile}map2pdb.log`


# get all stats
echo "go" | mapdump mapin ${tempfile}dumpme.map >! ${tempfile}mapdump
set mapSG   = `awk '/Space-group/{print $NF; exit}' ${tempfile}mapdump`
set mapCELL = `awk '/Cell dimensions/{print $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump`
set GRID    = `awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump`
set LIMITS  = `awk '/Fast, medium, slow/{o[$(NF-2)]=1;o[$(NF-1)]=2;o[$NF]=3; print b[o["X"]],e[o["X"]], b[o["Y"]],e[o["Y"]], b[o["Z"]],e[o["Z"]]; exit} /Start and stop points/{b[1]=$(NF-5); e[1]=$(NF-4); b[2]=$(NF-3); e[2]=$(NF-2); b[3]=$(NF-1); e[3]=$NF}' ${tempfile}mapdump`
set AXIS    = `awk '/Fast, medium, slow/{print $(NF-2), $(NF-1), $NF}' ${tempfile}mapdump | awk '! /[^ XYZ]/'`

set voxels = `awk '/Number of columns, rows, sections/{print $7*$8*$9; exit}' ${tempfile}mapdump`
set header = `ls -l ${tempfile}dumpme.map | awk -v voxels=$voxels '{print $5-4*voxels}'`
echo "map header is $header bytes"

# note: new map will need to have the same cell as 2nd PDB file!


#========================================================================
#    round of gnuplot jobs to evaluate bend function at every map grid point
#    note the negative sign below because we are first bending coordinates
#    and then evaluating the map.

set n = 1
foreach x ( x y z )
  @ n = ( $n + 1 )

  gnuplot << EOF &
i = sqrt(-1)
c(x) = column(x)
load 'func.gnuplot'
load 'fitparams.gnuplot'

set table "${tempfile}atom_${x}.txt"
plot 'mapdump.frac' using (c(0)+1):(c($n)-real(d${x}(c(2),c(3),c(4))))
EOF
end
wait

awk '$3=="i"{++n;print n,$2,"x"}' ${tempfile}atom_x.txt >! ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"y"}' ${tempfile}atom_y.txt >> ${tempfile}newpos.txt
awk '$3=="i"{++n;print n,$2,"z"}' ${tempfile}atom_z.txt >> ${tempfile}newpos.txt

cat ${tempfile}newpos.txt mapdump.frac |\
awk '$3=="x"{x[$1]=$2;next}\
     $3=="y"{y[$1]=$2;next}\
     $3=="z"{z[$1]=$2;next}\
 {++n;an=$1;ID=substr($0,51);\
  printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%1s\n", \
       an, x[n], y[n], z[n], B[n], o[n], ID}' |\
cat >! ${tempfile}probe.frac

echo "preparing to peek with mapman"
coordconv XYZIN ${tempfile}probe.frac \
         XYZOUT ${tempfile}probe.pdb << EOF-conv >> $logfile
CELL $mapCELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv

########################################################
#  note:
#  the bend function moves points in space from the bendme
#  frame to the reference frame.
#  we are now going to sample the map at these moved points
#  then move those values back onto the map grid
#
#  so, $mapfile had better be in the reference frame!
#

# extend the map a little for cubic splines
cat ${tempfile}probe.frac |\
awk '$2<minx{minx=$2} $2>maxx{maxx=$2}\
     $3<miny{miny=$3} $3>maxy{maxy=$3}\
     $4<minz{minz=$4} $4>maxz{maxz=$4}\
     END{print minx+0,maxx+0,miny+0,maxy+0,minz+0,maxz+0}' |\
cat >! ${tempfile}fraclimits.txt
set fracLIMITS = `cat ${tempfile}fraclimits.txt`
set gridLIMITS = `echo $fracLIMITS $GRID | awk '{gx=$7;gy=$8;gz=$9;print $1*gx,$2*gx,$3*gy,$4*gy,$5*gz,$6*gz}'`
set extLIMITS = `echo $gridLIMITS 5 | awk '{print int($1-$NF),int($2+$NF),int($3-$NF),int($4+$NF),int($5-$NF),int($6+$NF)}'`

echo "extending map a bit for cubic spline interpolation in mapman"
mapmask mapin ${tempfile}dumpme.map mapout ${tempfile}expanded.map << EOF >> $logfile
XYZLIM $extLIMITS
EOF


echo "getting density values from expanded.map using mapman"
setenv MAPSIZE `ls -l ${tempfile}expanded.map | awk '{printf "%d", $5/3.5}'`
$mapman -b mapsize $MAPSIZE << mapman-end >&! ${tempfile}mappeek.log
read map1 ${tempfile}expanded.map ccp4
peek value map1 ${tempfile}probe.pdb /dev/null spline ;
quit
mapman-end
set test = `grep -l PEEK-A-BOO ${tempfile}mappeek.log | wc -l`
if( ! $test ) then
    set BAD = "MAPMAN failed "
    goto exit
endif


# need to create a map header with the right cell
echo "1 1 1 1 1" >! ${tempfile}dummy.hkl
f2mtz hklin ${tempfile}dummy.hkl hklout ${tempfile}dummy.mtz << EOF > $logfile
CELL $CELL2
SYMM $mapSG
labou H K L F PHI
ctypo H H H F P
EOF
fft hklin ${tempfile}dummy.mtz mapout ${tempfile}dummy.map << EOF > $logfile
labin F1=F PHI=PHI
GRID $GRID
EOF
mapmask mapin ${tempfile}dummy.map mapout ${tempfile}header.map << EOF > $logfile
XYZLIM $LIMITS
AXIS $AXIS
EOF
echo "go" | mapdump mapin ${tempfile}header.map >! ${tempfile}mapdump 
set voxels = `awk '/Number of columns, rows, sections/{print $7*$8*$9; exit}' ${tempfile}mapdump`
set header = `ls -l ${tempfile}dumpme.map | awk -v voxels=$voxels '{print $5-4*voxels}'`
echo "new map header is $header bytes"


echo "assembing final bent map"
head -c ${header} ${tempfile}header.map >! ${tempfile}new.map
awk '/PEEK-A-BOO/{print $NF}' ${tempfile}mappeek.log |\
floatgen >> ${tempfile}new.map

# test to see if map interpolation worked
set test = `ls -l ${tempfile}new.map | awk '{print ( $5 > 100 )}'`
if("$test" != "1") then
    stat ${tempfile}expanded.map
    set BAD = "unable to interpolate ${tempfile}expanded.map. can $mapman read the file ? "
    goto exit
endif

# force rewrite of header stats
echo scale factor 1 |\
mapmask mapin ${tempfile}new.map mapout bent${order}.map > $logfile
ls -l $mapfile bent${order}.map


set aligner = "$pdb1"
if(-e "$user_fitparams") set aligner = "$user_fitparams"
echo "bent${order}.map is now a version of $mapfile aligned with $aligner"

# we are not looping if user specified params
if(-e "$user_fitparams") goto cleanup

#========================================================================
#    end of big loop over orders
end

cleanup:

set deltaT = `date +%s | awk -v s=$starttime '{print $1-s}'`
echo ""
echo "=================  all done at $deltaT s"
echo ""

if(! $?debug) then
    rm -f ${tempfile}* >& /dev/null
endif

if(-e busy) rm busy

exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

exit

































###############################################################################################
###############################################################################################
###############################################################################################
#
#	Bringing our own lunch: option to deploy/compile most of the 3rd-part programs
#
###############################################################################################
###############################################################################################
###############################################################################################





deploy_scripts:

setenv DEPLOYED


cat << EOF-script >! rmsd
#! `which awk` -f
#! /usr/bin/awk -f
#
#   Calculate RMSD of atoms with the same name in two PDB files		James Holton 9-5-11
#
#	The PDB feild:
#           |<--- here -->|
#ATOM      1  N   ALA A 327      40.574  34.523  43.012  1.00 34.04
#
#	is used to determine if two atoms are a "pair"
#
BEGIN {
if(! atom) atom = "CA"
maxXYZ = maxdB = maxdO = 0
max_atom_XYZ = max_atom_dB = max_atom_dO = 0
maxXYZ_ID = maxdB_ID = maxdO_ID = "-"
max_atom_XYZ_ID = max_atom_dB_ID = max_atom_dO_ID = "-"
}

/^ATOM|^HETATM/{
    # read in values (watching for duplications)
    ID = substr(\$0,12,15)
    ++count[ID]

    if(count[ID] == 1)
    {
	# not initialized yet
	X[ID] = substr(\$0, 31, 8)+0
	Y[ID] = substr(\$0, 39, 8)+0
	Z[ID] = substr(\$0, 47, 8)+0
	O[ID] = substr(\$0, 55, 6)+0
	B[ID] = substr(\$0, 61, 6)+0
    }
    
    if(count[ID] == 2)
    {
	++pairs
    
	# seen this before, subtract values
	dX     = X[ID] - substr(\$0, 31, 8)
	dY     = Y[ID] - substr(\$0, 39, 8)
	dZ     = Z[ID] - substr(\$0, 47, 8)
	dO[ID] = O[ID] - substr(\$0, 55, 6)
	dB[ID] = B[ID] - substr(\$0, 61, 6)
	
	# get drift (and add up squares of drifts)
	sqrD   = dX*dX + dY*dY + dZ*dZ
	dXYZ[ID] = sqrt(sqrD)

	# remember maximum shifts
	if(dXYZ[ID] > maxXYZ) {maxXYZ  = dXYZ[ID]; maxXYZ_ID = ID }
	if(dO[ID]*dO[ID] > maxdO*maxdO) {maxdO = dO[ID]; maxdO_ID = ID }
	if(dB[ID]*dB[ID] > maxdB*maxdB) {maxdB = dB[ID]; maxdB_ID = ID }

	# maintain mean-square sums	
	sumXYZ += sqrD
	sumO   += dO[ID]*dO[ID]
	sumB   += dB[ID]*dB[ID]

	# separate stats for special atom type
	split(ID,word)
	if(word[1] == atom)
	{
	    ++atom_pairs

	    # maintain separate mean-square sums
	    sum_atom_XYZ += sqrD
	    sum_atom_O   += dO[ID]*dO[ID]
	    sum_atom_B   += dB[ID]*dB[ID]
	    
	    # remember maximum drifts too
	    if(dXYZ[ID] > max_atom_XYZ) {max_atom_XYZ  = dXYZ[ID]; max_atom_XYZ_ID = ID }
	    if(dO[ID]*dO[ID] > max_atom_dO*max_atom_dO) {max_atom_dO = dO[ID]; max_atom_dO_ID = ID }	    
	    if(dB[ID]*dB[ID] > max_atom_dB*max_atom_dB) {max_atom_dB = dB[ID]; max_atom_dB_ID = ID }	    
	}
	# debug output
	if(debug)
	{
	    printf("%s moved %8.4f (XYZ) %6.2f (occ) %6.2f (B) at %s\\n", ID, dXYZ[ID],dO[ID], dB[ID], "cen_x " substr(\$0,31,24))
	}	
    }

    if(count[ID] > 2)
    {
	print "WARNING: " ID " appeared more than twice! "
    }
}


END{
    
    if(pairs+0 == 0) 
    {
	print "no atom pairs found"
	exit
    }
    rmsXYZ = sqrt(sumXYZ/pairs)
    rmsO = sqrt(sumO/pairs)
    rmsB = sqrt(sumB/pairs)
    if(atom_pairs+0 != 0)
    {
	rms_atom_XYZ = sqrt(sum_atom_XYZ/atom_pairs)
	rms_atom_O = sqrt(sum_atom_O/atom_pairs)
	rms_atom_B = sqrt(sum_atom_B/atom_pairs)
    }
    

    if(! xlog) 
    {
	print pairs " atom pairs found"
	print "RMSD("atom" )= " rms_atom_XYZ " ("atom_pairs, atom " pairs)"
	print "RMSD(all)= " rmsXYZ " ("pairs" atom pairs)"
	print "RMSD(Bfac)= " rmsB
	if(maxdO_ID != "-") print "RMSD(occ)= " rmsO

	print "MAXD(all)= " maxXYZ "\\tfor " maxXYZ_ID
	print "MAXD(Bfac)= " maxdB "\\tfor " maxdB_ID
	if(maxdO_ID != "-") print "MAXD(occ)= " maxdO "\\tfor " maxdO_ID
	
	# final check for orphan atoms 
	for(id in count)
	{
	    if(count[id]<2) print "WARNING: " id " only found once"
	}
    }
    else
    {
#	printf "%10.8f %10.8f %10.5f %8.4f %10.4f %8.3f %8.3f    %s %s %s\\n", rms_atom_XYZ, rmsXYZ, rmsO, rmsB, maxXYZ, maxO,maxdB, maxXYZ_ID,maxdB_ID,maxdO_ID
	printf "%10.8f %10.8f %10.5f %8.4f %10.4f %8.3f %8.3f\\n", rms_atom_XYZ, rmsXYZ, rmsO, rmsB, maxXYZ, maxdO,maxdB
    }
}
EOF-script
chmod a+x rmsd




















cat << EOF-script >! map2pdb.com
#! `which tcsh` -f
#! /bin/tcsh -f
#
#	convert a CCP4 map into a PDB file with "rho" encoded as the B factor
#
#	-James Holton 9-20-17
#
#
set mapfile = "\$1"

if(! -e "\$mapfile") then
    cat << EOF
usage: \$0 mapfile.map

where: 
mapfile.map   is a CCP4 format map file

output will be a PDB file with one atom per map voxel.  Electron
density will be encoded in the B factors, as well as off the end of
the cannonical 80-byte record.

EOF
    exit 9
endif

set tempfile = \${CCP4_SCR}/map2pdb\$\$
mkdir -p \${CCP4_SCR} >& /dev/null

mapmask mapin \$mapfile mapout \${tempfile}dumpme.map << EOF | tee \${tempfile}mapdump.log | grep dens
axis X Y Z
EOF
cat \${tempfile}mapdump.log |\\
awk '/Grid sampling on x, y, z/{gx=\$8;gy=\$9;gz=\$10}\\
     /Cell dimensions /{xs=\$4/gx;ys=\$5/gy;zs=\$6/gz}\\
     /Maximum density /{max=\$NF}\\
     /Minimum density /{min=\$NF}\\
     /Start and stop points on cols/{cs=\$9+0;rs=\$11+0;ss=\$13+0}\\
     /Number of columns, rows, sections/{nc=\$7;nr=\$8;ns=\$9}\\
 END{print xs,ys,zs,cs,rs,ss,nc,nr,ns,gx,gy,gz,max,min}' >! \${tempfile}mapstuff.txt

awk '{print "grid points every:",\$1,\$2,\$3}' \${tempfile}mapstuff.txt


set CELL = \`awk '/Cell dimensions /{print \$4,\$5,\$6,\$7,\$8,\$9;exit}' \${tempfile}mapdump.log\`
set SG   = \`awk '/Space-group /{print \$NF;exit}' \${tempfile}mapdump.log\`


set size = \`awk '{print 4*\$7*\$8*\$9}' \${tempfile}mapstuff.txt\`
set head = \`ls -l \${tempfile}dumpme.map | awk -v size=\$size '{print \$5-size}'\`
echo "map header is \$head bytes"
od -vf -w4 -j \$head \${tempfile}dumpme.map | awk 'NF==2{print \$2}'|\\
cat \${tempfile}mapstuff.txt - |\\
awk 'NR==1{cs=\$4;rs=\$5;ss=\$6;nc=\$7;nr=\$8;ns=\$9;gx=\$10;gy=\$11;gz=\$12;max=\$13;min=\$14;chain=65;\\
   if(max+0<-min)max=-min;\\
   if(max+0==0)max=0.01;\\
   next}\\
 {printf("%5d%10.5f%10.5f%10.5f%10.5f%5.2f    1%10d E   RHO %c    %15.12f\\n",\\
   ++a,(c+cs)/gx,(r+rs)/gy,(s+ss)/gz,\$1,\$1/max,++n,chain,\$1)}\\
  a>=99999{a=0}\\
  n>=9999 && chain==32{n=0;chain=65}\\
  n>=9999 && chain==90{n=0;chain=97}\\
  n>=9999 && chain==126{n=0;chain=33} n>=9999{n=0;++chain}\\
  {++c} c>=nc{c=0;++r} r>=nr{r=0;++s}' |\\
tee mapdump.frac | awk -v size=\$size '{++n;pct=int(n/size*400)} ! seen[pct]{++seen[pct];printf("%d%% done\\r",pct)}'
echo ""
echo "converting mapdump.frac into mapdump.pdb"
coordconv XYZIN mapdump.frac XYZOUT mapdump.pdb << EOF-conv > /dev/null
CELL \$CELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv

echo "decorating..."
cat mapdump.frac mapdump.pdb |\\
awk 'substr(\$0,71,3)=="RHO"{++n;v[n]=substr(\$0,77);next}\\
 ! /^ATOM/{print}\\
   /^ATOM/{++m;print \$0,v[m]}' |\\
cat >! \${tempfile}mapdump.pdb
mv \${tempfile}mapdump.pdb mapdump.pdb

if(! -s mapdump.pdb) exit 9

rm -f \${tempfile}dumpme.map \${tempfile}mapstuff.txt \${tempfile}mapdump.log \${tempfile}pdbset.log
rm -f \${tempfile}mapdump.pdb
EOF-script
chmod a+x map2pdb.com




















cat << EOF-origins.com >! origins.com
#! `which tcsh` -f
#! /bin/tcsh -f
#
#    origins.com                    - James Holton 9-15-17
#
#    script for translating one PDB to a number of alternative origins
#       and indexing conventions
#    and checking if the symmetry-expanded atoms are close to those of
#    another PDB.
#
#    Atoms with the same name and resiude number are compared
#
set awk = awk
\$awk 'BEGIN{print}' >& /dev/null
if(\$status) set awk = gawk
alias awk \$awk


set SG        = ""
set right_pdb = ""
set wrong_pdb = ""
set outfile = neworigin.pdb


mkdir -p \${CCP4_SCR} >&! /dev/null
set tempfile = \${CCP4_SCR}/origins_temp\$\$
if(\$?debug) set tempfile = ./origins_temp
################################################################################
goto Setup
Help:
cat << EOF

usage: \$0 right_origin.pdb wrong_origin.pdb \$SG [correlate] [nochains]

where: right_origin.pdb    is relative to your "desired" origin
       wrong_origin.pdb    is relative to another origin
       \$SG        is your space group

wrong_origin.pdb will be moved to every possible origin in \${SG}
and then checked to see if the moved, and symmetry-expanded atoms
line up with the ones in right_origin.pdb.

The results of the origin choice that give the best agreement
to the reference pdb will be copied to \$outfile

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
if(\$?CORRELATE) set SCORE = " CC "

# get unit cell from either pdb
set CELL = \`awk '/^CRYST/{print \$2,\$3,\$4,\$5,\$6,\$7}' \$right_pdb \$wrong_pdb | head -1\`
if("\$CELL" == "") then
    echo "ERROR: no CRYST card in \$right_pdb"
    goto Help
endif
# just in case these are different?
set right_cell = \`awk '/^CRYST/{print \$2,\$3,\$4,\$5,\$6,\$7}' \$right_pdb | head -1\`
set wrong_cell = \`awk '/^CRYST/{print \$2,\$3,\$4,\$5,\$6,\$7}' \$wrong_pdb | head -1\`
if("\$right_cell" == "") set right_cell = "\$CELL"
if("\$wrong_cell" == "") set wrong_cell = "\$right_cell"
# make the wrong_cell the default, so we are moving in allowed space
set CELL = ( \$wrong_cell )

# try to get space group?
if("\$SG" == "") then
    set pdbSG = \`awk '/^CRYST/{print substr(\$0,56,12)}' \$right_pdb \$wrong_pdb | head -1\`
    set SG = \`awk -v pdbSG="\$pdbSG" -F "[\\047]" 'pdbSG==\$2{print;exit}' \${CLIBD}/symop.lib | awk '{print \$4}'\`
endif
if("\$SG" == "") then
    # hmm, throw an error or go on? ...
    set SG = P1
endif
set right_latt = \`echo \$SG | awk '{print substr(\$0,1,1)}'\`

# check for H3/R3 weirdness
set test = \`echo \$CELL | awk '{print ( \$4+0==90 && \$5+0==90 && \$6+0==120 )}'\`
if("\$right_latt" == "R" && "\$test" == "1" ) then
    # probably should use hexagonal system
    set SG = \`echo \$SG | awk '{gsub("R","H");print}'\`
endif
if("\$right_latt" == "H" && "\$test" == "0" ) then
    # probably should use rhombohedral system
    set SG = \`echo \$SG | awk '{gsub("H","R");print}'\`
endif

# check if wrong pdb has different info
set wrong_pdbSG = \`awk '/^CRYST/{print substr(\$0,56,12)}' \$right_pdb \$wrong_pdb | head -1\`
set wrong_SG = \`awk -v pdbSG="\$wrong_pdbSG" -F "[\\047]" 'pdbSG==\$2{print;exit}' \${CLIBD}/symop.lib | awk '{print \$4}'\`
if("\$wrong_SG" == "") set wrong_SG = \$SG
set wrong_latt = \`echo \$wrong_SG | awk '{print substr(\$0,1,1)}'\`
set test = \`echo \$wrong_cell | awk '{print ( \$4+0==90 && \$5+0==90 && \$6+0==120 )}'\`

if("\$wrong_latt" == "R" && "\$test" == "1" ) then
    # probably should use hexagonal system
    set wrong_SG = \`echo \$wrong_SG | awk '{gsub("R","H");print}'\`
endif
if("\$wrong_latt" == "H" && "\$test" == "0" ) then
    # probably should use rhombohedral system
    set wrong_SG = \`echo \$wrong_SG | awk '{gsub("H","R");print}'\`
endif
set wrong_latt = \`echo \$wrong_SG | awk '{print substr(\$0,1,1)}'\`



# space group for correlation calculations
set CCSG = \$SG

# preemtive cleanup
rm -f \${tempfile}_right[0-9][0-9][0-9].pdb >& /dev/null
rm -f \${tempfile}_wrong[0-9][0-9][0-9].pdb >& /dev/null


set reindexings = X,Y,Z

if(! \$?ALTINDEX) goto otherhand

    echo "finding all possible re-indexing operators..."
    echo "from: \$wrong_cell latt= \$wrong_latt"
    echo "  to: \$right_cell"
    # find all possible reindexing operations
    othercell << EOF | tee \${tempfile}othercell.log > /dev/null
\$wrong_cell \$wrong_latt
\$right_cell
EOF
    set reciprocal_symops = \`awk '/close to target:/{getline;while(/]\$/){gsub("[][]","");print \$NF;getline}}' \${tempfile}othercell.log | sort -u\`
    # gather all symmetry operators from all possible space groups for this cell
    cat \${tempfile}othercell.log |\\
    awk '\$1~/^</ && \$NF~/>\$/' | sort -u |\\
    awk '{split(\$0,w,"[<>]");print "SPACEGROUP",w[2]}' |\\
    cat - \${CLIBD}/symop.lib |\\
    awk '/^SPACEGROUP/{pdbSG=substr(\$0,12);\\
              #print "SG:",pdbSG;\\
              ++seen[pdbSG];next}\\
         /[XYZ]/ && ! /[PCIFRH]/ && p==1 {print \$1 \$2 \$3 \$4}\\
        \$5 ~ /^PG/{p=0} {split(\$0,w,"\\047");pdbSG=w[2];if(seen[pdbSG])p=1}' |\\
    sort -u >! \${tempfile}all_symops

    # and eliminate ones that are already in the space group to be used
    cat \${CLIBD}/symop.lib |\\
    awk -v SG=\$CCSG '/[XYZ]/ && ! /[PCIFRH]/ && p==1 {print \$1 \$2 \$3 \$4}\\
        \$5 ~ /^PG/{p=0} \$4 == SG{p=1}' |\\
    awk '{print "GOT",\$0}' >! \${tempfile}SG_symops

    echo "+X,+Y,+Z" >! \${tempfile}reindexings.txt
    cat \${tempfile}SG_symops \${tempfile}all_symops |\\
    awk '/^GOT/{++got[substr(\$0,4)];next}\\
        ! got[\$0]{print}' |\\
    awk '{gsub(" ","");print}' |\\
    awk -F "," '{for(i=1;i<=NF;++i){print \$i};print ""}' |\\
    awk '/^[^-]/{\$0="+"\$0} {print}' |\\
    awk 'NF==0{print substr(op,2);op=""} NF>0{op=op"," \$0}' |\\
    cat >> \${tempfile}reindexings.txt

    # get list of all alternative indexings to try
    set reindexings = \`cat \${tempfile}reindexings.txt\`
    echo "\$#reindexings alternative indexing operators."
    rm -f \${tempfile}symops
    rm -f \${tempfile}all_symops
    rm -f \${tempfile}SG_symops
    rm -f \${tempfile}othercell.log


# filter out equivalent operators?
echo "filtering out redundant re-indexing operators"

# do a coarser map grid sampling than normal (~ 3A resolution)
set reso = 3
set fakeCELL = \`echo \$wrong_cell \$reso | awk '{print \$1/\$NF/2,\$2/\$NF/2,\$3/\$NF/2,\$4,\$5,\$6}'\`
set BADD     = \`echo \$reso | awk '{printf "%d", 79*(\$1/3)^2}'\`
echo "CELL \$fakeCELL" | pdbset xyzin \$wrong_pdb xyzout \${tempfile}.pdb > /dev/null

# make an all-carbon version of this model
cat \${tempfile}.pdb |\\
awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
cat >! \${tempfile}mapme.pdb
sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}right.map << EOF-sfall > /dev/null
MODE ATMMAP
CELL \$fakeCELL
SYMM \$CCSG
EOF-sfall
# recover coarse grid spacing 
set coarseGRID = \`echo "GO" | mapdump MAPIN \${tempfile}right.map | awk '/Grid sampling/{print \$(NF-2), \$(NF-1), \$NF; exit}'\`
rm -f \${tempfile}.pdb >& /dev/null

# now that we have a grid, make a safe version of whole file

echo "CELL \$wrong_cell" | pdbset xyzin \$wrong_pdb xyzout \${tempfile}.pdb > /dev/null
cat \${tempfile}.pdb |\\
awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
cat >! \${tempfile}reindexme.pdb

set n = 0
foreach reindexing ( \$reindexings )
    @ n = ( \$n + 1 )
    #echo "applying \$reindexing to \$wrong_pdb"
    pdbset xyzin \${tempfile}reindexme.pdb xyzout \${tempfile}reindexed.pdb << EOF > /dev/null
    symgen \$reindexing
    BFACTOR 80
EOF

    sfall xyzin \${tempfile}reindexed.pdb mapout \${tempfile}test_\${n}.map << EOF > /dev/null
    mode atmmap
    resolution 3
    GRID \$coarseGRID
    symm \$SG
EOF
end
echo "maps done"

# now look for maps that are prefectly correlated, these indicate one of the ops is redundant
set redundant = ""
rm -f CC_vs_ij.log
foreach i ( \`seq 1 \$n\` )
  @ k = ( \$i + 1 )
  foreach j ( \`seq \$k \$n\` )
    set test = \`echo \$i \$j \$redundant | awk '{for(k=3;k<=NF;++k){if(\$1==\$k || \$2==\$k){print 1;exit}}}'\`
    if("\$test" == "1") continue
    echo "correlate section" |\\
    overlapmap mapin1 \${tempfile}test_\${i}.map \\
            mapin2 \${tempfile}test_\${j}.map mapout /dev/null |\\
    awk '/Total correlation/{print \$NF}' >! \${tempfile}correlation
    set CC = \`awk '{print \$1}' \${tempfile}correlation\`
    #echo "\$CC \$i \$j" | tee -a CC_vs_ij.log
    set test = \`echo \$CC | awk '{print ( \$1 > 0.95 )}'\`
    if(\$test) then
       # the jth operator is redundant with teh ith operator
       set redundant = ( \$redundant \$j )
    endif
  end
end
foreach dupe ( \$redundant )
    set reindexings[\$dupe] = ""
end
set reindexings = ( \$reindexings )
echo "\$#reindexings re-indexing operators remain"
rm -f \${tempfile}correlation >& /dev/null
rm -f \${tempfile}test_*.map >& /dev/null
rm -f CC_vs_ij.log >& /dev/null


otherhand:
if(\$?OTHERHAND) then
    echo " \$reindexings " |\\
    awk '{for(i=1;i<=NF;++i)print \$i}' |\\
    awk '{print;\\
          gsub("+","q");gsub("-","+");gsub("q","-");\\
          print}' |\\
    cat >! \${tempfile}newops.txt
    set reindexings = \`cat \${tempfile}newops.txt\`
    rm -f \${tempfile}newops.txt
endif


if (! \$?BYFILE) then
    # break up "right" file into its respective chains
    cat \$right_pdb |\\
    awk '/^ATOM|^HETAT/{resnum=substr(\$0, 23, 4)+0;segid = substr(\$0, 22, 1); \\
         if(last_segid=="")last_segid=segid;\\
             if(resnum<last_resnum || (segid!=last_segid)) print "BREAK";\\
                last_resnum=resnum;last_segid=segid} {print}' |\\
    awk -v tempfile=\${tempfile}_right 'BEGIN{chain="001"} /^ATOM|^HETAT/{segid = substr(\$0, 22, 1)}\\
     /^BREAK/{print "REMARK chain", segid > tempfile chain ".pdb";\\
              chain=sprintf("%03d", chain+1)}\\
      /^ATOM|^HETATM/{print "ATOM  "substr(\$0,7,66) > tempfile chain ".pdb"}\\
          END{print "REMARK chain", segid > tempfile chain ".pdb"}'
    # chains are now separated into files
    # \${tempfile}_right###.pdb

    # break up "wrong" file into its respective chains
    cat \$wrong_pdb |\\
    awk '/^ATOM|^HETAT/{resnum=substr(\$0, 23, 4)+0;segid = substr(\$0, 22, 1); \\
         if(last_segid=="")last_segid=segid;\\
             if(resnum<last_resnum || (segid!=last_segid)) print "BREAK";\\
                last_resnum=resnum;last_segid=segid} {print}' |\\
    awk -v tempfile=\${tempfile}_wrong 'BEGIN{chain="001"} /^ATOM|^HETAT/{segid = substr(\$0, 22, 1)}\\
     /^BREAK/{print "REMARK chain", segid > tempfile chain ".pdb";\\
              chain=sprintf("%03d", chain+1)}\\
      /^ATOM|^HETAT/{print "ATOM  "substr(\$0,7,66) > tempfile chain ".pdb"}\\
          END{print "REMARK chain", segid > tempfile chain ".pdb"}'
    # chains are now separated into files
    # \${tempfile}_wrong###.pdb

else
    echo "REMARK CHAIN _" >! \${tempfile}_right001.pdb
    cat \$right_pdb |\\
    awk '/^ATOM|^HETAT/{print "ATOM  "substr(\$0,7,66)}' |\\
    cat >> \${tempfile}_right001.pdb
    echo "REMARK CHAIN _" >> \${tempfile}_right001.pdb

    echo "REMARK CHAIN _" >! \${tempfile}_wrong001.pdb
    cat \$wrong_pdb |\\
    awk '/^ATOM|^HETAT/{print "ATOM  "substr(\$0,7,66)}' |\\
    cat >> \${tempfile}_wrong001.pdb
    echo "REMARK CHAIN _" >> \${tempfile}_wrong001.pdb
    
endif

# add the proper headers
foreach file ( \${tempfile}_right[0-9][0-9][0-9].pdb \${tempfile}_wrong[0-9][0-9][0-9].pdb )

    echo "END" >> \$file

    # calculate the center of mass
    pdbset xyzin \$file xyzout \${tempfile}.pdb << EOF >! \${tempfile}.log
    CELL \$CELL
    CHAIN " "
    COM
EOF
    egrep "^REMARK" \$file >! \${tempfile}new.pdb
    egrep -v "REMARK" \${tempfile}.pdb >> \${tempfile}new.pdb
    mv \${tempfile}new.pdb \${tempfile}.pdb
    
    # make a pdb file of an atom at the COM of this file
    awk '/^CRYST/ || /^SCALE/' \${tempfile}.pdb >! \${tempfile}COM.pdb
    cat \${tempfile}.log |\\
    awk '\$1=="Center"{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\\n", \$4, \$5, \$6;exit}' |\\
    cat >> \${tempfile}COM.pdb
    echo "END" >> \${tempfile}COM.pdb
    
    # get fractional version of the COM too
    coordconv XYZIN \${tempfile}COM.pdb XYZOUT \${tempfile}.xyz << EOF >& /dev/null
    CELL \$CELL
    INPUT PDB
    OUTPUT FRAC
EOF

    # put this center as a remark in the PDB file
    cat \${tempfile}.log |\\
    awk '\$1=="Center"{print "COM", \$4, \$5, \$6}' |\\
    cat - \${tempfile}.xyz |\\
    awk '{printf "REMARK " \$0; getline; print "  ", \$2, \$3, \$4}' |\\
    cat >! \$file
    cat \${tempfile}.pdb |\\
    awk '/^REMARK/ && \$NF=="chain"{print "REMARK chain _";next}\\
          /^REMARK/{print}' >> \$file
    cat \${tempfile}.pdb |\\
    awk '! /REMARK/{print substr(\$0,1,66)}' >> \$file
    rm -f \${tempfile}.pdb \${tempfile}COM.pdb \${tempfile}.xyz \${tempfile}.log >& /dev/null

    # now "\$file" has its center of mass in its header
end

# display stats
echo -n "reference chains: "
foreach file ( \${tempfile}_right[0-9][0-9][0-9].pdb )
    awk '/^REMARK chain/{printf "%s ", \$3;exit}' \$file
end
echo "    ( "\`basename \$right_pdb\`" )"
echo -n "  subject chains: "
foreach file ( \${tempfile}_wrong[0-9][0-9][0-9].pdb )
    awk '/^REMARK chain/{printf "%s ", \$3;exit}' \$file
end
echo "    ( "\`basename \$wrong_pdb\`" )"
echo "CELL \$CELL"
echo "SG   \$SG"
echo ""







################################
# get tranformations from a space group

# get symmetry operations from space group
cat \${CLIBD}/symop.lib |\\
awk -v SG=\$SG '/[XYZ]/ && ! /[PCIFRH]/ && p==1 {print \$1 \$2 \$3 \$4}\\
\$5 ~ /^PG/{p=0} \$4 == SG{p=1}' |\\
cat >! \${tempfile}symops
set symops = \`cat \${tempfile}symops\`

# get origin-shift operators from a space group (at the bottom of this script)
set table = \$0
cat \$table |\\
awk '/TABLE OF ALLOWED ORIGIN SHIFTS/,/END OF THE TABLE/' |\\
awk -v SG=\$SG '/^[PCIFRH]/ && \$0 ~ SG" "{getline; while(NF>0){\\
    print \$1 "," \$2 "," \$3;getline};exit}' |\\
cat >! \${tempfile}origins
# format: dX,dY,dZ
set origins = \`cat \${tempfile}origins\`
#set origins = \`awk 'BEGIN{for(y=0;y<1;y+=0.01){print "0,"y",0";print "1/2,"y",0";print "0,"y",1/2";print "1/2,"y",1/2"}}'\`

if(\$?NO_ORIGINS || "\$origins" == "") then
    set origins = "0,0,0"
    set best_reindex_origin = "+X,+Y,+Z_0,0,0"
    set CCSG = 1
endif

rm -f \${tempfile}origins
rm -f \${tempfile}symops



# preemtive:
rm -f \${tempfile}all_reindex_origins.txt

# use deconvolution to find optimal shift
foreach reindexing ( \$reindexings )
    echo "applying \$reindexing to \$wrong_pdb"
    echo "symgen \$reindexing" |\\
      pdbset xyzin \${wrong_pdb} xyzout \${tempfile}reindexed.pdb > /dev/null

    echo "deconvoluting maps..."

    # do a coarser map grid sampling than normal, ~ 2A resolution
    set reso = 2

    # make a "mask" of all possible origin shifts
    echo -n "" >! \${tempfile}xyz.txt
    if("\$SG" == "P1") then
        echo "0 0 0" >> \${tempfile}xyz.txt
    else
        foreach origin ( \$origins )
            # break up the origin string
            set Xf = \`echo "\$origin" | awk -F "," '{print \$1}'\`
            set Yf = \`echo "\$origin" | awk -F "," '{print \$2}'\`
            set Zf = \`echo "\$origin" | awk -F "," '{print \$3}'\`
            set X = \`echo "\$Xf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`
            set Y = \`echo "\$Yf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`
            set Z = \`echo "\$Zf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`
        
            echo "\$X \$Y \$Z \$Xf \$Yf \$Zf" |\\
            awk '{print \$1+0,\$2+0,\$3+0}\\
             \$4=="x"{for(x=0;x<1;x+=0.001)print x,\$2+0,\$3+0}\\
             \$5=="y"{for(y=0;y<1;y+=0.001)print \$1+0,y,\$3+0}\\
             \$6=="z"{for(z=0;z<1;z+=0.001)print \$1+0,\$2+0,z}\\
             \$4=="x=y=z"{for(x=0;x<2;x+=0.001)print x,x,x}' |\\
            cat >> \${tempfile}xyz.txt
        end
    endif
    cat \${tempfile}xyz.txt |\\
    awk 'NF>2{++n; printf "%5d%10.5f%10.5f%10.5f%10.5f%5.2f%5d%10d%2s%3s%3s %1s\\n", \\
               n, \$1, \$2, \$3, 80, 1, "38", n, "H", "", "IUM", " "}' |\\
    cat >! \${tempfile}.frac
    coordconv XYZIN \${tempfile}.frac \\
          XYZOUT \${tempfile}originmask.pdb << EOF-conv >& /dev/null
CELL \$CELL
INPUT FRAC
OUTPUT PDB ORTH 1
END
EOF-conv
    sfall xyzin \${tempfile}originmask.pdb hklout \${tempfile}originmask.mtz << EOF-sfall > /dev/null
MODE sfcalc xyzin
CELL \$CELL
SYMM 1
RESO \$reso
EOF-sfall
    fft hklin \${tempfile}originmask.mtz mapout \${tempfile}originmask.map << EOF >! \${tempfile}.log
    labin F1=FC PHI=PHIC
    scale F1 1 80
    reso \$reso
    symm P1
EOF
    # turn it into a binary mask
    set scale = \`awk '/Maximum density/ && \$NF>0{print 1/\$NF}' \${tempfile}.log\`
    rm -f \${tempfile}.log
    if("\$scale" == "") set scale = 1
    echo "scale factor \$scale 0" |\\
      mapmask mapin \${tempfile}originmask.map mapout \${tempfile}new.map > /dev/null
    mv \${tempfile}new.map \${tempfile}originmask.map
    echo "scale factor 0 1" |\\
      mapmask mapin \${tempfile}originmask.map mapout \${tempfile}one.map > /dev/null

#    echo "mask cut 0" |\\
#    mapmask mapin \${tempfile}originmask.map mskout \${tempfile}originmask.msk > /dev/null
#    echo "maps mult" |\\
#    mapmask mapin \${tempfile}one.map mskin \${tempfile}originmask.msk \\
#     mapout \${tempfile}originmask.map > /dev/null
    if("\$SG" == "P1") then
        # P1 origin is good everywhere, so no mask
        cp \${tempfile}one.map \${tempfile}originmask.map > /dev/null
    endif
    rm -f \${tempfile}one.map \${tempfile}originmask.mtz
    rm -f \${tempfile}originmask.pdb \${tempfile}.frac \${tempfile}xyz.txt


    # make an all-carbon version of this model
    echo "CELL \$CELL" | pdbset xyzin \$right_pdb xyzout \${tempfile}.pdb >! \${tempfile}.log
    if(\$status) exit 9
    egrep "^CRYST1|^ATOM|^HETAT" \${tempfile}.pdb |\\
    awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
    cat >! \${tempfile}mapme.pdb
    sfall xyzin \${tempfile}mapme.pdb hklout \${tempfile}right.mtz << EOF-sfall >! \${tempfile}.log
    MODE sfcalc xyzin
    CELL \$CELL
    SYMM \$SG
    RESO \$reso
EOF-sfall
    if(\$status) exit 9
    echo "CELL \$CELL" | pdbset xyzin \${tempfile}reindexed.pdb xyzout \${tempfile}.pdb > /dev/null
    egrep "^CRYST1|^ATOM|^HETAT" \${tempfile}.pdb |\\
    awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
    cat >! \${tempfile}mapme.pdb
    sfall xyzin \${tempfile}mapme.pdb hklout \${tempfile}wrong.mtz << EOF-sfall >! \${tempfile}.log
    MODE sfcalc xyzin
    CELL \$CELL
    SYMM \$SG
    RESO \$reso
EOF-sfall
    if(\$status) exit 9

    rm -f \${tempfile}del.mtz
    sftools << EOF > /dev/null
read \${tempfile}right.mtz
read \${tempfile}wrong.mtz
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
write \${tempfile}del.mtz col Fdel PHIdel
y
stop
EOF
    fft hklin \${tempfile}del.mtz mapout \${tempfile}del.map << EOF >! \${tempfile}.log
    labin F1=Fdel PHI=PHIdel
    reso \$reso
    symm P1
EOF
    # make sure that we define "sigma" for the unmasked map
    echo "scale sigma 1 0" |\\
    mapmask mapin \${tempfile}del.map mapout \${tempfile}zscore.map  > /dev/null
    mapmask mapin1 \${tempfile}zscore.map mapin2 \${tempfile}originmask.map \\
     mapout \${tempfile}pickme.map << EOF >! \${tempfile}.log
mode mapin1 mapin2
maps mult
EOF
    peakmax mapin \${tempfile}pickme.map xyzfrc \${tempfile}peak.txt << EOF >! \${tempfile}.log
    output frac
#    threshold 3
    numpeaks 10
EOF
    set scale = \`awk '/peaks higher than the threshold/{print \$(NF-1)/\$9}' \${tempfile}.log\`
    echo "fractional origin shift   Z-score"
    cat \${tempfile}peak.txt |\\
    awk -v scale=\$scale '/ATOM|^HETAT/{printf "%7.4f %7.4f %7.4f   %s\\n", \$3,\$4,\$5,\$6/scale}' |\\
    awk 'NR==1{max=\$4} \$4>max/3 || \$4>3' |\\
    tee \${tempfile}likely_origins.txt

    # clean up a bit
    rm -f \${tempfile}right.mtz \${tempfile}wrong.mtz
    rm -f \${tempfile}del.mtz \${tempfile}del.map
    rm -f \${tempfile}originmask.map \${tempfile}pickme.map \${tempfile}zscore.map
    rm -f \${tempfile}peak.txt \${tempfile}.log

    # now integrate these shifts with the "official" origin list
    echo "ORIGINS \$origins" |\\
    cat - \${tempfile}likely_origins.txt |\\
    awk '/^ORIGINS/{\\
      for(i=2;i<=NF;++i){\\
        ++n;\\
        split(\$i,xyz,",");\\
        for(k in xyz){\\
          if(xyz[k]~/\\//){\\
            split(xyz[k],w,"/");\\
            xyz[k]=w[1]/w[2];\\
          }\\
        }\\
        x0[n]=xyz[1];y0[n]=xyz[2];z0[n]=xyz[3];\\
      }\\
    }\\
    NF==4{x=\$1;y=\$2;z=\$3;score=\$4;\\
      minfd=999;\\
      for(i in x0){\\
        xp=x0[i];yp=y0[i];zp=z0[i];\\
        if(xp=="x")xp=x;\\
        if(yp=="y")yp=y;\\
        if(zp=="z")zp=z;\\
        if(xp=="x=y=z")xp=yp=zp=x;\\
        origin[i]= xp","yp","zp;\\
        fdx=sqrt((x-xp)^2);if(fdx>0.9)fdx-=1;\\
        fdy=sqrt((y-yp)^2);if(fdy>0.9)fdy-=1;\\
        fdz=sqrt((z-zp)^2);if(fdz>0.9)fdz-=1;\\
        fd=sqrt(fdx^2+fdy^2+fdz^2);\\
        if(fd<minfd){minfd=fd;best=i};\\
      };\\
      print score,1-minfd,origin[best], x,y,z;\\
    }' |\\
    awk -v op="\$reindexing" '{print \$0,op}' |\\
    tee -a \${tempfile}all_reindex_origins.txt > /dev/null
    # format: score rel_height x,y,z  x y z  reindex 

    rm -f \${tempfile}likely_origins.txt

    set goodenough = \`awk '\$1>50' \${tempfile}all_reindex_origins.txt | wc -l\`
    if(("\$goodenough" == "1") && (\$?SPEEDUP)) then
        echo "thats good enough..."
        break
    endif
end

sort -gr \${tempfile}all_reindex_origins.txt >! \${tempfile}neworigins.txt
set reindex_origins = \`awk '{ro=\$7"_"\$3} ! seen[ro]{print ro} {++seen[ro]}' \${tempfile}neworigins.txt\`
set best_reindex = \`awk '{print \$7;exit}' \${tempfile}neworigins.txt\`
#rm -f \${tempfile}all_reindex_origins.txt \${tempfile}neworigins.txt



########################################
# now run through all origins, and chain pairings
echo "chain       origin "
echo "r    s      reindex        dX  dY  dZ  symop                \$SCORE"

again:
echo -n "" >! \${tempfile}scores
foreach right_chain ( \${tempfile}_right[0-9][0-9][0-9].pdb )

# get fractional coordinate limits that cover the "right" chain
coordconv xyzin \$right_chain xyzout \${tempfile}.xyz << EOF >> /dev/null
INPUT PDB
OUTPUT FRAC
END
EOF
cat \${tempfile}.xyz |\\
awk -v del=0.5 'BEGIN{xmin=ymin=zmin=99999999} \\
\$2<xmin{xmin=\$2} \$2>xmax{xmax=\$2}\\
\$3<ymin{ymin=\$3} \$3>ymax{ymax=\$3}\\
\$4<zmin{zmin=\$4} \$4>zmax{zmax=\$4}\\
END{print xmin-del, xmax+del, ymin-del, ymax+del, zmin-del, zmax+del}' |\\
cat >! \${tempfile}.xyzlim
set xyzlim = \`cat \${tempfile}.xyzlim\`
rm -f \${tempfile}.xyzlim
rm -f \${tempfile}.xyz


if(\$?CORRELATE) then
    # we want to do comparison with maps, not atoms

    if(! \$?coarseGRID) then
        # do a coarser map grid sampling than normal (~ 3A resolution)
        set reso = 3
        set fakeCELL = \`echo \$CELL \$reso | awk '{print \$1/\$NF/2,\$2/\$NF/2,\$3/\$NF/2,\$4,\$5,\$6}'\`
        set BADD     = \`echo \$reso | awk '{printf "%d", 79*(\$1/3)^2}'\`
        echo "CELL \$fakeCELL" | pdbset xyzin \$right_chain xyzout \${tempfile}.pdb > /dev/null

        # make an all-carbon version of this model
        cat \${tempfile}.pdb |\\
        awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
        cat >! \${tempfile}mapme.pdb
        sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}right.map << EOF-sfall > /dev/null
        MODE ATMMAP
        CELL \$fakeCELL
        SYMM \$CCSG
EOF-sfall
        # recover coarse grid spacing 
        set coarseGRID = \`echo "GO" | mapdump MAPIN \${tempfile}right.map | awk '/Grid sampling/{print \$(NF-2), \$(NF-1), \$NF; exit}'\`
        rm -f \${tempfile}.pdb >& /dev/null
    endif
    set GRID = "\$coarseGRID"

    # re-calculate the "right" map with reduced grid
    cat \$right_chain |\\
    awk '/^ATOM|^HETAT/{++n;\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
    awk '/^ATOM|^HETAT/{++n;}\\
     n==1{++n;\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\\
         {print}' |\\
    cat >! \${tempfile}mapme.pdb
    sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}right.map << EOF-sfall > /dev/null
    MODE ATMMAP
    CELL \$CELL
    SYMM \$CCSG
    GRID \$GRID
    BADD \$BADD
EOF-sfall
    rm -f \${tempfile}mapme.pdb >& /dev/null
    # determine grid spacing to use for "wrong" map
    set GRID = \`echo "GO" | mapdump MAPIN \${tempfile}right.map | awk '/Grid sampling/{print \$(NF-2), \$(NF-1), \$NF; exit}'\`
#    set xyzlim = "0 1 0 1 0 1"
endif

# now loop over the wrong chains
foreach wrong_chain ( \${tempfile}_wrong[0-9][0-9][0-9].pdb )


# try every possible origin choice
foreach reindex_origin ( \$reindex_origins )

# break up the origin string
set reindex = \`echo \$reindex_origin | awk -F "_" '{print \$1}'\`
set origin = \`echo \$reindex_origin | awk -F "_" '{print \$2}'\`
set Xf = \`echo "\$origin" | awk -F "," '{print \$1}'\`
set Yf = \`echo "\$origin" | awk -F "," '{print \$2}'\`
set Zf = \`echo "\$origin" | awk -F "," '{print \$3}'\`
set X = \`echo "\$Xf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`
set Y = \`echo "\$Yf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`
set Z = \`echo "\$Zf" | awk -F "/" 'NF==2{\$1/=\$2} {printf "%.10f",\$1}'\`

# apply this re-indexing to the "wrong" pdb
pdbset xyzin \${wrong_chain} xyzout \${tempfile}wrong_reindexed.pdb << EOF > /dev/null
symgen \$reindex
EOF

# get the center of mass in orthogonal coordinates
pdbset xyzin \${tempfile}wrong_reindexed.pdb xyzout \${tempfile}.pdb << EOF >! \${tempfile}.log
COM
EOF
set wrong_COM = \`awk '\$1=="Center"{print \$4,\$5,\$6}' \${tempfile}.log\`

# make a pdb file of an atom at the COM of this file
awk '/^CRYST/ || /^SCALE/' \${tempfile}.pdb >! \${tempfile}COM.pdb
echo \$wrong_COM |\\
awk '{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\\n", \$1, \$2, \$3;exit}' |\\
cat >> \${tempfile}COM.pdb
echo "END" >> \${tempfile}COM.pdb

# get fractional version of the COM too
coordconv XYZIN \${tempfile}COM.pdb XYZOUT \${tempfile}.xyz << EOF >& /dev/null
CELL \$CELL
INPUT PDB
OUTPUT FRAC
EOF
set wrong_COM_frac = \`awk '{print \$2,\$3,\$4}' \${tempfile}.xyz\`


# retrieve info from file headers
set right_COM = \`awk '/^REMARK COM/{print \$3, \$4, \$5;exit}' \$right_chain\`
#set wrong_COM = \`awk '/^REMARK COM/{print \$3, \$4, \$5;exit}' \$wrong_chain\`
set shift_COM = \`echo "\$wrong_COM \$right_COM" | awk '{print \$4-\$1,\$5-\$2,\$6-\$3}'\`
set right_COM_frac = \`awk '/^REMARK COM/{print \$6, \$7, \$8;exit}' \$right_chain\`
#set wrong_COM_frac = \`awk '/^REMARK COM/{print \$6, \$7, \$8;exit}' \$wrong_chain\`
set shift_COM_frac = \`echo "\$wrong_COM_frac \$right_COM_frac" | awk '{print \$4-\$1,\$5-\$2,\$6-\$3}'\`
set right_chain_ID = \`awk '/^REMARK chain/{print \$3}' \$right_chain\`
set wrong_chain_ID = \`awk '/^REMARK chain/{print \$3}' \$wrong_chain\`
if("\$right_chain_ID" == "") set right_chain_ID = "_"
if("\$wrong_chain_ID" == "") set wrong_chain_ID = "_"

# add translation along any polar axes
#set X = \`echo "\$opt_frac_shift \$Xf \$X" | awk '/x/{\$NF=\$1} {print \$NF}'\`
#set Y = \`echo "\$opt_frac_shift \$Yf \$Y" | awk '/y/{\$NF=\$2} {print \$NF}'\`
#set Z = \`echo "\$opt_frac_shift \$Zf \$Z" | awk '/z/{\$NF=\$3} {print \$NF}'\`

# calculate how to "level" the centers of the atom constellations
# NOTE: this is not a good idea if the chains are not rigid-body related
#set polar_slip = \`echo "\$opt_frac_shift \$Xf \$Yf \$Zf" | awk '! /x/{\$1=0} ! /y/{\$2=0} ! /z/{\$3=0} {print \$1,\$2,\$3}'\`

# now see which (if any) symops/cell translations will bring the 
# shifted COM of "wrong_chain" anywhere near the COM of the "right" chain

# apply current origin shift to "wrong" COM
set new_COM = \`echo "\$wrong_COM_frac \$X \$Y \$Z" | awk '{print \$1+\$4, \$2+\$5, \$3+\$6}'\`
gensym << EOF >! \${tempfile}.log
CELL \$CELL
SYMM \$SG
atom X \$new_COM
XYZLIM \$xyzlim
EOF

# only use symops that showed up in the gensym result
cat \${tempfile}.log |\\
awk '/List of sites/,/atoms generated from/' |\\
awk 'NF>10{print \$NF, \$2, \$3, \$4, \$5, \$6, \$7}' |\\
sort -u -k1,1 >! \${tempfile}symop_center
# format: symop xf yf zf X Y Z

if(\$?CORRELATE && ! \$?best_reindex_origin) then
    echo "1 \$right_COM_frac \$new_COM" >! \${tempfile}symop_center
endif

# loop over the symmetry operations
foreach line ( \`awk '{print NR}' \${tempfile}symop_center\` )

    # skip all other origins if a good one has been found
    if(("\$good_origin" != "")&&("\$good_origin" != "\$X \$Y \$Z") && "\$CCSG" != "1") then
        # echo "skipping: \$good_origin   \$X \$Y \$Z"
        continue
    endif
    # skip rest of "right_chain" if we have already found a match
    if("\$good_right" == "\$right_chain") then
        # echo "skipping: \$good_right == \$right_chain"
        continue
    endif
    if("\$good_wrong" == "\$wrong_chain") then
        # echo "skipping: \$good_wrong == \$wrong_chain"
        continue
    endif

    # retrieve symmetry operation
    set symop = \`awk -v line=\$line 'NR==line{print \$1}' \${tempfile}symop_center\`
    set symop = \$symops[\$symop]
    
    # progress meter
    echo "\$right_chain_ID \$wrong_chain_ID  \$reindex  \$Xf \$Yf \$Zf \$symop" |\\
     awk '{printf "%s vs %s by %-15s @ %3s %3s %3s %-20s", \$1, \$2, \$3, \$4,\$5,\$6, \$7}'
    
    # move the reindexed wrong PDB to the new position
    pdbset xyzin \${tempfile}wrong_reindexed.pdb \\
          xyzout \${tempfile}symmed.pdb << EOF >! \${tempfile}.log
    CELL \$wrong_cell
    SYMGEN \$symop
    SHIFT FRAC \$X \$Y \$Z
    COM
EOF
    set new_COM = \`awk '\$1=="Center"{print \$4, \$5, \$6}' \${tempfile}.log\`
    # get fractional version too
    awk '/^CRYST/ || /^SCALE/' \${tempfile}wrong_reindexed.pdb >! \${tempfile}COM.pdb
    cat \${tempfile}.log |\\
    awk '\$1=="Center"{printf "ATOM      1  CA  GLY     1    %8.3f%8.3f%8.3f  1.00 80.00\\n", \$4, \$5, \$6;exit}' |\\
    cat >> \${tempfile}COM.pdb
    echo "END" >> \${tempfile}COM.pdb
    coordconv XYZIN \${tempfile}COM.pdb XYZOUT \${tempfile}.xyz << EOF > /dev/null
    CELL \$wrong_cell
    INPUT PDB
    OUTPUT FRAC
EOF
    set new_COM_frac = \`awk '{print "  ", \$2, \$3, \$4}' \${tempfile}.xyz\`
    
    # now find the nearest cell translation to bring these COMs together
    set cell_shift = \`echo "\$right_COM_frac \$new_COM_frac" | awk '{printf "%.0f %.0f %.0f", \$1-\$4+0, \$2-\$5+0, \$3-\$6+0}'\`
    
    # move the "symmed" "wrong" chain into the right cell
    pdbset xyzin \${tempfile}symmed.pdb \\
          xyzout \${tempfile}moved.pdb << EOF >! \${tempfile}.log
    CELL \$wrong_cell
    SHIFT FRAC \$cell_shift
EOF
 
    # calclulate RMS distance between chains
    set score = \`awk -f \${tempfile}rmsd.awk \$right_chain \${tempfile}moved.pdb | awk '/RMSD\\(all/{print -\$2}'\`

    if(\$?CORRELATE) then
        # make sure all atoms are residues are readabel
        cat \${tempfile}moved.pdb |\\
        awk '/^ATOM|^HETAT/{\$0="ATOM      1  CA  ALA     1    "substr(\$0,31,25)" 1.00 80.00"} {print}' |\\
            awk '/^ATOM|^HETAT/{++n;}\\
         n==1{++n;\\
          print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\\
          print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\\
             {print}' |\\
        cat >! \${tempfile}mapme.pdb
        # create a map of these atoms
        sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}moved.map << EOF-sfall > /dev/null
        MODE ATMMAP
        CELL \$CELL
        SYMM \$CCSG
        GRID \$GRID
        BADD \$BADD
EOF-sfall
        rm -f \${tempfile}mapme.pdb >& /dev/null
        
        # compare to the "right" map
        echo "correlate section" |\\
        overlapmap mapin1 \${tempfile}right.map \\
               mapin2 \${tempfile}moved.map mapout /dev/null |\\
        awk '/Total correlation/{print \$NF}' >! \${tempfile}correlation
        set score = \`awk '{print \$1}' \${tempfile}correlation\`
        
        # display CC
        echo "\$score" | awk '\$1!~/[^0-9]/{print "  -";exit} {printf "%5.2f\\n", \$1}'
    else
        # display rmsd
        echo "\$score" | awk '\$1~/[0-9]/{printf "%5.1f\\n", -\$1; exit} {print "no identities"}'
    endif
    rm -f \${tempfile}symmed.pdb \${tempfile}moved.pdb
    if("\$score" == "") continue
    
    # see if fit is "good enough"
    set goodenough = \`echo \$score | awk '{print (\$1*\$1<1)}'\`
    if(\$?CORRELATE) set goodenough = \`echo \$score | awk '{print (\$1>0.8)}'\`
    if((\$goodenough == 1)&&(\$?SPEEDUP)) then
        echo "thats good enough..."
        set good_origin = "\$X \$Y \$Z"
        set good_right = "\$right_chain"
        set good_wrong = "\$wrong_chain"
    endif
    
    # make a running log
    echo "\$right_chain \$wrong_chain \$reindex \$Xf \$Yf \$Zf \$X \$Y \$Z \$cell_shift \$symop \$score " >> \${tempfile}scores

end
end
end
end

echo ""

# select out the "best" matches (on the same origin)
sort -nr -k14 \${tempfile}scores |\\
awk '{ros=\$3" "\$4" "\$5" "\$6;chain_ros=\$2" "ros} ! seen[chain_ros]{++seen[chain_ros];\\
        score[ros]+=\$NF;++count[ros]}\\
     END{for(ros in score) print score[ros], count[ros], ros}' |\\
sort -nr >! \${tempfile}best_oriscores 
# format: score n X,Y,Z dX dY dZ
set best_reindex_origin = \`head -1 \${tempfile}best_oriscores | awk '{print \$3, \$4, \$5, \$6}'\`

if(\$?CORRELATE && "\$CCSG" != "1" && ! \$?NO_ORIGINS) then
    # need to go back and do this again for symops
    set reindex_origins = \`echo \$best_reindex_origin | awk '{print \$1"_"\$2","\$3","\$4}'\`
    set CCSG = 1
    set good_right = ""
    set good_wrong = ""
    echo "now checking symops"
    goto again
endif

# extract the best operators consistent with this origin
sort -nr -k14 \${tempfile}scores |\\
awk -v best_reindex_origin="\$best_reindex_origin" '\$3" "\$4" "\$5" "\$6==best_reindex_origin{print}' |\\
awk '! seen[\$2]{print; seen[\$2]=1}' |\\
sort >! \${tempfile}best_scores 
#format r_chain w_chain X,Y,Z 1/2 1/2 1/2 0.000 0.000 0.000  1 0 1  X,Y,Z  rmsd


echo -n "" >! \${tempfile}out.pdb
foreach line ( \`awk '{print NR}' \${tempfile}best_scores\` )
    # retrieve chain info
    set right_chain = \`awk -v line=\$line 'NR==line{print \$1}' \${tempfile}best_scores\`
    set wrong_chain = \`awk -v line=\$line 'NR==line{print \$2}' \${tempfile}best_scores\`

    set right_chain_ID = \`awk '/^REMARK chain/{print \$3}' \$right_chain\`
    set wrong_chain_ID = \`awk '/^REMARK chain/{print \$3}' \$wrong_chain\`
    if("\$right_chain_ID" == "") set right_chain_ID = "_"
    if("\$wrong_chain_ID" == "") set wrong_chain_ID = "_"
    

    # retrieve transofmation info
    set reindex = \`awk -v chain=\$wrong_chain '\$2==chain{print \$3}' \${tempfile}best_scores\`
    set X     = \`awk -v chain=\$wrong_chain '\$2==chain{printf "%.2f", \$7}' \${tempfile}best_scores\`
    set Y     = \`awk -v chain=\$wrong_chain '\$2==chain{printf "%.2f", \$8}' \${tempfile}best_scores\`
    set Z     = \`awk -v chain=\$wrong_chain '\$2==chain{printf "%.2f", \$9}' \${tempfile}best_scores\`
    set x     = \`awk -v chain=\$wrong_chain '\$2==chain{print \$7+\$10}' \${tempfile}best_scores\`
    set y     = \`awk -v chain=\$wrong_chain '\$2==chain{print \$8+\$11}' \${tempfile}best_scores\`
    set z     = \`awk -v chain=\$wrong_chain '\$2==chain{print \$9+\$12}' \${tempfile}best_scores\`
    set symop = \`awk -v chain=\$wrong_chain '\$2==chain{print \$13}' \${tempfile}best_scores\`

    # retrieve score
    set score  = \`awk -v chain=\$wrong_chain '\$2==chain{print -\$14}' \${tempfile}best_scores\`
    if(\$?CORRELATE) set score = \`echo \$score | awk '{print -\$1}'\`
    
    if("\$score" == "") then
        echo "unable to match \$wrong_pdb chain \$wrong_chain_ID to any chain in \$right_pdb"
        continue
    endif
        
    # apply it
    echo "\$reindex \$x \$y \$z \$symop" |\\
    awk '{printf "applying %-15s %5.2f %5.2f %5.2f %-15s\\n", \$1, \$2, \$3, \$4, \$5}'
    echo "to \$wrong_pdb chain \$wrong_chain_ID --> \$outfile chain \$right_chain_ID  (\$SCORE = \$score)"

    if("\$right_chain_ID" == "_") set right_chain_ID = " "
    pdbset xyzin \$wrong_chain xyzout \${tempfile}reindexed.pdb << EOF >> /dev/null
    SYMGEN \$reindex
EOF
    pdbset xyzin \${tempfile}reindexed.pdb xyzout \${tempfile}moved.pdb << EOF >> /dev/null
    CELL \$wrong_cell
    CHAIN "\$right_chain_ID"
    SYMGEN \$symop
    SHIFT FRAC \$x \$y \$z
EOF
    # update the output file
    egrep "^ATOM|^HETAT" \${tempfile}moved.pdb >> \${tempfile}out.pdb
    
    # warn about different origins
    if(! \$?last_reindex_origin) set last_reindex_origin
    if(("\$last_reindex_origin" != "\$reindex \$X \$Y \$Z")&&("\$last_reindex_origin" != "")) then
        echo "WARNING: \$wrong_pdb chain \$wrong_chain_ID is on a different origin! "
    endif
    set last_reindex_origin = "\$reindex \$X \$Y \$Z"
end
echo "END" >> \${tempfile}out.pdb

pdbset xyzin \${tempfile}out.pdb xyzout \$outfile << EOF > /dev/null
CELL \$wrong_cell
SPACE \$SG
EOF

if(\$?BYFILE) then
    # no point in messing around, just apply the shift to the original file
    pdbset xyzin \$wrong_pdb xyzout \${tempfile}out.pdb << EOF > /dev/null
CELL \$wrong_cell
symgen \$reindex
EOF
    pdbset xyzin \${tempfile}out.pdb xyzout \$outfile << EOF > /dev/null
CELL \$wrong_cell
SPACE \$SG
SYMGEN \$symop
SHIFT FRAC \$x \$y \$z
EOF

# now recover the operation in other conventions
lsqkab XYZIN1 \$outfile XYZIN2 \$wrong_pdb  << EOF >! \${tempfile}lsq.log
FIT ATOM 1 to 99999
MATCH    1 to 99999
END
EOF
set lsq_polar = \`awk '/OMEGA PHI CHI/ && \$1~/^SPHERICAL/{print \$(NF-2), \$(NF-1), \$NF}' \${tempfile}lsq.log\`
set lsq_trans = \`awk '/TRANSLATION VECTOR IN AS/{print \$(NF-2), \$(NF-1), \$NF}' \${tempfile}lsq.log\`
cat << EOF
general transformation for pdbset, maprot, ncsmask, dm, etc.:
ROTATE POLAR \$lsq_polar
TRANSLATE \$lsq_trans
EOF
endif

if(\$?SPEEDUP) goto exit

# report final, overall score
if(\$?CORRELATE) then
    # calculate the overall map correlation
    set logfile = /dev/null

    # create a simplified file for SFALL
    cat \$right_pdb |\\
    awk '/^CRYST|^END|^REM/{print}\\
            /^ATOM|^HETAT/{++n;printf("ATOM      1  CA  ALA %5d    %25s 1.00 80.00\\n",n,substr(\$0,31,25))}' |\\
    awk '/^ATOM|^HETAT/{++n;}\\
     0 && n==1{++n;\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\\
         {print}' |\\
    cat >! \${tempfile}mapme.pdb
    sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}right.map << EOF-sfall > \$logfile
    MODE ATMMAP
    CELL \$CELL
    SYMM \$SG
EOF-sfall
    # also record the "right" atom codes
    cat \$right_pdb |\\
    awk '/^ATOM|^HETAT/{++n;printf "code %04d %s\\n",n,\$0}' |\\
    cat >! \${tempfile}atomcodes.txt
    # recover grid spacing 
    set GRID = \`echo "GO" | mapdump MAPIN \${tempfile}right.map | awk '/Grid sampling/{print \$(NF-2), \$(NF-1), \$NF; exit}'\`

    # calculate the "label" map so we can assign atom-specific correlations
    sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}label.map << EOF-sfall >> \$logfile
    MODE ATMMAP RESMOD
    CELL \$CELL
    SYMM \$SG
    GRID \$GRID
EOF-sfall

    # calculate map for "moved" model
    cat \$outfile |\\
    awk '/^CRYST|^END|^REM/{print}\\
            /^ATOM|^HETATM/{++n;printf("ATOM      1  CA  ALA %5d    %25s 1.00 80.00\\n",n,substr(\$0,31,25))}' |\\
    awk '/^ATOM|^HETAT/{++n;}\\
     0 && n==1{++n;\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000  0.01 99.00";\\
      print "ATOM      0  CA  SHT     0       0.000   0.000   0.000 -0.01 99.00";}\\
         {print}' |\\
    cat >! \${tempfile}mapme.pdb
    # also record the "wrong" atom codes
    cat \$outfile |\\
    awk '/^ATOM|^HETAT/{++n;printf "code %04d %s\\n",n,\$0}' |\\
    cat >! \${tempfile}wrong_atomcodes.txt
    sfall xyzin \${tempfile}mapme.pdb mapout \${tempfile}moved.map << EOF-sfall > \$logfile
    MODE ATMMAP
    CELL \$CELL
    SYMM \$SG
    GRID \$GRID
EOF-sfall
# why do they make me do these things... sfall map can sometimes be too small
#dd if=/dev/zero bs=1024 count=10000 >>& \${tempfile}moved.map
#    sfall xyzin \${tempfile}mapme.pdb HKLOUT \${tempfile}moved.mtz << EOF-sfall > \$logfile
#    MODE SFCALC XYZIN
#    CELL \$CELL
#    SYMM \$CCSG
#EOF-sfall
#    fft hklin \${tempfile}moved.mtz mapout \${tempfile}moved.map << EOF-fft > \$logfile
#    LABIN F1=FC PHI=PHIC
#    GRID \$GRID
#EOF-fft


    rm -f \${tempfile}mapme.pdb >& /dev/null
    # compare to the "right" map
    echo "correlate residue" |\\
    overlapmap mapin1 \${tempfile}right.map \\
           mapin2 \${tempfile}moved.map \\
           mapin3 \${tempfile}label.map mapout /dev/null |\\
    awk '/Main-corr-coef/,/Total/{print}' | tee \${tempfile}correlation |\\
    awk '\$1+0>0{print}' |\\
    awk '\$2+0>0{print \$1,\$2;next} \$3+0>0{print \$1,\$3}' |\\
    awk '{print "ATOMCC",\$0}' >! \${tempfile}atom_corr.txt
    # format: ATOMCC rightatomnum CC

    # encode the atomic CCs as an occupancy in a new PDB file
    cat \${tempfile}atomcodes.txt \${tempfile}atom_corr.txt \$right_pdb |\\
    awk '/^code /{name[\$2+0]=substr(\$0,11,26);next} /^ATOMCC /{occ[\$2]=\$3;next} ! /^ATOM/ {print}\\
      /^ATOM|^HETATM/{++n;\\
    if(name[n]=="" || occ[n]+0<0.01){name[n]=substr(\$0,1,26);occ[x]=substr(\$0,55,6)+0};\\
    printf "%26s%28s%6.2f%s\\n", name[n],substr(\$0,27,28),occ[n],substr(\$0,61)}' |\\
    cat >! newlabel_\$outfile

    set score = \`awk '/Total/{print \$NF;print \$(NF-1)}' \${tempfile}correlation | sort -nr | head -1\`

    rm -f \${tempfile}right.map >& /dev/null
    rm -f \${tempfile}moved.map >& /dev/null
    rm -f \${tempfile}label.map >& /dev/null
    rm -f \${tempfile}correlation >& /dev/null    
    rm -f \${tempfile}atom_corr.txt >& /dev/null    
    rm -f \${tempfile}atomcodes.txt >& /dev/null    
    rm -f \${tempfile}wrong_atomcodes.txt >& /dev/null
else
    # the overall RMSD is the score
    set score = \`awk -f \${tempfile}rmsd.awk \$right_pdb \$outfile | awk '/RMSD\\(all/{print \$2} /no atom pairs/{print "no identities"}'\`
endif
echo "Overall \${SCORE}: \$score"



#clean up
exit:
if(\$?debug) exit
rm -f \${tempfile}rmsd.awk >& /dev/null
rm -f \${tempfile}out.pdb >& /dev/null
rm -f \${tempfile}moved.pdb >& /dev/null
rm -f \${tempfile}best_scores >& /dev/null
rm -f \${tempfile}scores >& /dev/null
rm -f \${tempfile}.log >& /dev/null
rm -f \${tempfile}.xyz >& /dev/null
rm -f \${tempfile}COM.pdb >& /dev/null
rm -f \${tempfile}symop_center >& /dev/null
rm -f \${tempfile}_wrong???.pdb >& /dev/null
rm -f \${tempfile}_right???.pdb >& /dev/null
rm -f \${tempfile}* >& /dev/null

exit

# Notes on hand flipping

cat rh.pdb |\\
awk -v SG=\$SG '/^CRYST1/{a=\$2;c=\$4;\\
     if(SG=="F4132"){dorx=dory=dorz=a*0.25};\\
     if(SG=="I41"){dorx=a*0.5};\\
     if(SG=="I4122"){dorx=a*0.5;dorz=c*0.25};\\
  } ! /^ATOM|^HETAT/{print} /^ATOM|^HETAT/{;\\
  x=substr(\$0,31,8)+0;y=substr(\$0,39,8)+0;z=substr(\$0,47,8)+0;\\
  printf("%s%8.3f%8.3f%8.3f%s\\n", substr(\$0,1,30),dorx-x,dory-y,dorz-z,substr(\$0,55))}' |\\
cat >! lh.pdb




################################################################################

Setup:
set good_origin = ""
set good_right = ""
set good_wrong = ""

foreach arg ( \$* )

    if(("\$arg" =~ *.pdb)||("\$arg" =~ *.brk)) then
    # warn about probable mispellings
    if(! -e "\$arg") then
        echo "WARNING: \$arg does not exist"
        continue
    endif
    # make sure its really a pdb file
    egrep -l "^ATOM |HETATM " "\$arg" >& /dev/null
    if(\$status) then
        echo "WARNING: \$arg contains no atoms! "
        continue
    endif
    
    if(-e "\$right_pdb") then
        set wrong_pdb = "\$arg"
    else
        set right_pdb = "\$arg"
    endif
    continue
    endif

    # space group
    if("\$arg" =~ [PpCcIiFfRrHh][1-6]*) then
    set temp = \`echo \$arg | awk '{print toupper(\$1)}'\`
    if(\$?CLIBD) then
        set temp = \`awk -v SG=\$temp '\$4 == SG {print \$4}' \$CLIBD/symop.lib | head -1\`
    endif
    if("\$temp" != "") then
        # add this SG to the space group list
        set SG = "\$temp"
    endif
    endif

    # flags and options
    if("\$arg" == "correlate") then
    # correlation instead of RMSD
    set CORRELATE
    endif
    if("\$arg" == "otherhand") then
        set OTHERHAND
    endif
    if("\$arg" == "altindex") then
        set ALTINDEX
    endif
    if("\$arg" == "nochains" || "\$arg" == "byfile") then
        # whole file is one chain
        set BYFILE
    endif
    if("\$arg" == "noorigins") then
        # just symmetry search
        set NO_ORIGINS
    endif
    if("\$arg" == "fast") then
        # skip several steps
        set SPEEDUP
    endif
end

if((! -e "\$right_pdb")||(! -e "\$wrong_pdb")) then
    goto Help
endif


# deploy rmsd awk program
cat << EOF-script >! \${tempfile}rmsd.awk
#! \$awk -f
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
    ID = substr(\\\$0,12,15)
    ++count[ID]

    if(count[ID] == 1)
    {
    # not initialized yet
    X[ID] = substr(\\\$0, 31, 8)+0
    Y[ID] = substr(\\\$0, 39, 8)+0
    Z[ID] = substr(\\\$0, 47, 8)+0
    B[ID] = substr(\\\$0, 61, 6)+0
    }
    
    if(count[ID] == 2)
    {
    ++pairs
    
    # seen this before, subtract values
    dX     = X[ID] - substr(\\\$0, 31, 8)
    dY     = Y[ID] - substr(\\\$0, 39, 8)
    dZ     = Z[ID] - substr(\\\$0, 47, 8)
    dB[ID] = B[ID] - substr(\\\$0, 61, 6)
    
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
        printf("%s moved %8.4f (XYZ) %6.2f (B)\\\\n", ID, dXYZ[ID], dB[ID])
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

    print "MAXD(all)= " maxXYZ "\\\\tfor " maxXYZ_ID
    print "MAXD(Bfac)= " maxdB "\\\\tfor " maxdB_ID
    
    # final check for orphan atoms 
    for(id in count)
    {
        if(count[id]<2) print "WARNING: " id " only found once"
    }
    }
    else
    {
    printf "%10.8f %10.8f %10.5f %10.8f %8.2f \\\\n", rms_atom_XYZ, rmsXYZ, rmsB, maxXYZ, maxdB
    }
}
EOF-script
chmod a+x \${tempfile}rmsd.awk


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
EOF-origins.com
chmod a+x origins.com




















cat << EOF-script >! floatgen.c
/*
**
**      floatgen.c  v0.1                                             James Holton 5-19-11
**
**      converts text number on stdin to a 4-byte float on stdout
**
**      compile this file with:
**              gcc -o floatgen floatgen.c -lm -static
**
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argv, char** argc)
{
    float fourbytes;
    char text[1024];
    char *token;
    const char delimiters[] = " \\t,;:!";
    const char numberstuf[] = "0123456789-+.EGeg";

    while ( fgets ( text, sizeof text, stdin ) != NULL ) {

        token = text;
        token += strspn(token,delimiters);
        if(strcmp(token,"\\n")==0) {
            //printf("blank\\n");
            continue;
        }

        fourbytes=atof(token);
 
	fwrite(&fourbytes,sizeof(float),1,stdout);
    }
}
EOF-script
gcc -o floatgen floatgen.c -lm






















# make sure this happens! 
set path = ( $path . )
rehash

goto return_from_deploy


