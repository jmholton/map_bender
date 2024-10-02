#! /bin/tcsh -f
#
#    re-compile mapman so that it works on 64-bit machines
#
#
# retrieve the Uppsala Software Factory distribution
wget http://xray.bmc.uu.se/usf/usf_distribution_kit.tar.gz
# in case above is down, this might work instead
#wget http://bl831.als.lbl.gov/~jamesh/pickup/usf_distribution_kit.tar.gz

# check it
md5sum usf_distribution_kit.tar.gz
# 0d01d5e0cb5072822870004bb29e331e  usf_distribution_kit.tar.gz

# unpack it
tar xzvf usf_distribution_kit.tar.gz
cd usf_export

# fix the bug
patch mapman/mapman.f << EOF
28,29c28,29
<       integer iaptr,ibptr
<       integer fmalloc
---
>       integer*8 iaptr,ibptr
>       integer*8 fmalloc
EOF

# fix other bug
patch gklib/fmalloc.c << EOF
24c24
< typedef int address_type;
---
> typedef long address_type;
EOF

# run the re-compilation script
./make_all.csh mapman -64 -static

# if it doesnt work, try removing static
#./make_all.csh mapman -64


cd mapman

# test it
wget https://raw.githubusercontent.com/fraser-lab/holton_scripts/master/map_bender/mapman_regression_test.csh
chmod a+x mapman_regression_test.csh
./mapman_regression_test.csh ../bin/mapman

exit

# if things go wrong, try:

yum -y install tcsh
yum -y install patch
yum -y install wget
yum -y install gcc-gfortran
yum -y install libgfortran-static
yum -y install gcc*
yum -y install libgfortran4
yum -y install glibc-static
yum -y install lib\*static
yum -y install lapack-*


# rebuilding

cd ../gklib
./make_fresh_gklib.csh Linux 64
cd ../mapman

rm mapman *.o
make -f Makefile_linux
cat mapman.in | ./mapman


