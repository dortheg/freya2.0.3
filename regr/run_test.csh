#!/bin/csh
#
# run_test [1|2|3|4|5|6]
#
# default is to run all 6
#
# Jerome Verbeke
# mods by Doug Wright

set CC  = gcc
set CXX = g++
set FC  = gfortran

#....choose linker based on availability of fortran compiler
which $FC > /dev/null && set LD=$FC || set LD=$CXX

setenv LD_LIBRARY_PATH ../lib
setenv DYLD_LIBRARY_PATH ${LD_LIBRARY_PATH} # MacOSX
setenv INCLUDE '-I../include -L../lib'
setenv LIBS '-lFission -lstdc++ -lm'
# setenv EXTRAFCFLAGS '-nofor_main'
setenv EXTRAFCFLAGS ''

#....process command line arguments (select which test to run)
if ($#argv == 0 ) then
    set test1 test2 test3 test4 test5 test6
else
    set test$1
endif

#....compile and run the tests
if ($?test1) then
   echo Compiling and linking testSpNuDist
   $CC testSpNuDist.c -c -g $INCLUDE
   $LD testSpNuDist.o -g -o testSpNuDist $EXTRAFCFLAGS $INCLUDE $LIBS

   echo "running testSpNuDist (takes 13 seconds on aztec)"
   if (-e testSpNuDist.res) rm -f testSpNuDist.res
   time ./testSpNuDist
   diff -s testSpNuDist.res testSpNuDist.goldref
   diff -s testSpspec.res testSpspec.goldref
   echo
endif

if ($?test2) then
   echo Compiling and linking testNuDist
   $CC  testNuDist.c -g -c $INCLUDE
   $LD  testNuDist.o -g -o testNuDist $EXTRAFCFLAGS $INCLUDE $LIBS

   echo "running testNuDist (takes 38 seconds on aztec)"
   if (-e testNuDist.res) rm -f testNuDist.res
   if (-e testspec.res) rm -f testspec.res
   time ./testNuDist
   diff -s testNuDist.res testNuDist.goldref
   diff -s testspec.res testspec.goldref
   echo
endif

if ($?test3) then
   echo Compiling and linking testEngCons
   $CC  testEngCons.c -g -c $INCLUDE
   $LD  testEngCons.o -g -o testEngCons $EXTRAFCFLAGS $INCLUDE $LIBS

   echo "running testEngCons (takes 29 seconds on aztec)"
   if (-e testEngConsNudist.res) rm -f testEngConsNudist.res
   if (-e testEngConsSpec.res) rm -f testEngConsSpec.res
   if (-e testEngConsSpecTotal.res) rm -f testEngConsSpecTotal.res
   if (-e testEngConsTotal.res) rm -f testEngConsTotal.res
   time ./testEngCons
   diff -s testEngConsNudist.res testEngConsNudist.goldref
   diff -s testEngConsSpec.res testEngConsSpec.goldref
   diff -s testEngConsSpecTotal.res testEngConsSpecTotal.goldref
   diff -s testEngConsTotal.res testEngConsTotal.goldref
   echo
endif

if ($?test4) then
   echo Compiling and linking testEngConsAllActinides
   $CC  testEngConsAllActinides.c -g -c $INCLUDE
   $LD  testEngConsAllActinides.o -g -o testEngConsAllActinides $EXTRAFCFLAGS $INCLUDE $LIBS

   echo "running testEngConsAllActinides (takes 56 seconds on aztec)"
   if (-e testEngConsAllActNudist.res) rm -f testEngConsAllActNudist.res
   if (-e testEngConsAllActSpec.res) rm -f testEngConsAllActSpec.res
   if (-e testEngConsAllActSpecTotal.res) rm -f testEngConsAllActSpecTotal.res
   if (-e testEngConsAllActTotal.res) rm -f testEngConsAllActTotal.res
   time ./testEngConsAllActinides
   diff -s testEngConsAllActNudist.res testEngConsAllActNudist.goldref
   diff -s testEngConsAllActSpec.res testEngConsAllActSpec.goldref
   diff -s testEngConsAllActSpecTotal.res testEngConsAllActSpecTotal.goldref
   diff -s testEngConsAllActTotal.res testEngConsAllActTotal.goldref
endif

if ($?test5) then
   echo Compiling and linking testPhotofission
   $CC  testPhotofission.c -g -c $INCLUDE
   $LD  testPhotofission.o -g -o testPhotofission $EXTRAFCFLAGS $INCLUDE $LIBS

   echo "running testPhotofission (takes 39 seconds on aztec)"
   if (-e testPhotofissionNuDist.res) rm -f testPhotofissionNuDist.res
   if (-e testPhotofissionSpec.res) rm -f testPhotofissionSpec.res
   time ./testPhotofission
   diff -s testPhotofissionNuDist.res testPhotofissionNuDist.goldref
   diff -s testPhotofissionSpec.res testPhotofissionSpec.goldref
   echo
endif

if ($?test6) then
   which $FC > /dev/null || echo Need $FC to link FREYA test && exit
   echo Compiling and linking testFreya
   $CC  testFreya.c -g -O3 -c $INCLUDE
   $LD  testFreya.o -g -O3 -o testFreya $EXTRAFCFLAGS $INCLUDE $LIBS

   setenv FREYADATAPATH ../data_freya
   echo "set FREYADATAPATH=$FREYADATAPATH"
   echo "running testFreya (takes 7 seconds on aztec)"
   if (-e testFreyaSpec.res) rm -f testFreyaSpec.res
   if (-e testFreyaSpecTotal.res) rm -f testFreyaSpecTotal.res
   if (-e testFreyaTotal.res) rm -f testFreyaTotal.res
   time ./testFreya
   diff -s testFreyaNuDist.res testFreyaNuDist.goldref
   diff -s testFreyaSpec.res testFreyaSpec.goldref
   diff -s testFreyaSpecTotal.res testFreyaSpecTotal.goldref
   diff -s testFreyaTotal.res testFreyaTotal.goldref
   echo
endif
