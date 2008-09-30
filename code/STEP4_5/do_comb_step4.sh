#!/bin/sh
if [ $# -lt 2 ]
then echo "Usage: $0 year month1 month2";
     echo "       to add some months to SBBX.HadR2"
     exit; fi

fortran_compile=$FC
if [ "$FC" = "" ]
then echo "set an environment variable FC to the fortran_compile_command like f90"
     echo "or do all compilation first and comment the compilation lines"
     exit
fi

if [ ! -s input/SBBX.HadR2 ]
then echo "file input/SBBX.HadR2 not found" ; exit ; fi

mo=$2 ; if [ mo -lt 10 ] ; then mo='0'${mo} ; fi
if [ ! -s input/oiv2mon.$1$mo ]
then echo "file input/oiv2mon.$1$mo not found" ; exit ; fi

mo2=$2 ; if [ $# -gt 2 ]; then mo2=$3 ; fi

$FC code/STEP4_5/convert1.HadR2_mod4.f -o bin/convert1.HadR2_mod4.exe
bin/convert1.HadR2_mod4.exe $1 $2 $mo2
# Replacing trimSBBX work/SBBX.HadR2.upd

ln work/SBBX.HadR2.upd fort.10
${fortran_compile} code/STEP4_5/trimSBBX.f -o bin/trimSBBX.exe
bin/trimSBBX.exe
mv work/fort.11 work/SBBX.HadR2.upd.trim
rm -f work/fort.10


echo 'If all went ok, type "mv work/SBBX.HadR2.upd.trim input/SBBX.HadR2"'
echo '          then  type: code/STEP4_5/do_comb_step5.sh 100 0'
