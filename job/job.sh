#!/bin/bash

source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh &&
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/SETAnalysis/lib/libSETAnalysis.so

rm -rf job${2}_${3}
mkdir job${2}_${3}
cd ./job${2}_${3}

cp ${1} .
filename=`ls *.slcio`

Marlin /afs/desy.de/user/d/dudarboh/SETAnalysis/xml/steer.xml --global.LCIOInputFiles="${filename}" --MySETAnalysis.smearing=${3}
mv SETAnalysis_RENAME.root ${2}_${3}.root
rm -f ../${2}_${3}.root
mv ${2}_${3}.root ..

cd ..
rm -rf job${2}_${3}
