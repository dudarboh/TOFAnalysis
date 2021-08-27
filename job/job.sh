#!/bin/bash

source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh &&
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/TOFAnalysis/lib/libSETAnalysis.so

rm -rf job${2}
mkdir job${2}
cd ./job${2}

cp ${1} .
filename=`ls *.slcio`

Marlin /afs/desy.de/user/d/dudarboh/TOFAnalysis/xml/steer.xml --global.LCIOInputFiles="${filename}"
mv SETAnalysis_RENAME.root ${2}.root
rm -f ../${2}.root
mv ${2}.root ..

cd ..
rm -rf job${2}
