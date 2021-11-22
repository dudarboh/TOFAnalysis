#!/bin/bash
source source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh && \
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/TOFAnalysis/lib/libTOFAnalysis.so && \
rm -rf * && cp ${1} . && \
filename=`ls *.slcio` && \
Marlin /afs/desy.de/user/d/dudarboh/TOFAnalysis/xml/steer.xml --global.LCIOInputFiles="${filename}" && \
mv -f rename.root ../output/${2}.root && \
rm -f *.slcio
