#!/bin/bash
. /cvmfs/ilc.desy.de/sw/x86_64_gcc82_sl6/v02-02/init_ilcsoft.sh
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/lib/libTOFAnalysis.so

cd /nfs/dust/ilc/user/dudarboh/job_files
Marlin /afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/xml/TOFAnalysis.xml --global.LCIOInputFiles=${1} --TOFAnalysis.output_filename="${2}.root"

# $1 1st command line argument --slcio file name
# $2 2nd command line argument --output root file name
