#!/bin/bash
. /cvmfs/ilc.desy.de/sw/x86_64_gcc82_sl6/v02-02/init_ilcsoft.sh
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/lib/libTOFAnalysis.so

# $1 is a job number. It defines name of the output root file and list of slcio files.
# It takes 35 files per job. So 12k files in total fit in 350 jobs

rm -rf job${1}
mkdir job${1}
cd ./job${1}

#n files per job
n_files=35
start_file=$((${n_files} * ${1} + 1))
end_file=$((${n_files} * (${1} + 1)))

slcio_files=`sed -n "$start_file , $end_file p" /afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/job/2f_Z_hadronic.txt`
slcio_files=`ls -1 $slcio_files | tr '\n' ' '`

Marlin /afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/xml/steer.xml --global.LCIOInputFiles="${slcio_files}" --TOFAnalysis.output_filename="${1}.root"

rm -f ../${1}.root
mv ${1}.root ..
cd ..
rm -rf job${1}
