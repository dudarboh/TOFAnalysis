#!/bin/bash
. /cvmfs/ilc.desy.de/sw/x86_64_gcc82_sl6/v02-02/init_ilcsoft.sh
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/lib/libTOFAnalysis.so

# $1 is a job number. It defines name of the output root file and list of slcio files.
# It takes 35 files per job. So 12k files in total fit in 350 jobs

#n files per job
n_files=35
start_job=0
job_number=`echo "${start_job} + ${1}" | bc`

rm -rf job${job_number}
mkdir job${job_number}
cd ./job${job_number}

start_file=`echo "$n_files * $job_number + 1" | bc`
echo "Start_file: $start_file"

end_file=`echo "$n_files * ($job_number + 1)" | bc `
echo "End_file: $end_file"

slcio_files=`sed -n "$start_file , $end_file p" /afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/job/2f_Z_hadronic.txt`
slcio_files=`ls -1 $slcio_files | tr '\n' ' '`

sleep `echo "0.01 * $job_number" | bc`

Marlin /afs/desy.de/user/d/dudarboh/ILCSoft/TOFAnalysis/xml/steer.xml --global.LCIOInputFiles="${slcio_files}" --TOFAnalysis.output_filename="${job_number}.root"

rm -f ../${job_number}.root
mv ${job_number}.root ..
cd ..
rm -rf job${job_number}
