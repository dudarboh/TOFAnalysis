#!/bin/bash
for i in 1 2 3 4 5 6
do
    cd /afs/desy.de/user/d/dudarboh/dust/job_files/jobs
    echo "Freshly baked ROOT files from the jobber:" `ls -1 *.root | wc -l`
    echo "Merge 500 of them to the result${i}.root file"
    find . -maxdepth 1 -type f |head -500|xargs mv -t /afs/desy.de/user/d/dudarboh/dust/final_files/2f_Z_hadronic/

    cd /afs/desy.de/user/d/dudarboh/dust/final_files/2f_Z_hadronic/

    hadd result${i}.root *.root
    mv result${i}.root ./results
    rm *.root
done
