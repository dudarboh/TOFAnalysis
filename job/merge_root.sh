#!/bin/bash
for i in 1
do
    cd /afs/desy.de/user/d/dudarboh/dust/job_files/jobs
    mkdir /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp
    echo "Freshly baked ROOT files from the jobber:" `ls -1 *.root | wc -l`
    echo "Merge 500 of them to the result${i}.root file"
    find . -maxdepth 1 -type f |head -500|xargs mv -t /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp

    cd /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp
    

    hadd result${i}.root *.root
    mv result${i}.root ..
    rm *.root
    cd .. && rmdir tmp
done
