#!/bin/bash
cd /afs/desy.de/user/d/dudarboh/dust/job_files/jobs &&
mkdir /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp &&

echo "Freshly baked ROOT files from the jobber:" `ls -1 *.root | wc -l`
echo "Merge 1000 of them to the root file"

find . -maxdepth 1 -type f |head -1000|xargs mv -t /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp &&

cd /afs/desy.de/user/d/dudarboh/dust/final_files/SET/tmp &&
hadd rename.root *.root &&
mv rename.root .. &&
rm *.root &&
cd .. &&
rmdir tmp
