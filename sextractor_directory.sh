#!/bin/bash

cd $1

for fname in `\ls *.txt`; do
    bname=`basename $fname .txt`
    echo "sextractoring directory for $bname"
    bash /home/jim/code/ngs_project/spotproc.sh $bname
done
