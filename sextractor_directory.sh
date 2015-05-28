#!/bin/bash

cd $1

for fname in `\ls *.txt`; do
    bname=`basename $fname .txt`
    echo $bname
    spotproc.sh $bname
done
