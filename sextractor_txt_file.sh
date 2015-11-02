#!/bin/bash

bname=`basename $1 .txt`
echo $bname
spotproc.sh $bname
