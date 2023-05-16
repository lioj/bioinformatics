#!/bin/bash

RUNF=$1
WDIR=work/$RUNF

echo "clean"
./scr/prod/clean.sh $WDIR

echo "trim input reads"
./scr/prod/run_trim.sh in/$RUNF $WDIR

echo "map reads"
export PATH=$PATH:/Users/ramil.mintaev/bioinf/bin/bowtie2-2.4.5-macos-x86_64/
./scr/prod/run_mapping.sh $WDIR

echo "diff calling"
python3 scr/prod/diff_calling_v2.py $WDIR > $WDIR/txt/indel_info.txt

echo "analyze"
./scr/prod/Ftest.R $WDIR

