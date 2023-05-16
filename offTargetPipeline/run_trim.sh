#!/bin/bash

FOLDER="/Users/ramil.mintaev/bioinf/crispr"
INDIR=$1
OUTDIR=$2
READ1suf="_R1"
READ2suf="_R2"

for NAME in '1' '2' 'c' 'c2'; do
   java -jar /Users/ramil.mintaev/bioinf/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE $FOLDER/$INDIR/$NAME$READ1suf.fq.gz $FOLDER/$INDIR/$NAME$READ2suf.fq.gz $FOLDER/$OUTDIR/in/trim/$NAME-R1.fq.gz $FOLDER/$OUTDIR/in/trim/$NAME-R1_unpaired.fq.gz $FOLDER/$OUTDIR/in/trim/$NAME-R2.fq.gz $FOLDER/$OUTDIR/in/trim/$NAME-R2_unpaired.fq.gz MINLEN:30
done
