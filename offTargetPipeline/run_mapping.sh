#!/bin/bash

FOLDER="/Users/ramil.mintaev/bioinf/crispr"
WDIR=$FOLDER/$1
READFOLDER="trim"
READ1suf="-R1"
READ2suf="-R2"

rm -rf $WDIR/sam
mkdir $WDIR/sam

for NAME in '1' '2' 'c' 'c2'; do
    echo $NAME
    bowtie2 -x $FOLDER/index/sites --threads 4 -1 $WDIR/in/$READFOLDER/$NAME$READ1suf.fq.gz -2 $WDIR/in/$READFOLDER/$NAME$READ2suf.fq.gz -S $WDIR/sam/$NAME.sam
    samtools view -bS $WDIR/sam/$NAME.sam > $WDIR/bam/$NAME.bam
    samtools sort $WDIR/bam/$NAME.bam > $WDIR/bam/$NAME"_s".bam
    samtools index $WDIR/bam/$NAME"_s".bam
    rm -r $WDIR/bam/$NAME.bam
done

rm -r $WDIR/sam
