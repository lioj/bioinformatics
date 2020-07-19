#!/usr/bin/env python

from Bio import SeqIO

LEN_TH = 400

in_fq_fn = "fq/38.fq"
out_fa_fn = "fa/38.fa"

def f1():
    out_rs = []
    count = 0
    sele = 0
    in_rs = list(SeqIO.parse(in_fq_fn, "fastq"))
    fcount = len(in_rs)
    for r in in_rs:
        count += 1
        if len(r.seq) > LEN_TH:
            sele += 1
            out_rs.append(r)
        print sele, count, fcount
    SeqIO.write(out_rs, out_fa_fn, "fasta")



f1()
