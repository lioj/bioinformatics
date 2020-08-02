#!/usr/bin/env python

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

bam_fn = "bam/38_s.bam"
cons_fn = "fa/38_cons.fa"

def bases(pileups):
    d={}
    for pileupread in pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            b = pileupread.alignment.query_sequence[pileupread.query_position]
            if not d.has_key(b):
                d[b] = 0
            d[b] += 1
    return d

def opt_base(dic):
    b = "-"
    v = 0
    for k in dic.keys():
        if dic[k] > v:
            v = dic[k]
            b = k
    return b

def f1():
    samfile = pysam.AlignmentFile(bam_fn, "rb")
    segs={}
    for pileupcolumn in samfile.pileup():
        d = bases(pileupcolumn.pileups)
        if not segs.has_key(pileupcolumn.reference_name):
            segs.setdefault(pileupcolumn.reference_name, {})
        segs[pileupcolumn.reference_name].setdefault(pileupcolumn.reference_pos, opt_base(d))
    out_rs = []
    for sname in segs.keys():
        seg = segs[sname]
        poss = seg.keys()
        poss.sort()
        print sname, poss
        seq = ""
        i = 0
        for p in poss:
            if p != i:
                i += 1
                seq += "-"
            seq += seg[p]
            i += 1
        out_rs.append(SeqRecord(Seq(seq), id = sname,
            description = "consensus was made by python scr"))
    SeqIO.write(out_rs, cons_fn, "fasta")

f1()
