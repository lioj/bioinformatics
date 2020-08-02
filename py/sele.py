#!/usr/bin/env python

from Bio import SeqIO

LEN_TH = 200

target = "38"

sele_fn = "txt/sele_" + target + "_ids.txt"
bl_out_fn = "txt/" + target + ".txt"
in_fq_fn = "fq/" + target + ".fq"
out_fq_fn = "fq/sele_" + target + ".fq"

def sele():
    sele_ids = {}
    bl_sd_ids = {}
    out_rs = []
    for l in open(sele_fn):
        sele_ids.setdefault(l[:-1], '')
    for l in open(bl_out_fn):
        a = l[:-1].split('\t')
        if sele_ids.has_key(a[1]) and int(a[3]) >= LEN_TH:
            bl_sd_ids.setdefault(a[0], '')
    print len(bl_sd_ids.keys())
    for r in SeqIO.parse(in_fq_fn, "fastq"):
        if bl_sd_ids.has_key(r.id):
            out_rs.append(r)
    print len(out_rs)
    SeqIO.write(out_rs, out_fq_fn, "fastq")



sele()
