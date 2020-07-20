#!/usr/bin/env python

from Bio import SeqIO
from math import log2

fn = "HA_B_aln.fa"
ref = "CY115151"
sep = ","
aas=['V','I','L','F','M','W','Y','C','A','T', 'D', 'E','G','P','R','K','H','N','Q','S','-']



rs = ""
msa=[]
count = 0.0
for r in SeqIO.parse(fn, "fasta"):
    ids = r.id.split('|')
    s = str(r.seq)
    count += 1
    if ids[0] == ref:
        rs = s
    for i, aa in enumerate(s):
        if len(msa) < i + 1:
            msa.append({})
        d = msa[i]
        if not aa in d:
            d[aa] = 0
        d[aa] += 1
print('#' + sep + 'aa' + sep + 'Entropy' + sep + '2^Entropy' + sep + sep.join(aas))
en=0
for i, aa in enumerate(rs):
    if aa != '-':
        en += 1
    col = msa[i]
    aafs=[]
    s = 0
    for ar in aas:
        q = 0
        if ar in col:
            q = col[ar] / count
            s -= q * log2(q)
        aafs.append(str(q))
    print(str(en) + sep + aa + sep + str(s) + sep + str(2**s) + sep + sep.join(aafs))
