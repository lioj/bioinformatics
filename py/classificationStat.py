#!/usr/bin/env python

from Bio import SeqIO
import numpy as np

LEN_TH = 200

bl_out_fn = "txt/38.txt"


def calc(a):
    return np.amin(a), np.median(a), np.amax(a)

def stat():
    bl_sd_ids = {}
    in_list = list(open(bl_out_fn))
    count = 0
    for l in in_list:
        a = l[:-1].split('\t')
        l = int(a[3])
        if l >= LEN_TH:
            sacca = a[1].split('_')
            l1 = len(sacca)
            sacc = sacca[l1 - 2] + '_' + sacca[l1 - 1]
            eval = float(a[2])
            pidt = float(a[4])
            if not bl_sd_ids.has_key(sacc):
                bl_sd_ids.setdefault(sacc, [[],[],[]])
            d = bl_sd_ids[sacc]
            d[0].append(l)
            d[1].append(pidt)
            d[2].append(eval)
            count += 1
    print len(in_list), count, len(bl_sd_ids.keys())
    print "acc", "count", "len_min", "len_med", "len_max", "pidt_min", "pidt_med", "pidt_max",\
        "eval_min", "eval_med", "eval_max"
    for acc in bl_sd_ids.keys():
        a = bl_sd_ids[acc]
        lc = calc(a[0])
        pdtc = calc(a[1])
        evalc = calc(a[2])
        print acc, len(a[0]), lc[0], lc[1], lc[2], pdtc[0], pdtc[1], pdtc[2], evalc[0], evalc[1], evalc[2]



stat()
