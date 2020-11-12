#!/usr/bin/env python

from selenium import webdriver
import re
import numpy as np
from Bio import SeqIO
from math import log2

VAXIJEN_TARGET = "virus"                  # You need to replace the field
VAXIJEN_THRESHOLD = 0.5                   # You need to replace the field

BEPIPRED_RESULT_FILENAME = "BP_HA_B.csv"  # You need to replace the field
BEPIPRED_RESULT_HAS_HEADLINE = True
BEPIPRED_THRESHOLD = 0.5                  # You need to replace the field


MSA_FILENAME = "HA_B_aln.fa"              # You need to replace the field
MSA_REF_ID = "CY115151"                   # You need to replace the field


AAS=['V','I','L','F','M','W','Y','C','A','T', 'D', 'E','G','P','R','K','H','N','Q','S','-']
CONSERVATION_THRESHOLD = 1.3

SEP = ","

def vaxijen_test(aa_seq):
    browser = webdriver.Chrome()
    browser.get("http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html")
    seq_f = browser.find_element_by_name("seq")
    seq_f.send_keys(aa_seq)
    target = browser.find_element_by_name("Target")
    all_options = target.find_elements_by_tag_name("option")
    for option in all_options:
        if option.get_attribute("value") == VAXIJEN_TARGET:
            option.click()
    threshold = browser.find_element_by_name("threshold")
    threshold.send_keys(str(VAXIJEN_THRESHOLD))
    browser.find_element_by_name("submit").click()
    r = re.findall(r'=\<b\> ([-\.0-9]+) \<\/b\>', browser.page_source)
    browser.close()
    return float(r[0])

def parse_bepipred_file():
    in_list = list(open(BEPIPRED_RESULT_FILENAME))
    r=[]
    for i,l in enumerate(in_list):
        if i > 0 or not BEPIPRED_RESULT_HAS_HEADLINE:
            a = l[:-1].split(SEP)
            r.append([a[2], float(a[-1])])
    return r

def msa_stat():
    ra=[]
    rs = ""
    msa=[]
    count = 0.0
    for r in SeqIO.parse(MSA_FILENAME, "fasta"):
        ids = r.id.split('|')
        s = str(r.seq)
        count += 1
        if ids[0] == MSA_REF_ID:
            rs = s
        for i, aa in enumerate(s):
            if len(msa) < i + 1:
                msa.append({})
            d = msa[i]
            if not aa in d:
                d[aa] = 0
            d[aa] += 1
    for i, aa in enumerate(rs):
        col = msa[i]
        aafs=[]
        s = 0
        for ar in aas:
            q = 0
            if ar in col:
                q = col[ar] / count
                s -= q * log2(q)
            aafs.append(str(q))
        if aa != '-':
            ra.append([aa, (2**s) <= CONSERVATION_THRESHOLD])
    return ra

def predict_simple_eps():
    bp = parse_bepipred_file()
    consv = parse_msa_stat_file()
    eps=[]
    ep = ""
    fp = 0
    bps = []
    for i, e in enumerate(bp):
        cv = consv[i]
        pos = i + 1
        aa = e[0]
        bp_prob = e[1]
        is_consv = cv[1]
        l = len(ep)
        if bp_prob >= BEPIPRED_THRESHOLD:
            if l == 0:
                fp = pos
            if is_consv:
                ep += aa
            else:
                ep += aa.lower()
            bps.append(bp_prob)
        elif l > 0:
            if l > 6:
                eps.append([fp, ep, np.median(bps), vaxijen_test(ep)])
            ep = ""
            bps = []
    print("pos"+SEP+"len"+SEP+"seq"+SEP+"VaxiJen"+SEP+"BepiPredMedian")
    for e in eps:
        s = e[1]
        print(str(e[0]) + SEP + str(len(s)) + SEP + s + SEP + str(e[2]) + SEP + str(e[3]))

def opt_vaxijen_seq(seq, s):
    l = len(seq)
    if l > 6 and s > 6:
        for i in range(0, l):
            e = i + s
            if e <= l:
                s1 = seq[i:e]
                v = vaxijen_test(s1)
                if v >= VAXIJEN_THRESHOLD:
                    return s1, i, round(v, 4)
        return opt_vaxijen_seq(seq, s - 1)
    return "", 0, 0

def bp_median(a, fp, s):
    bps = []
    p = fp - 1
    for i in range(p, p + s):
        bps.append(a[i][1])
    return round(np.median(bps), 4)

def predict_opt_eps():
    bp = parse_bepipred_file()
    consv = msa_stat()
    eps=[]
    ep = ""
    fp = 0
    for i, e in enumerate(bp):
        cv = consv[i]
        pos = i + 1
        aa = e[0]
        bp_prob = e[1]
        is_consv = cv[1]
        l = len(ep)
        if bp_prob >= BEPIPRED_THRESHOLD and is_consv:
            if l == 0:
                fp = pos
            ep += aa
        elif l > 0:
            if l > 6:
                vs = opt_vaxijen_seq(ep, len(ep))
                if vs[0] != '':
                    fp1 = fp + vs[1]
                    s1 = vs[0]
                    vj = vs[2]
                    bpm = bp_median(bp, fp1, len(s1))
                    eps.append([fp1, s1, vj, bpm])
            ep = ""
            bps = []
    print("pos"+SEP+"len"+SEP+"seq"+SEP+"VaxiJen"+SEP+"BepiPredMedian")
    for e in eps:
        s = e[1]
        print(str(e[0]) + SEP + str(len(s)) + SEP + s + SEP + str(e[2]) + SEP + str(e[3]))

predict_opt_eps()
