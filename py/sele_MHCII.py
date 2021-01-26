#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from selenium import webdriver
import re
from Bio import SeqIO
from math import log2

VAXIJEN_TARGET = "virus"     # You need to replace the field
VAXIJEN_THRESHOLD = 0.5      # You need to replace the field

RESULT_FILENAME = "NP_A.txt" # You need to replace the field


MSA_FILENAME = "NP_A_aln.fa" # You need to replace the field
MSA_REF_ID = "CY115151"      # You need to replace the field

AAS=['V','I','L','F','M','W','Y','C','A','T', 'D', 'E','G','P','R','K','H','N','Q','S','-']
CONSERVATION_THRESHOLD = 1.3

SEP = ";"

def chSep(v):
    a = v.split('.')
    return ','.join(a)

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
        for ar in AAS:
            q = 0
            if ar in col:
                q = col[ar] / count
                s -= q * log2(q)
            aafs.append(str(q))
        if aa != '-':
            ra.append([aa, (2**s) <= CONSERVATION_THRESHOLD])
    return ra

def pep_consv_count(st, end, consv):
    c = 0
    st1 = int(st) - 1
    end1 = int(end)
    count = float(end1 - st1)
    for i in range(st1, end1):
        if consv[i][1]:
            c += 1
    return int(c / count * 100)

def f1():
    d = {}
    d1 = {}
    consv = msa_stat()
    vjd = {}
    for l in open(RESULT_FILENAME).readlines():
        ta = l[:-2].split('\t')
        seq = ta[4]
        if not seq in d:
            d[seq] = []
        d[seq].append(ta)
    for k in d.keys():
        a = d[k]
        core = a[0][7]
        if not core in d1:
            d1[core] = []
        d1[core].append(k)
    coreOrdered = {}
    for core in d1.keys():
        pepsOredered = {}
        mic50 = 0
        peps = d1[core]
        for pep in peps:
            idtPeps = d[pep]
            ic50 = 0
            for a in idtPeps:
                ic50 += float(a[8])
            ic50 /= len(idtPeps)
            mic50 += ic50
            if not ic50 in pepsOredered:
                pepsOredered[ic50] = []
            pepsOredered[ic50].append(idtPeps)
        mic50 /= len(peps)
        if not mic50 in coreOrdered:
            coreOrdered[mic50] = []
        coreOrdered[mic50].append(pepsOredered)
    num = 1
    for mic50 in sorted(coreOrdered.keys()):
        a = coreOrdered[mic50]
        for pepsOredered in a:
            print ('')
            for ic50 in sorted(pepsOredered.keys()):
                idtPepsL = pepsOredered[ic50]
                for idtPeps in idtPepsL:
                    print ('')
                    a1 = idtPeps[0]
                    cv = pep_consv_count(int(a1[1]), int(a1[2]), consv)
                    pep = a1[4]
                    vj = '0'
                    if cv == 100:
                        if pep in vjd:
                            vj = vjd[pep]
                        else:
                            vj = str(vaxijen_test(pep))
                            vjd[pep] = vj
                        vj = chSep(vj)
                    print (str(num) + SEP + a1[0] + SEP + a1[1] + SEP + a1[2] + SEP + a1[3] + SEP + pep + SEP + chSep(a1[5]) + SEP + chSep(a1[6]) + SEP + a1[7] + SEP + chSep(a1[8]) + SEP + chSep(a1[9]) + SEP + str(cv) + SEP + vj)
                    num += 1
                    a3 = idtPeps[1:]
                    for a2 in a3:
                        print (str(num) + SEP + a2[0] + SEP +SEP + SEP + SEP + SEP + chSep(a2[5]) + SEP + chSep(a2[6]) + SEP + SEP + chSep(a2[8]) + SEP + chSep(a2[9]))
                        num += 1

f1()
