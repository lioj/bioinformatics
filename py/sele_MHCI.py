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

RESULT_FILENAME = "txt/mhci_ns5_denv4.txt" # You need to replace the field


MSA_FILENAME = "aln/st/denv4_NS5_aln.fasta" # You need to replace the field
MSA_REF_ID = "NP_740325.1"      # You need to replace the field

AAS=['V','I','L','F','M','W','Y','C','A','T', 'D', 'E','G','P','R','K','H','N','Q','S','-']
CONSERVATION_THRESHOLD = 1.3

SEP = ";"

def chSep(v):
    a = v.split('.')
    return ','.join(a)

def vaxijen_test(aa_seq, browser):
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
    browser = webdriver.Chrome()
    for i,l in enumerate(open(RESULT_FILENAME).readlines()):
        ta = l[:-2].split('\t')
        if i > 0:
            seq = ta[5]
            cv = pep_consv_count(int(ta[2]), int(ta[3]), consv)
            vj = 0
            if cv == 100:
                if seq in vjd:
                    vj = vjd[seq]
                else:
                    vj = str(vaxijen_test(seq, browser))
                    vjd[seq] = vj
            vj = ','.join(str(vj).split('.'))
            print(ta[0]+';'+ta[1]+';'+ta[2]+';'+ta[3]+';'+ta[4]+';'+ta[5]+';'+ta[6]+';'+ta[7]+';'+ta[8]+';'+ta[9]+';' + str(cv) + ';' + vj)
        else:
            print(ta[0]+';'+ta[1]+';'+ta[2]+';'+ta[3]+';'+ta[4]+';'+ta[5]+';'+ta[6]+';'+ta[7]+';'+ta[8]+';'+ta[9]+';consv;vj')
    browser.close()

f1()

#print(msa_stat())
