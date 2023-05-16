#!/usr/bin/python

import pysam
from gzip import open as gzopen
from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
from sys import argv
from multiprocessing import Pool, Queue
import re

script_name, WDIR = argv

WDIR += "/"

def init_pool_processes(the_queue):
    global queue
    queue = the_queue

def load_bam(sample_name):
    return pysam.AlignmentFile(WDIR + "bam/" + sample_name + "_s.bam",\
     "rb")

def gids(id):
    ids = id.split('|')[0].split('_')
    i = 0
    if len(ids) > 1:
        i = int(ids[1]) - 1
    return ids[0], i

def window_seq(seq, st, wd):
    s = ""
    ls = len(seq)
    stp = wd[0] - st + 1
    ep = wd[1] - st + 1
    if stp >= 0 and ls > 0 and ep <= ls:
        s = seq[stp:ep]
    return s

def find_diff(seq):
    a = []
    for i, s in enumerate(seq):
        if s != s.upper():
            a.append(i)
    return a

def parse_cigar(cigar):
    ciga = []
    a = re.findall(r'\d+[MIDNSHP=X]?', cigar)
    for cig in a:
        cig_value = int(cig[:-1])
        cig_type = cig[-1:]
        if cig_type == "S":
            raise Exception("parse_cigar: There is an unexpected soft clipping in cigar string!")
        ciga.append([cig_value, cig_type])
    return ciga

def find_read_sq(r, wd, p):
    ciga = parse_cigar(r.cigarstring)
    rd_st, rd_e = -1, -1
    ref_st, ref_e = r.reference_start - 1, r.reference_start - 1
    for cp in ciga:
        if cp[1] == 'M':
            rd_st = rd_e + 1
            rd_e = rd_st + cp[0]
            ref_st = ref_e + 1
            ref_e = ref_st + cp[0]
        elif cp[1] == 'I':
            rd_st = rd_e + 1
            rd_e = rd_st + cp[0]
            #ref_st = ref_e + 1
            #ref_e = ref_st + 1
        elif cp[1] == 'D':
            ref_st = ref_e + 1
            ref_e = ref_st + 1
        else:
            raise Exception("Unexpected cigar: " + cp[1])
        pos = wd[0] + p + 1
        if pos >= ref_st and pos < ref_e:
            pos = rd_st + r.query_alignment_start + pos - ref_st
            #print("Yes")
            return r.query_sequence[pos:pos+1],\
                r.query_qualities[pos:pos+1]
    #print("B", pos, ref_st, ref_e, r.reference_start)
    return "B", [-1]
    raise Exception("NOT FOUND: rd_st:" + str(rd_st) + " rd_e:"  + str(rd_e) + " ref_st:" + str(ref_st) + " ref_e:" + str(ref_e) + " wd[0]:" + str(wd[0]) + " wd[1]:" + str(wd[1]) + " p:" + str(p) + " cigar:" + r.cigarstring)

def wd_read_sq(r, wd, dual_diffs):
    ciga = parse_cigar(r.cigarstring)
    stp = wd[0] - r.reference_start + 1 + r.query_alignment_start
    ep = wd[1] - r.reference_start + 1 + r.query_alignment_start
    t, q = "", []
    for p in range(wd[1] - wd[0]):
        sq = find_read_sq(r, wd, p)
        t += sq[0]
        try:
            q.append(sq[1][0])
        except:
            q.append(15)
            print("!!!!!!", sq, r, wd, dual_diffs, file=sys.stderr)
    #print(t, q)
    #print("cigar", ciga, dual_diffs, q[dual_diffs[0]])
    for di in dual_diffs:
        if q[di] < 10:
            return False
    #print("reference_start", r.reference_start)
    #print("is_forward", r.is_forward)
    #print("stp ep", stp, ep)
    #print("wd", wd)
    return True
    

def has_mut(bam, r, wd):
    if not r.has_tag('MD'):
        raise Exception("has_indel: There is no MD tag in the read!")
    ref_st1 = r.reference_start
    wd_seq1 = window_seq(r.get_reference_sequence(), ref_st1, wd)
    #quals1 = r.query_qualities
    if wd_seq1 != "" and wd_seq1 == wd_seq1.upper():
        return False
    diffs1 = find_diff(wd_seq1)
    #
    dual_diffs = []
    mr = bam.mate(r)
    ref_st2 = mr.reference_start
    wd_seq2 = window_seq(mr.get_reference_sequence(), ref_st2, wd)
    if wd_seq2 != "" and wd_seq2 == wd_seq2.upper()\
     or (r.mapping_quality <= 10 and mr.mapping_quality <= 10):
        return False
    diffs2 = find_diff(wd_seq2)
    for di1 in diffs1:
        if di1 in diffs2:
            dual_diffs.append(di1)
    if len(dual_diffs) == 0:
        return False
    #NEED CONTINUE
    #print("!")
    #print("wd_seq1", wd_seq1)
    #print("wd_seq2", wd_seq2)
    if not wd_read_sq(r, wd, dual_diffs) or\
        not wd_read_sq(mr, wd, dual_diffs):
        return False
    return True

def norm_size(size1, size2, thr):
    s1=size1
    s2=size2
    if size1 > thr:
        s1 = thr
    if size2 > thr:
        s2 = thr
    if s1 < thr:
        s2 += thr - s1
    if s2 < thr:
        s1 += thr - s2
    return int(s1), int(s2)

def apply_thrs(d):
    rd = {}
    for site in d.keys():
        td = d[site]
        size_1 = 0
        size_2 = 0
        size_c = 0
        size_c2 = 0
        if "1" in td:
            size_1 = td["1"]
        if "2" in td:
            size_2 = td["2"]
        if "c" in td:
            size_c = td["c"]
        if "c2" in td:
            size_c2 = td["c2"]
        size = size_1 + size_2
        c_size = size_c + size_c2
        if size > c_size:
            size = c_size
        size /= 2
        size_1, size_2 = norm_size(size_1, size_2, size)
        size_c, size_c2 = norm_size(size_c, size_c2, size)
        rd.setdefault(site, {"1" : size_1, "2" : size_2,\
            "c" : size_c, "c2" : size_c2})
    return rd

class mainCl():
    def __init__(self):
        self.window_size = 5
        self.sample_names = ["1", "2", "c", "c2"]
        self.load_off_target_info()

    def load_off_target_info(self):
        # extr sites
        self.sites = {}
        for l in open("seq/sites.txt", 'r').readlines():
            a = l[:-1].split("\t")
            if not a[0] in self.sites:
                self.sites.setdefault(a[0], [])
            self.sites[a[0]].append(a[1].upper())
        # extr double break info
        self.dbs = {}
        for r in SeqIO.parse("seq/site_regions.fa", "fasta"):
            spacer, i = gids(r.id)
            site = self.sites[spacer][i]
            r.seq = Seq(str(r.seq).upper())
            seq = r.seq
            pos = seq.find(site)
            is_forward = True
            is_cas9 = True
            if r.id[0] == "c":
                is_cas9 = False
            if pos == -1:
                is_forward = False
                seq = seq.reverse_complement()
                pos = len(seq) - seq.find(site)
                if is_cas9:
                    pos -= 17
                else:
                    pos -= len(site) + 1
            else:
                if is_cas9:
                    pos += 15
                else:
                    pos += 22
            self.dbs.setdefault(r.id, [pos, r, is_forward, is_cas9])

    def start(self):
        self.reads_balancing()
        stats = self.count_mutations()
        print("site indels regular c_indels c_regular")
        for site in stats.keys():
            sdic = stats[site]
            t = ""
            for sample_names in [self.sample_names[:2], self.sample_names[2:]]:
                mut_s = 0
                wt_s = 0
                for sample_name in sample_names:
                    if sample_name in sdic:
                        vdic = sdic[sample_name]
                        mut_s += vdic["mutated"]
                        wt_s += vdic["total"] - vdic["mutated"]
                mut_s /= 2
                wt_s /= 2
                t += str(mut_s) + " " + str(wt_s) + " "
            print(site, t[:-1])

    def reads_balancing(self):
        print("reads_balancing", file=sys.stderr)
        queue = Queue()
        with Pool(initializer=init_pool_processes, initargs=(queue,)) as p:
            p.map(self.sele_reads, self.sample_names)
        counted = {}
        while not queue.empty():
            sample_name, d = queue.get()
            for site in d.keys():
                if site not in counted:
                    counted.setdefault(site, {})
                counted[site].setdefault(sample_name, d[site])
        queue.close()
        queue.join_thread()
        self.thrs = apply_thrs(counted)
        for site in self.thrs:
            t = ""
            d = self.thrs[site]
            for sample in ["1", "2", "c", "c2"]:
                if sample in d:
                    t += str(d[sample]) + " "
                else:
                    t += "None "
            print(site, t, file=sys.stderr)

    def sele_reads(self, sample_name):
        print("started:", sample_name, file=sys.stderr)
        bam = load_bam(sample_name)
        d = {}
        for r in bam.fetch():
            #if r.reference_name != "R_3|KIAA1614":
            #    continue
            if r.is_mapped and r.mate_is_mapped:
                if not r.reference_name in d:
                    d.setdefault(r.reference_name, 0)
                d[r.reference_name] += 1
        print("ended:", sample_name, file=sys.stderr)
        queue.put((sample_name, d))

    def count_mutations(self):
        print("mutations counting", file=sys.stderr)
        queue = Queue()
        with Pool(initializer=init_pool_processes, initargs=(queue,)) as p:
            p.map(self.anal, self.sample_names)
        print("parsing", file=sys.stderr)
        counted = {}
        while not queue.empty():
            sample_name, d = queue.get()
            for site in d.keys():
                if site not in counted:
                    counted.setdefault(site, {})
                counted[site].setdefault(sample_name, d[site])
        queue.close()
        queue.join_thread()
        return counted

    def anal(self, sample_name):
        print("started:", sample_name, file=sys.stderr)
        bam = load_bam(sample_name)
        d = {}
        cd = {}
        for r in bam.fetch():
            #if r.reference_name != "R_3|KIAA1614":
            #    continue
            if r.is_mapped and r.mate_is_mapped:
                thr_size = self.thrs[r.reference_name][sample_name]
                if not r.reference_name in cd:
                    cd.setdefault(r.reference_name, 0)
                    d.setdefault(r.reference_name,\
                        {"mutated":0, "total":0})
                if cd[r.reference_name] > thr_size:
                    continue
                cd[r.reference_name] += 1
                #
                td = d[r.reference_name]
                td["total"] += 1
                sa = r.get_cigar_stats()[0]
                if len(sa) == 11 and (sa[1] > 0 or sa[2] > 0):
                    db = self.dbs[r.reference_name]
                    wd = (db[0] - self.window_size, db[0] + self.window_size)
                    if has_mut(bam, r, wd):
                        td["mutated"] += 1
                        #break
        queue.put((sample_name, d))
        print("ended:", sample_name, file=sys.stderr)


if __name__ == '__main__':
    mc = mainCl()
    mc.start()
