#!/usr/bin/python 

import pysam
from gzip import open as gzopen
from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
from sys import argv
from multiprocessing import Pool, Queue

script_name, WDIR = argv

WDIR += "/"

def init_pool_processes(the_queue):
    global queue
    queue = the_queue

def load_bam(sample_name):
    return pysam.AlignmentFile(WDIR + "bam/" + sample_name + "_s.bam", "rb")

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
    if r.mate_is_mapped:
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
    elif r.mapping_quality <= 10:
        return False
    #NEED CONTINUE
    return True

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
        queue = Queue()
        with Pool(initializer=init_pool_processes, initargs=(queue,)) as p:
            p.map(self.anal, self.sample_names)
        p.join()
        stats = {}
        print("parsing", file=sys.stderr)
        while not queue.empty():
            sample_name, d = queue.get()
            for site in d.keys():
                if site not in stats:
                    stats.setdefault(site, {})
                stats[site].setdefault(sample_name, d[site])
        queue.close()
        queue.join_thread()
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

    def anal(self, sample_name):
        print("started:", sample_name, file=sys.stderr)
        bam = load_bam(sample_name)
        d = {}
        for r in bam.fetch():
            if r.is_mapped:
                if not r.reference_name in d:
                    d.setdefault(r.reference_name, {"mutated":0, "total":0})
                td = d[r.reference_name]
                td["total"] += 1
                sa = r.get_cigar_stats()[0]
                if len(sa) == 11 and (sa[1] > 0 or sa[2] > 0):
                    db = self.dbs[r.reference_name]
                    wd = (db[0] - self.window_size, db[0] + self.window_size)
                    if has_mut(bam, r, wd):
                        td["mutated"] += 1
        queue.put((sample_name, d))
        print("ended:", sample_name, file=sys.stderr)


if __name__ == '__main__':
    mc = mainCl()
    mc.start()

