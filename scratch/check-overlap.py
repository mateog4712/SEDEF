# 786 

import pandas as pd
import numpy as np
import subprocess as sp
import sys, re, os, glob
from collections import *

def system(x):
    return sp.check_output(x, shell=True, executable="/bin/bash").strip()

tab_file = sys.argv[1]
path = sys.argv[2]
log = sys.argv[3]
if len(sys.argv) > 4:
    chrom1 = sys.argv[4]
    chrom2 = sys.argv[5]
    strand = sys.argv[6]
    strand = '_' if strand == 'y' else '+'
else:
    chrom1 = ''

df = pd.read_table(tab_file)
if chrom1 != '':
    if chrom1 != chrom2 or strand == '_':
        df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand)]
    else:
        df = df[(df.chrom == chrom1) & (df.otherChrom == chrom2) & (df.strand == strand) & (df.chromStart < df.otherStart)]
df['chromSize'] = df.chromEnd - df.chromStart
# print 'Loaded {} hits from WGAC'.format(df.shape[0])
hits = {}
name_to_coor = {}
with open(path + '_temp.bed', 'w') as f:
    for _, r in df.iterrows():
        if '_' in r.chrom or '_' in r.otherChrom:
            continue
        A = [r.chrom, r.chromStart, r.chromEnd,
            r.otherChrom, r.otherStart, r.otherEnd,
            r.alignfile, 0, '+', '-' if r.strand == '_' else '+']
        print >>f, '\t'.join(map(str, A))
        hits[r.alignfile] = list()
        name_to_coor[r.alignfile] = r #(r.chromStart, r.chromEnd, r.otherStart, r.otherEnd)
print 'WGAC:     {:6,} ({:,} lines)'.format(len(hits), df.shape[0])

def overlap(sa, ea, sb, eb):
    return max(0, min(ea, eb) - max(sa, sb))

def match(sW, eW, sS, eS):
    oo = overlap(sW, eW, sS, eS)
    wW = oo
    dW = 100.0 * wW / float(eW - sW)
    nE = (eS-sS+eW-sW-wW) / float(eS-sS) # how much extention is needed in pct of sedef's size
    return [(dW, wW, eW - sW, nE)]

def diff(wgac, sedef): # how much wgac is off sedef
    if not tuple(wgac[0:3]) < tuple(wgac[3:6]):
        wgac = wgac[3:6] + wgac[0:3] + wgac[6:]
    if not tuple(sedef[0:3]) < tuple(sedef[3:6]):
        sedef = sedef[3:6] + sedef[0:3] + sedef[6:]
    return match(*(wgac[1:3] + sedef[1:3])) + match(*(wgac[4:6] + sedef[4:6]))

def process_path(path, pnew):
    with open(pnew, 'w') as fw:
        with open(path) as f:
            for l in f:
                l = l.strip().split('\t')
                # s1, e1, s2, e2 = map(int, l[1:3] + l[4:6])
                # w = max(e1 - s1, e2 - s2) * 1
                # o = min(15000, w * 5)
                # l[1:3] = [max(1, s1 - o), e1 + o]
                # l[4:6] = [max(1, s2 - o), e2 + o]
                print >>fw, '\t'.join(map(str, l[:15]))

# path = path + "____"
process_path(path, path + '#')
if True or not os.path.exists(path + '_temp_diff.bed'):
    system("bedtools pairtopair -a {0}_temp.bed -b {0}# -type both > {0}_temp_diff.bed".format(path))
# os.unlink(path)
# print 'Processing {:,} bedtools lines'.format(system("wc -l {}_temp_diff.bed".format(path)))

fo = open(log, 'w')

with open(path + '_temp_diff.bed') as f:
    for l in f:
        l = l.strip().split()
        q = 0
        A = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5]), l[q+8]+l[q+9]) # WGAC
        q = 10
        B = (l[q], int(l[q+1]), int(l[q+2]), l[q+3], int(l[q+4]), int(l[q+5]), l[q+7]+l[q+8]) # SEDEF
        name = l[6]
        hits[name].append((diff(A, B), A, B))

missed_bases = 0
part_missed_bases = 0

try:
    tm = sum(1 for k, vs in hits.items() if len(vs) == 0)
    print 'Missed:  {:7,} hits ({:4.1f}%)'.format(tm, 100.0*tm/len(hits))
    for name, vs in sorted(hits.items()):
        if len(vs) == 0:
            r = name_to_coor[name]
            print >>fo, 'MISS\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tlen1={:,}\tlen2={:,}\thttp://humanparalogy.gs.washington.edu/build37/{}'.format(
                r.chrom, r.chromStart, r.chromEnd, r.otherChrom, r.otherStart, r.otherEnd, r.strand,
                r.chromEnd - r.chromStart, r.otherEnd - r.otherStart,
                name
            )
            missed_bases += -r.chromStart+r.chromEnd
            missed_bases += -r.otherStart+r.otherEnd

    partials = defaultdict(list)
    for name, h in hits.items():
        if len(h) == 0: 
            continue

        A = h[0][1]
        oqcov = np.zeros(A[2] - A[1])
        orcov = np.zeros(A[5] - A[4])

        for (((p1, n1, t1, e1), (p2, n2, t2, e2)), A, B) in sorted(h):
            over_qs = max(B[1], A[1])
            over_qe = min(B[2], A[2])
            over_rs = max(B[4], A[4])
            over_re = min(B[5], A[5])
            if over_qs <= over_qe and over_rs <= over_re:
                oqcov[over_qs - A[1]:over_qe - A[1]] = 1
                orcov[over_rs - A[4]:over_re - A[4]] = 1

        p1 = np.sum(oqcov)/len(oqcov)
        p2 = np.sum(orcov)/len(orcov)

        if round(p1*100,0) < 80 or round(p2*100,0) < 80:
            partials[name] += [(p1, p2), h]

        part_missed_bases += len(oqcov) - np.sum(oqcov)
        part_missed_bases += len(orcov) - np.sum(orcov)
        # ok = any(v[0][0][0] >= 100 and v[0][1][0] >= 100 for v in h)
        # if not ok: 

    tm = len(partials)
    print 'Partial: {:7,} hits ({:4.1f}%)'.format(tm, 100.0*tm/len(hits))
    for k in sorted(partials.keys(), key=lambda y: sum(partials[y][0])):
        p1, p2 = partials[k][0]
        print >>fo, 'PART\t{:.2f}\t{:.2f}\thttp://humanparalogy.gs.washington.edu/build37/{}'.format(
            p1*100,p2*100,k)
    tm = sum(1 for k, v in hits.iteritems() if k not in partials and len(v) > 0)
    print 'Full:    {:7,} hits ({:4.1f}%)'.format(tm, 100.0*tm/len(hits))

    fo.close()
except IOError:
    pass
