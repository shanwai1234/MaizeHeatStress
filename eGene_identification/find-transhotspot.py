import sys
import operator 
from bx.intervals.intersection import Intersecter, Interval
import numpy as np

fh = open(sys.argv[1],'r') # trans-eQTL results from matrixeQTL

mdict = {}
for line in fh:
    new = line.strip().split('\t')
    hap = 'S'+new[0]+':'+new[1]+'-'+new[2]
    nlist = []
    if hap not in mdict:
        mdict[hap] = {}
    for x in new[-1].split(';'):
        if x not in mdict[hap]:
            mdict[hap][x] = 0
        mdict[hap][x] += 1

ndict = {}
for h in mdict:
    if len(mdict[h]) < 3:continue # filtering out less than 3 SNPs inside the haplotype block
    ndict[h] = len(mdict[h])

sorted_x = sorted(ndict.items(), key=operator.itemgetter(1), reverse=True) # sorting haplotype based on its containing SNPs

sh = open(sys.argv[2],'r') # TF information per gene
tdict = {}

for line in sh:
    new = line.strip().split('\t')
    tdict[new[0]] = new[1]

sh.close()

th = open(sys.argv[3],'r') # gene annotation file in the bed format
th.readline()

intersect_dict = {}

for line in th:
    new = line.strip().split('\t')
    chrom = new[0]
    st = int(new[1]) - 2000
    if st < 0:
        start = 0
    else:
        start = st
    stop = int(new[2]) + 2000
    value = new[3]
    if chrom not in intersect_dict:
        intersect_dict[chrom] = Intersecter()
    intersect_dict[chrom].add_interval(Interval(start, stop, value=value))
th.close()

for x in sorted_x:
    n = x[0].split(':')
    chrom = n[0].replace('S','')
    m = n[1].split('-')
    start = m[0]
    stop = m[1]
    a = intersect_dict[chrom].find(int(start), int(stop))
    glist = []
    glist1 = []
    if len(a) > 0:
        for b in a:
            glist.append(b.value)
            if b.value not in tdict:continue
            glist1.append(tdict[b.value])
    else:
        glist.append('No Known Gene')
        glist1.append('No Known Gene')
    nlist = []
    for y in mdict[x[0]]:
        if y in tdict:
            nlist.append(tdict[y])
        else:
            nlist.append(y)
    if len(nlist) == 0:continue
    if len(glist1) == 0:
        glist1 = ['NA']
    else:
        glist1 = glist1
    if len(nlist) < 10:continue
    print (chrom+'\t'+start+'\t'+stop+'\t'+';'.join(glist)+'\t'+';'.join(glist1)+'\t'+';'.join(nlist))

fh.close()
