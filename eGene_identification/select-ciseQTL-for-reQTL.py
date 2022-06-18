import sys
import operator
from scipy.stats import pearsonr
import numpy as np

fh = open(sys.argv[1], 'r')
# SNP     gene    beta    t-stat  p-value FDR
# S9_126066333    Zm00001d047300  -0.882832011410394      -25.6020892583907       2.8578889610346e-37     2.95687861360795e-30
# S9_126066333    Zm00001d047300  -0.882832011410394      -25.6020892583907       2.8578889610346e-37     2.18487345904983e-27    cis
cdict = {} # control condition
for line in fh:
    new = line.strip().split('\t')
    if float(new[-2]) >= 0.01:
        continue
    if new[1] not in cdict:
        cdict[new[1]] = {}
    if new[0] not in cdict[new[1]]:
        cdict[new[1]][new[0]] = float(new[-2])

fh.close()

sh = open(sys.argv[2], 'r')

hdict = {} # heat condition
for line in sh:
    new = line.strip().split('\t')
    if float(new[-2]) >= 0.01:
        continue
    if new[1] not in hdict:
        hdict[new[1]] = {}
    if new[0] not in hdict[new[1]]:
        hdict[new[1]][new[0]] = float(new[-2])

sh.close()

# cdict1 and hdict1 are going to pull top significants SNPs per gene
cdict1 = {}
cgene = set([])
for c in cdict:
    sorted_c = sorted(cdict[c].items(), key=operator.itemgetter(1))
    top = sorted_c[0][1]
    gene = c
    cgene.add(gene)
    for x in sorted_c:
        if x[1] == top:
            snp = x[0]
            if gene not in cdict1:
                cdict1[gene] = {}
            if snp not in cdict1[gene]:
                cdict1[gene][snp] = top

hdict1 = {}
hgene = set([])
for h in hdict:
    sorted_h = sorted(hdict[h].items(), key=operator.itemgetter(1))
    top = sorted_h[0][1]
    gene = h
    hgene.add(gene)
    for x in sorted_h:
        if x[1] == top:
            snp = x[0]
            if gene not in hdict1:
                hdict1[gene] = {}
            if snp not in hdict1[gene]:
                hdict1[gene][snp] = top

shared = cgene.intersection(hgene)
cuniq = cgene - hgene
huniq = hgene - cgene

snp = open(sys.argv[3], 'r') # the genotype file
snp.readline()

sdict = {}
for line in snp:
    new = line.strip().split('\t')
    if new[0] not in sdict:
        sdict[new[0]] = list(map(int, new[1:]))
snp.close()

for x in sorted(list(shared)):
    csnp = set(sorted(list(cdict1[x].keys())))
    hsnp = set(sorted(list(hdict1[x].keys())))
    snpshare = csnp.intersection(hsnp)
    if len(snpshare) != 0:
        fsnp = list(snpshare)[0]
        n = fsnp.split('_')
        chrom = 'Chr' + n[0].replace('S', '')
        pos = n[1]
        print (x + '\t' + chrom + '\t' + pos + '\t' + fsnp)
    else:
        corr, _ = pearsonr(sdict[list(csnp)[0]], sdict[list(hsnp)[0]])
        n1 = list(csnp)[0].split('_')
        chrom1 = 'Chr' + n1[0].replace('S', '')
        pos1 = n1[1]
        n2 = list(hsnp)[0].split('_')
        chrom2 = 'Chr' + n2[0].replace('S', '')
        pos2 = n2[1]
        if corr**2 > 0.8:
            print (x + '\t' + chrom1 + '\t' + pos1 + '\t' + list(csnp)[0])
        else:
            print (x + '\t' + chrom1 + '\t' + pos1 + '\t' + list(csnp)[0])
            print (x + '\t' + chrom2 + '\t' + pos2 + '\t' + list(hsnp)[0])

def uniq_gene(gset, odict):
    for x in list(gset):
        snp = set(list(odict[x].keys()))
        n = list(snp)[0].split('_')
        chrom = 'Chr' + n[0].replace('S', '')
        pos = n[1]
        print (x + '\t' + chrom + '\t' + pos + '\t' + list(snp)[0])

uniq_gene(cuniq, cdict1)
uniq_gene(huniq, hdict1)
