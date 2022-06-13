#!/usr/bin/env python

## calculate N50 from fasta file
## N50 = contig length such that half of the contigs are longer and 1/2 of contigs are shorter

#import commands
import sys
import os
from itertools import groupby
import numpy
# from Bio import SeqIO


lengths = []

fasta = sys.stdin
## parse each sequence by header: groupby(data, key)
faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

for record in faiter:
    ## join sequence lines
    seq = "".join(s.strip() for s in next(faiter))
    #print(seq)
    lengths.append(len(seq))

## sort contigs longest>shortest
all_len=sorted(lengths, reverse=True)
csum=numpy.cumsum(all_len)

print("N: %d" % int(sum(lengths)))
n2=int(sum(lengths)/2)

# get index for cumsum >= N/2
csumn2=min(csum[csum >= n2])
ind=numpy.where(csum == csumn2)

n50 = all_len[int(ind[0])]
print("N50: %s" % n50)

if __name__ == "__main__":
    input = sys.stdin
    output = sys.stdout
    N50 = str(n50)
    output.write(f'N50:{N50}, \n')