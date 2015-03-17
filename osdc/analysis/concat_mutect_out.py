#!/usr/bin/python

import sys

if(sys.argv[1] is None):
    print "Usage: base_filename (i.e. 2011-123_2012-345)"
    exit()

base_filename = sys.argv[1]

chrs = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')

header_flag = 0
for chrom in chrs:
    dot_out = open(base_filename + '.' + chrom + '.out.keep', 'r')
    for line in dot_out:
        line = line.strip()
        if(line.startswith('#') and header_flag == 0):
            print line
        if(line.startswith('contig') and header_flag == 0):
            print line
            header_flag = 1
        elif(not line.startswith('#')):
            print line


