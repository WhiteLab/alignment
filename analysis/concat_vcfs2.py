#!/usr/bin/python
"""
Usage: ./concat_vcfs2.py <fai> <dir>

Arguments:
<fai> fasta index
<dir>  directory vcfs

Options:
-h
"""
from docopt import docopt

args = docopt(__doc__)
chrs = open(args['<fai>'], 'r')
header_flag = 0
for entry in chrs:
    entry = entry.rstrip('\n')
    chrom = entry.split('\t')
    vcf = open(args['<dir>'] + '/' + chrom[0] + '.vcf', 'r')
    for line in vcf:
        line = line.strip()
        if line.startswith('#') and header_flag == 0:
            print line
        elif not line.startswith('#'):
            print line
    header_flag = 1
