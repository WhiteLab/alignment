#!/usr/bin/env python3
"""
Usage: ./concat_vcfs2.py <fai> <list>

Arguments:
<fai> fasta index
<list>  vcfs list

Options:
-h
"""
from docopt import docopt

args = docopt(__doc__)
chrs = open(args['<fai>'], 'r')
header_flag = 0
header = []
vcf_entries = {}
order = []
import re
import sys

for entry in chrs:
    entry = entry.rstrip('\n')
    chrom = entry.split('\t')
    vcf_entries[chrom[0]] = []
    order.append(chrom[0])

vcf_list = open(args['<list>'], 'r')
for line in vcf_list:
    line = line.strip()
    m = re.search('(chr[A-Z,a-z,0-9,_]+)', line)
    cur = m.group(1)
    vcf = open(line, 'r')
    for variant in vcf:
        variant = variant.rstrip('\n')
        if variant.startswith('#') and header_flag == 0:
            header.append(variant)
        elif not variant.startswith('#'):
            vcf_entries[cur].append(variant)
    header_flag = 1
    vcf.close()
vcf_list.close()
print('\n'.join(header))
for chrom in order:
    sys.stdout.write('\n'.join(vcf_entries[chrom]))

