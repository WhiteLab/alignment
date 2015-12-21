#!/usr/bin/env python
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

for entry in chrs:
    entry = entry.rstrip('\n')
    chrom = entry.split('\t')
    vcf_entries[chrom[0]] = []
    order.append(chrom[0])

vcf = open(args['<list>'], 'r')
for line in vcf:
    line = line.strip()
    m = re.search('(chr\S+)\.', line)
    cur = m.group(1)
    for variant in line:
        if variant.startswith('#') and header_flag == 0:
            header.append(variant)
        elif not variant.startswith('#'):
            vcf_entries[cur].append(variant)
    header_flag = 1
vcf.close()
print header
for chrom in order:
    print vcf_entries[chrom]

