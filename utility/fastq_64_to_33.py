#!/usr/bin/env python3
"""
Usage: ./fastq_file <fastq>

Options:
-h

Arguments:
<fastq> fastq file

"""
import gzip
import os
import sys

from Bio import SeqIO

sys.path.append('/lustre/beagle2/mbrown/PDX_TEST/alignment')
sys.path.append('/lustre/beagle2/mbrown/PDX_TEST/utility')

from docopt import docopt

args = docopt(__doc__)
fn = args['<fastq>']
fastq = gzip.open(fn, 'rb')
indir = os.path.dirname(fn)
bn = os.path.basename(fn)
if not os.path.isdir(indir + '/converted'):
    os.mkdir(indir + '/converted')
out = gzip.open(indir + '/converted/' + bn, 'wb')

# code snippets obtained from https://github.com/vpiro/readtools/blob/master/PHRED_converter.py
quali = 'fastq-illumina'
qualo = 'fastq-sanger'
SeqIO.convert(fastq, quali, out, qualo)
fastq.close()
out.close()
