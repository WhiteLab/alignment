#!/usr/bin/env python
'''
Usage: ./fastq_file <fastq>

Options:
-h

Arguments:
<fastq> fastq file

'''
import sys
import os
import gzip
from Bio import SeqIO
sys.path.append('/lustre/beagle2/mbrown/PDX_TEST/alignment')
sys.path.append('/lustre/beagle2/mbrown/PDX_TEST/utility')

from docopt import docopt
args = docopt(__doc__)
fn = args['<fastq>']
fastq = gzip.open(fn,'rb')
indir = os.path.dirname(fn)
if len(indir) > 0:
    indir = indir + '/'
bn = os.path.basename(fn)
out = gzip.open(indir + 'converted/' + bn,'wb')

# code snippets obtained from https://github.com/vpiro/readtools/blob/master/PHRED_converter.py
quali='fastq-illumina'
qualo='fastq-sanger'
SeqIO.convert(fastq, quali, out, qualo)
fastq.close()
out.close()
