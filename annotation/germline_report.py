#!/usr/bin/env python
import re
from pysam import VariantFile
import sys
import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time


def gen_report(vcf):
    vcf_in = VariantFile(vcf)
    for record in vcf_in.fetch():
        sys.stdout.write('CHROM\tPOS\tREF\tSNP\tTR\tTC\tSYMBOL\tCADD\n')
        pdb.set_trace()
        horse = 'nay'


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='VEP annotated variant file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    vcf = inputs.infile
gen_report(vcf)
