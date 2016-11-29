#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from subprocess import call
from utility.log import log


def filter_indel(bedtools, bed_ref, vcf):
    cmd = bedtools + ' intersect -header -v -a ' + vcf + ' -b ' + bed_ref + ' > '

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Helper script for filtering indel calls against a bed file with dust mask coords')
    parser.add_argument('-b', '--bedtools', action='store', dest='bedtools',
                        help='Location of bedtools')
    parser.add_argument('-r', '--bed_ref', action='store', dest='bed_ref',
                        help='Bed file for filtering')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='vcf file to filter')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (bedtools, bed_ref, vcf) = (inputs.bedtools, inputs.ref, inputs.vcf)
    filter_indel(bedtools, bed_ref, vcf)
