#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call
from utility.log import log


def filter_indel(bedtools, bed_ref, pair, loc):
    cmd = bedtools + ' intersect -header -v -a somatic.indel.vcf -b ' + bed_ref + ' > ' + pair \
          + '.somatic.indel.dustmasked.vcf'
    log(loc, date_time() + cmd + '\n')
    check = call(cmd, shell=True)
    return check


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Helper script for filtering indel calls against a bed file with dust mask coords')
    parser.add_argument('-b', '--bedtools', action='store', dest='bedtools', help='Location of bedtools')
    parser.add_argument('-r', '--bed_ref', action='store', dest='bed_ref', help='Bed file for filtering')
    parser.add_argument('-sp', '--sample_pair', action='store', dest='pair', help='vcf file to filter')
    parser.add_argument('-l', '--log', action='store', dest='loc', help='Log file output location')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (bedtools, bed_ref, pair, loc) = (inputs.bedtools, inputs.ref, inputs.pair, inputs.loc)
    filter_indel(bedtools, bed_ref, pair, loc)
