#!/usr/bin/python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from subprocess import Popen
from subprocess import call
from utility.job_manager import job_manager


def calc_coverage(bedtools2_tool, sample, bedfile):


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Bedtools coverage calculation module.  Typically run last in pipeline.  See coverage parameter.')
    parser.add_argument('-bt', '--bedtools', action='store', dest='bedtools2_tool', help='Location of bedtools2 tool.')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample list')
    parser.add_argument('-bf', '--bed_file', action='store', dest='bed_file', help='Bed file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    bedtools2_tool = inputs.bedtools2_tool
    sample = inputs.sample
    bedfile = inputs.bed_file

    calc_coverage(bedtools2_tool, sample, bedfile)