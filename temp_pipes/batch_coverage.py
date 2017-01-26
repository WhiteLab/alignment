#!/usr/bin/python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
import re
import subprocess
from utility.job_manager import job_manager


def get_bam_name(bnid, src_cmd, cont):
    list_cmd = src_cmd + 'swift list ' + cont + ' --prefix ALIGN/' + bnid + '/BAM/'
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    bam = ''
    bai = ''
    dl_cmd = ''
    for fn in re.findall('(.*)\n', flist):
        # depending on software used to index, .bai extension may follow bam
        # test = re.match('^\S+_\w*\d+\.rmdup.srt.ba[m|i]$', fn) or re.match('^\S+_\w*\d+\.rmdup.srt.bam.bai$', fn)
        test = re.match(bnid + '.merged.final.ba[m|i]$', fn) or re.match(bnid + '.merged.final.bam.bai$', fn)
        if test:
            dl_cmd += src_cmd + 'swift download ' + cont + ' ' + fn + ';'
            if fn[-3:] == 'bam':
                bam = fn
            else:
                bai = fn
    return dl_cmd, bam, bai


def calc_coverage(bedtools2_tool, sample, bedfile, cont):
    src_cmd = '. /home/ubuntu/.novarc;'
    job_list = []
    for bnid in open(sample):
        bnid = bnid.rstrip('\n')
        (dl_cmd, bam, bai) = get_bam_name(bnid, src_cmd, cont)
        bed_cmd = bedtools2_tool + ' coverage -hist -abam ' + bam + ' -b ' + bedfile + ' > ' + bnid + '.hist;'
        cleanup = 'rm ' + bam + ' ' + bai + ';'
        final = dl_cmd + bed_cmd + cleanup
        job_list.append(final)
    job_manager(job_list, '8')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Bedtools coverage calculation module.  Typically run last in pipeline.  See coverage parameter.')
    parser.add_argument('-c', '--container', action='store', dest='cont', help='Swift container location')
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
    cont = inputs.cont

    calc_coverage(bedtools2_tool, sample, bedfile, cont)
