#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import re
import subprocess
from utility.job_manager import job_manager


def get_bam_name(bnid, src_cmd, cont, obj):
    list_cmd = src_cmd + 'swift list ' + cont + ' --prefix ' + obj + '/' + bnid + '/BAM/'
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    bam = ''
    bai = ''
    dl_cmd = ''
    bam2 = []
    bai2 = []
    dl_cmd2 = []
    # flag to see if merged found, if not, use rmdup singleton files
    mflag = 0
    dl_ind = 0
    dl_ct = 1
    for fn in re.findall('(.*)\n', flist):
        # depending on software used to index, .bai extension may follow bam
        test = re.search(bnid + '.merged.final.ba[m|i]$', fn) or re.search(bnid + '.merged.final.bam.bai$', fn)
        if test:
            dl_cmd += src_cmd + 'swift download ' + cont + ' ' + fn + ';'
            mflag = 1
            if fn[-3:] == 'bam':
                bam = fn
            else:
                bai = fn
        test = re.match('^\S+_\w*\d+\.rmdup.srt.ba[m|i]$', fn) or re.match('^\S+_\w*\d+\.rmdup.srt.bam.bai$', fn)
        if test:
            if fn[-3:] == 'bam':
                bam = fn
                bam2.append(bam)

            else:
                bai = fn
                bai2.append(bai)
            if not (dl_ct % 2 == 0):
                dl_cmd2.append(src_cmd + 'swift download ' + cont + ' ' + fn + ';')
            else:
                dl_cmd2[dl_ind] += src_cmd + 'swift download ' + cont + ' ' + fn + ';'
                dl_ind += 1
            dl_ct += 1
    if mflag == 1:
        return dl_cmd, bam, bai
    else:
        return dl_cmd2, bam2, bai2


def calc_coverage(bedtools2_tool, sample, bedfile, cont, obj):
    src_cmd = '. /home/ubuntu/.novarc;'
    job_list = []
    for bnid in open(sample):
        bnid = bnid.rstrip('\n')
        (dl_cmd, bam, bai) = get_bam_name(bnid, src_cmd, cont, obj)
        if isinstance(bam, str):
            bed_cmd = bedtools2_tool + ' coverage -hist -abam ' + bam + ' -b ' + bedfile + ' > ' + bnid + '.hist;'
            cleanup = 'rm ' + bam + ' ' + bai + ';'
            final = dl_cmd + bed_cmd + cleanup
            job_list.append(final)
        else:
            for i in range(len(bam)):
                bed_cmd = bedtools2_tool + ' coverage -hist -abam ' + bam[i] + ' -b ' + bedfile + ' > ' + bnid \
                          + '_' + str(i) + '.hist;'
                cleanup = 'rm ' + bam[i] + ' ' + bai[i] + ';'
                final = dl_cmd[i] + bed_cmd + cleanup
                job_list.append(final)
    job_manager(job_list, '8')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Bedtools coverage calculation module.  Typically run last in pipeline.  See coverage parameter.')
    parser.add_argument('-c', '--container', action='store', dest='cont', help='Swift container location')
    parser.add_argument('-o', '--object', action='store', dest='obj', help='Swift object prefix')
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
    obj = inputs.obj

    calc_coverage(bedtools2_tool, sample, bedfile, cont, obj)
