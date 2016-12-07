#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from analysis.variant_annot_pipe import *


def cleanup_bams():
    for bnid in open('sample_list.txt'):
        bnid = bnid.rstrip('\n')
        cmd = 'rm ' + bnid + '.merged.*;'
        call(cmd, shell=True)
    cmd = 'rm sample_list.txt'
    call(cmd, shell=True)


def temp_indel_pipe(config_file, sample_pairs, ref_mnt):
    mk_dir = 'mkdir BAM LOGS ANALYSIS ANNOTATION'
    call(mk_dir, shell=True)
    sys.stderr.write(date_time() + 'Splitting up sample pairs list\n')
    m = 4
    x = 0
    f = 0
    temp = open('temp_pairs.txt', 'w')
    for pair in open(sample_pairs):
        f = 0
        x += 1
        temp.write(pair)
        cmd = 'mkdir ' + pair
        call(cmd, shell=True)
        if x % m == 0:
            f = 1
            temp.close()
            x = 0
            check = scalpel_indel(temp, 'LOGS/', config_file, ref_mnt)
            if check != 0:
                sys.stderr.write(date_time() + 'scalpel failed around ' + pair + '\n')
            vep(config_file, temp, ref_mnt, '.indel.vcf', '.somatic.indel.vep.vcf', 'scalpel')
            # leave only vcfs and reports for manual upload, clear out merged bams and sample_list, scalpel creates
            # one each time
            cleanup_bams()
            temp = open('temp_pairs.txt', 'w')
    temp.close()
    if f == 0:
        # do remaining is total number of pairs left is less than m
        check = scalpel_indel(temp, 'LOGS/', config_file, ref_mnt)
        if check != 0:
            sys.stderr.write(date_time() + 'scalpel failed at last set \n')
        vep(config_file, temp, ref_mnt, '.indel.vcf', '.somatic.indel.vep.vcf', 'scalpel')
        # leave only vcfs and reports for manual upload, clear out merged bams and sample_list, scalpel creates
        # one each time
        cleanup_bams()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mini temp pipe to run a part of variant_annot_pipe')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, ref_mnt) = (
        inputs.sample_pairs, inputs.config_file, inputs.ref_mnt)
    temp_indel_pipe(config_file, sample_pairs, ref_mnt)
