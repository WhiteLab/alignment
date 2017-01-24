#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from analysis.variant_annot_pipe import *


def cleanup_bams():
    for bnid in open('temp_list.txt'):
        bnid = bnid.rstrip('\n')
        cmd = 'rm ' + bnid + '.merged.*;'
        call(cmd, shell=True)


def temp_germline_pipe(config_file, samples, ref_mnt):
    mk_dir = 'mkdir BAM LOGS ANALYSIS ANNOTATION'
    call(mk_dir, shell=True)
    sys.stderr.write(date_time() + 'Splitting up samples list\n')
    m = 4
    x = 0
    f = 0
    temp_name = 'temp_list.txt'
    temp = open(temp_name, 'w')

    for sample in open(samples):
        f = 0
        x += 1
        temp.write(sample)

        if x % m == 0:
            f = 1
            temp.close()
            x = 0
            check = get_merged_bams(config_file, temp_name)
            if len(check) != 0:
                sys.stderr.write('Can\'t find merged bams around ' + sample)
            check = platypus_germline(config_file, temp_name, 'LOGS/', 'n', ref_mnt)
            if check != 0:
                sys.stderr.write(date_time() + 'platypus failed around ' + sample + '\n')
            check += annot_platypus(config_file, temp_name, ref_mnt)
            if check == 0:
                sys.stderr.write(date_time() + 'Germ line call complete\n')
            else:
                sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
                exit(1)
            cleanup_bams()
            temp = open(temp_name, 'w')
    temp.close()
    if f == 0:
        # do remaining is total number of pairs left is less than m
        get_merged_bams(config_file, temp_name)
        check = platypus_germline(config_file, temp_name, 'LOGS/', 'n', ref_mnt)
        if check != 0:
            sys.stderr.write(date_time() + 'platypus failed a the end of list\n')
        check += annot_platypus(config_file, temp_name, ref_mnt)
        if check == 0:
            sys.stderr.write(date_time() + 'Germ line call complete\n')
        else:
            sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
            exit(1)
        cleanup_bams()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mini temp pipe to run a part of variant_annot_pipe')
    parser.add_argument('-s', '--sample-list', action='store', dest='samples',
                        help='Single sample list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, ref_mnt) = (
        inputs.samples, inputs.config_file, inputs.ref_mnt)
    temp_germline_pipe(config_file, sample_pairs, ref_mnt)
