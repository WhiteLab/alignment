#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from alignment.novosort_merge_pe import novosort_merge_pe
from alignment.check_for_merged_bams import check_for_merged_bams


def run_novosort(config_file, sample_list):
        check = novosort_merge_pe(config_file, sample_list)
        if check == 0:
            sys.stderr.write(date_time() + 'File merge complete!\n')

        else:
            sys.stderr.write(date_time() + 'File download and merge failed.\n')
            exit(1)


def preprocess_bams(config_file, sample_pairs):
    # create sample list
    sample_list = 'sample_list.txt'
    fh = open(sample_pairs, 'r')
    sl = open(sample_list, 'w')
    temp = {}
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        if len(cur) == 3:
            if cur[1] not in temp:
                sl.write(cur[1] + '\n')
                temp[cur[1]] = 1
            if cur[2] not in temp:
                sl.write(cur[2] + '\n')
                temp[cur[2]] = 1
        else:
            if cur[0] not in temp:
                sl.write(cur[1] + '\n')
                temp[cur[1]] = 1
    sl.close()
    fh .close()
    miss_list = check_for_merged_bams(config_file, sample_list)
    if len(miss_list) > 0:
        sys.stderr.write(date_time() + 'Missing files detected, merging lane files\n')
        temp_fn = 'temp_samp_list.txt'
        temp_fh = open(temp_fn, 'w')
        temp_fh.write('\n'.join(miss_list))
        temp_fh.close()
        run_novosort(config_file, temp_fn)
    else:
        sys.stderr.write(date_time() + 'All bams found. Ready for next step!\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pre-variant calling step to merge all bam files for each sample '
                                                 'before running')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list or single ID list if not paired')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file) = (inputs.sample_pairs, inputs.config_file)
    preprocess_bams(config_file, sample_pairs)
