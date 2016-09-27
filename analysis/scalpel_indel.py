#!/usr/bin/env python

import sys
import json
import os
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
from date_time import date_time
from subprocess import call
from log import log
from get_merged_bams import get_merged_bams


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['scalpel'], config_data['refs']['capture'], config_data['refs']['fa_ordered'], \
           config_data['params']['threads']


def create_sample_list(sample_pairs):
    sample_list = 'sample_list.txt'
    fh = open(sample_pairs, 'r')
    sl = open(sample_list, 'w')
    temp = {}
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        if cur[1] not in temp:
            sl.write(cur[1] + '\n')
            temp[cur[1]] = 1
        if cur[2] not in temp:
            sl.write(cur[2] + '\n')
            temp[cur[2]] = 1
    sl.close()
    fh .close()
    del temp


def scalpel_indel(pairs, log_dir, config_file, ref_mnt):
    (scalpel, bed, fasta, cpus) = parse_config(config_file)
    bed = ref_mnt + '/' + bed
    fasta = ref_mnt + '/' + fasta
    # use get_merged_bams api
    sample_list = 'sample_list.txt'
    if not os.path.isfile(sample_list):
        create_sample_list(pairs)
        sys.stderr.write(date_time() + 'Sample pairs list not created - creating one since this is being run likely '
                                       'outside of pipeline')
        get_merged_bams(config_file, sample_list)
    fh = open(pairs, 'r')
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        loc = log_dir + cur[0] + '.scalpel.log'
        tumor_bam = cur[1] + '.merged.final.bam'
        normal_bam = cur[2] + '.merged.final.bam'
        scalpel_cmd = scalpel + ' --somatic --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' + normal_bam \
                      + ' --bed ' + bed + ' --ref ' + fasta + ' 2>> ' + loc
        log(loc, date_time() + 'Starting indel calls for ' + cur[0] + ' with command:\n' + scalpel_cmd + '\n')
        check = call(scalpel_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'Indel calling failed for pair ' + cur[0] + ' with command:\n' +
                             scalpel_cmd + '\n')
        log(loc, date_time() + 'Indel calling complete for pair ' + cur[0])
    fh.close()
    sys.stderr.write(date_time() + 'Indel call completed\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Scalpel indel caller wrapper script.  Samples must have been aligned and bams merged '
                    'ahead of time')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, log_dir, config_file, ref_mnt) = (inputs.sample_pairs, inputs.log_dir, inputs.config_file
                                                     , inputs.ref_mnt)
    scalpel_indel(sample_pairs, log_dir, config_file, ref_mnt)
