#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import json
import os
from utility.date_time import date_time
from subprocess import call
from utility.log import log
from alignment.get_merged_bams import get_merged_bams
from dustmask_filter import filter_indel


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['scalpel'], config_data['tools']['bedtools'], config_data['refs']['capture'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads'], config_data['params']['dustmask_flag'],\
           config_data['refs']['dustmask'], config_data['params']['wg_flag']


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


def wg_mode(scalpel, tumor_bam, normal_bam, bed, fasta, cpus, pair):
    for coords in open(bed):
        cur = coords.rstrip('\n').split('\t')
        c_string = cur[0] + ':' + str(int(cur[1]) + 1) + '-' + str(int(cur[2]) + 1)
        loc = 'LOGS/' + pair + '_' + cur[0] + '.scalpel.log'
        cmd = scalpel + ' --somatic --logs --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' \
        + normal_bam + ' --window 600 --two-pass --bed ' + c_string + ' --ref ' + fasta + ' 2> ' + loc
        log(loc, date_time() + cmd + '\n')
        check = call(cmd, shell=True)
        if check != 0:
            return 1, cur[0], pair
        mv_cmd = 'mkdir ' + cur[0] + '; mv outdir/main/* ' + cur[0] + '; rm -rf outdir/main;'
        call(mv_cmd, shell=True)
    return 0



def scalpel_indel(pairs, log_dir, config_file, ref_mnt):
    (scalpel, bedtools, bed, fasta, cpus, dustmask_flag, dustmask_bed, wg) = parse_config(config_file)
    bed = ref_mnt + '/' + bed
    fasta = ref_mnt + '/' + fasta
    dustmask_bed = ref_mnt + '/' + dustmask_bed
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
        if wg == 'n':
            scalpel_cmd = scalpel + ' --somatic --logs --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' \
                          + normal_bam + ' --bed ' + bed + ' --ref ' + fasta + ' 2>> ' + loc
            sys.stderr.write(date_time() + 'Starting indel calls for ' + cur[0] + '\n')
            log(loc, date_time() + 'Starting indel calls for ' + cur[0] + ' in capture mode with command:\n'
                + scalpel_cmd + '\n')
            check = call(scalpel_cmd, shell=True)
            if check != 0:
                sys.stderr.write(date_time() + 'Indel calling failed for pair ' + cur[0] + ' with command:\n' +
                                 scalpel_cmd + '\n')
                exit(1)
        else:
            check = wg_mode(scalpel, tumor_bam, normal_bam, bed, fasta, cpus, cur[0])
            if check[0] != 0:
                sys.stderr.write('Scalpel failed for ' + cur[2] + ' at ' + cur[1] + '\n')
                exit(1)
        log(loc, date_time() + 'Indel calling complete for pair ' + cur[0] + ' moving output files\n')
        mv_cmd = 'mkdir ' + cur[0] + '; mv outdir/main/* ' + cur[0] + '; rm -rf outdir/main;'
        log(loc, date_time() + mv_cmd + '\n')
        call(mv_cmd, shell=True)
        sys.stderr.write(date_time() + 'Completed indel calls for ' + cur[0] + '\n')
        if dustmask_flag == 'Y':
            log(loc, date_time() + 'Filter dustmask flag given\n')
            check = filter_indel(bedtools, dustmask_bed, cur[0])
            if check != 0:
                sys.stderr.write(date_time() + 'Dustmask failed for ' + cur[0] + '\n')
                exit(1)
            else:
                log(loc, date_time() + 'Dustmask complete for ' + cur[0] + '\n')
    fh.close()
    sys.stderr.write(date_time() + 'Indel call completed\n')
    return 0


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
