#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['strelka_pipe'], \
           config_data['tools']['strelka_slurm_wrap']


def queue_strelka(config_file, sample_list):
    fh = open(sample_list, 'r')
    (cores, mem, strelka_pipe, strelka_slurm_wrap) = parse_config(config_file)
    for line in fh:
        (pair, tumor, normal) = line.rstrip('\n').split('\t')
        job_name = 'dnaseq-vep91-strelka-' + pair
        job_log = pair + '.vep91-strelka.log'
        batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + 'G -o ' + job_log \
                + ' --export=strelka_pipe="' + strelka_pipe + '",tumor="' + tumor + '",normal="' + normal + '",j="' \
                + config_file + '" ' + strelka_slurm_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + pair + ' failed! Check logs!\n')

    fh.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-l', '--lane-list', action='store', dest='lane_list', help='Lane list in '
                                        'format <sample_id_seprated_by_underscore><tab><tumor id> <tab><normal id>')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, lane_list) = (inputs.config_file, inputs.lane_list)
    queue_strelka(config_file, lane_list)
