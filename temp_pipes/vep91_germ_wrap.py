#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['germ_pipe'], \
           config_data['tools']['germ_slurm_wrap']


def queue_vep91(config_file, sample_list, skip):
    fh = open(sample_list, 'r')
    (cores, mem, vep91_pipe, vep91_slurm_wrap) = parse_config(config_file)
    for line in fh:
        sample = line.rstrip('\n')
        job_name = 'dnaseq-vep91-germ-' + sample
        job_log = sample + '.vep91-germ.log'
        batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + ' -o ' + job_log \
                + ' --export=germ_pipe="' + vep91_pipe + '",normal="' + sample_list + '",j="' + config_file + '",f="' \
                + skip + '" ' + vep91_slurm_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + sample + ' failed! Check logs!\n')

    fh.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-s', '--sample-list', action='store', dest='sample_list', help='Normal sample ID list')
    parser.add_argument('-f', '--skip', action='store', dest='skip', help='\'y\' to skip pass filter if already run')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_list, skip) = (inputs.config_file, inputs.sample_list, inputs.skip)
    queue_vep91(config_file, sample_list, skip)