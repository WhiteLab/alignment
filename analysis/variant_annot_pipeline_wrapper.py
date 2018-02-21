#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['variant_pipe'], \
           config_data['tools']['variant_slurm_wrap']


def variant_pipe_wrap(config_file, sample_pairs, estep):
    # create sample list
    fh = open(sample_pairs, 'r')
    (cores, mem, variant_pipe, variant_slurm_wrap) = parse_config(config_file)
    for line in fh:
        (sample_pair, tumor_id, normal_id) = line.rstrip('\n').split('\t')

        # quick check to see if just need to restart pipleine from mutect, or actually get merged bams
        job_name = 'dnaseq-annot-' + estep + '-' + sample_pair
        job_log = sample_pair + '.anno.log'
        batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + 'G -o ' + job_log \
                + ' --export=pipeline="' + variant_pipe + '",tumor="' + tumor_id + '",normal="' + normal_id \
                + '",j="' + config_file + '",e="' + estep + '" ' + variant_slurm_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + sample_pair + ' failed! Check logs!\n')

    fh.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations, USE FULL PATH!')
    parser.add_argument('-e', '--execute', action='store', dest='estep',
                        help='Steps to start at, valid entries are start, indel, snv-annot, snv-indel, germ-call, '
                             'germ-annot')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, estep) = (inputs.sample_pairs, inputs.config_file, inputs.estep)
    variant_pipe_wrap(config_file, sample_pairs, estep)
