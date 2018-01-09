#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import json
from subprocess import call


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['germ_pipe'], \
           config_data['tools']['germ_slurm_wrap']


def germ_pipe_wrap(config_file, samples, estep):
    (cores, mem, germ_pipe, germ_slurm_wrap) = parse_config(config_file)
    for sample in open(samples):
        job_name = 'dnaseq-germ-' + sample
        job_log = sample + '.anno.log'
        batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + ' -o ' + job_log \
                + ' --export=germ_pipe="' + germ_pipe + '",normal="' + sample \
                + '",j="' + config_file + '",e="' + estep + '" ' + germ_slurm_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + sample + ' failed! Check logs!\n')
    sys.stderr.write(date_time() + "Jobs submitted.  Check logs for any errors\n")

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mini temp pipe to run a part of variant_annot_pipe')
    parser.add_argument('-s', '--sample-list', action='store', dest='samples',
                        help='Single sample list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-e', '--estep', action='store', dest='estep',
                        help='Step to start at, \'start\' or \'annot\'')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (samples, config_file, estep) = (
        inputs.samples, inputs.config_file, inputs.estep)
    germ_pipe_wrap(config_file, samples, estep)
