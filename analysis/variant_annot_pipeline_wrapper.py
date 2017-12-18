#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['params']['ram'], config_data['refs']['variant_pipe'], \
           config_data['refs']['variant_slurm_wrap']


def variant_pipe_wrap(config_file, sample_pairs):
    # create eventual output location directories

    # create sample list
    fh = open(sample_pairs, 'r')

    for line in fh:
        (sample_pair, tumor_id, normal_id) = line.rstrip('\n').split('\t')
        (cores, mem, variant_pipe, variant_slurm_wrap) = parse_config(config_file)
        # likely not needed, will revisit
        # if kflag == 'y':
        #     # create bam list for ksort
        #     bam_list = 'bam_list.txt'
        #     blist_cmd = 'ls *.merged.final.bam > ' + bam_list
        #     call(blist_cmd, shell=True)
        #     check = ksort(config_file, bam_list, kflag, ref_mnt)
        #     if check == 0:
        #         sys.stderr.write(date_time() + 'Karyotypic reorder of BAM files completed\n')
        #     else:
        #         sys.stderr.write(date_time() + 'Karyotypic reorder of BAM files failed.\n')
        #         exit(1)
        # quick check to see if just need to restart pipleine from mutect, or actually get merged bams
        job_log = sample_pair + '.anno.log'
        batch = 'sbatch -c ' + cores + ' --mem ' + mem + ' -o ' + job_log \
                + ' --export=pipeline="' + variant_pipe + '",tumor="' + tumor_id + '",normal="' + normal_id \
                + '",j="' + config_file + '"' + ' ' + variant_slurm_wrap
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
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file) = (inputs.sample_pairs, inputs.config_file)
    variant_pipe_wrap(config_file, sample_pairs)
