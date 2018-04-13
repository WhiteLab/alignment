#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import subprocess
import os
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'], \
           config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['cnv_pipe'],\
           config_data['tools']['cnv_slurm']


def get_bam_name(bnid, project_dir, project, align_dir):
    bam_dir = project_dir + project + '/' + align_dir + '/' + bnid + '/BAM/'
    bam = bam_dir + bnid + '.merged.final.bam'
    bai = bam_dir + bnid + '.merged.final.bai'
    f = 0
    if not os.path.isfile(bam):
        sys.stderr.write(date_time() + 'Bam not found in ' + bam_dir + '\n')
        f = 1
        return f, bam, bai
    if not os.path.isfile(bai):
        bai = bam_dir + bnid + '.merged.final.bam.bai'
        if not os.path.isfile(bai):
            sys.stderr.write(date_time() + 'Bam index file for ' + bnid + ' not found!  Please index first\n')
            f = 1
            return f, bam, bai

    return f, bam, bai


def cnv_wrap(config_file, sample_pairs, project2, o_flag):
    (project_dir, project, align_dir, cores, mem, cnv_pipe, cnv_slurm) = parse_config(config_file)

    pair_list = []

    # calc coverage for all gene capture regions
    for pairs in open(sample_pairs):
        pair_set = pairs.rstrip('\n').split('\t')
        pair_list.append('\t'.join((pair_set[1], pair_set[2])))
        (flag, tum_bam, tum_bai) = get_bam_name(pair_set[1], project_dir, project, align_dir)
        # ensure bam is complete, otherwise try a different location
        if not (flag == 0 and os.path.getsize(tum_bam) > 0):
            sys.stderr.write(date_time() + 'Locating bam failed for ' + pair_set[1] + ', trying backup project\n')
            (flag, tum_bam, tum_bai) = get_bam_name(pair_set[1], project_dir, project2, align_dir)
            if not (flag == 0 and os.path.getsize(tum_bam) > 0):
                sys.stderr.write(date_time() + 'Suitable bam for ' + pair_set[1] + ' not found. Check parameters!\n')
        (flag, norm_bam, norm_bai) = get_bam_name(pair_set[2], project_dir, project, align_dir)
        # ensure bam is complete, otherwise try a different location
        if not (flag == 0 and os.path.getsize(norm_bam) > 0):
            sys.stderr.write(date_time() + 'Locating bam failed for ' + pair_set[2] + ', trying backup container\n')
            (flag, norm_bam, norm_bai) = get_bam_name(pair_set[2], project_dir, project2, align_dir)
            if not (flag == 0 and os.path.getsize(norm_bam) > 0):
                sys.stderr.write(date_time() + 'Suitable bam for ' + pair_set[2] + ' not found. Check parameters!\n')
        if flag == 0:
            job_name = 'cnv-est_' + pair_set[0]
            job_log = pair_set[0] + '.cnv.log'
            batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + 'G -o ' + job_log \
                    + ' --export=cnv_pipe="' + cnv_pipe + '",tumor="' + tum_bam + '",normal="' + norm_bam \
                    + '",j="' + config_file + '",o="' + o_flag + '",p="' + project2 + '" ' + cnv_slurm
            sys.stderr.write(date_time() + 'Submitted cnv est job for ' + pair_set[0] + '\n' + batch + '\n')
            subprocess.call(batch, shell=True)
        else:
            sys.stderr.write(date_time() + 'Skipping pair ' + pair_set[0]
                             + ', check params and file locations!\n')

    sys.stderr.write(date_time() + 'Submitted all cnv jobs!\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline wrapper to estimate CNV from tumor/normal '
                                                 'bams from custom capture')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-o', '--project2', action='store', dest='project2',
                        help='Backup project in case the first location does not work')
    parser.add_argument('-w', '--overwrite', action='store', dest='o_flag',
                        help='overwrite flag to determine whether to wrote over existing coverage files')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, project2, o_flag) = (inputs.sample_pairs, inputs.config_file, inputs.project2,
                                                     inputs.o_flag)
    cnv_wrap(config_file, sample_pairs, project2, o_flag)
