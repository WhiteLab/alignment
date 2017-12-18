#!/usr/bin/env python3

import json
import os
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project'], config_data['refs']['project_dir'], config_data['refs']['align']


def list_bam(project, project_dir, align_dir, sample):
    project_dir + project
    bam = project_dir + project + '/' + align_dir + '/' + sample + '/BAM/' + sample + '.merged.final.bam'
    check_file = os.path.isfile(bam)
    if not check_file:
        sys.stderr.write(date_time() + 'Merged bam ' + bam + ' not found.\n')
    return check_file


def check_for_merged_bams(config_file, sample_list):
    fh = open(sample_list, 'r')
    (project, project_dir, align_dir) = parse_config(config_file)
    missing = []
    for sample in fh:
        sample = sample.rstrip('\n')
        check = list_bam(project, project_dir, align_dir, sample)
        if not check:
            missing.append(sample)
    return missing


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample/project prefix list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file) = (inputs.sample_list, inputs.config_file)
    check_for_merged_bams(config_file, sample_list)
