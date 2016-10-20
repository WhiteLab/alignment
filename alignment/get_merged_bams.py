#!/usr/bin/env python
import json
import os
import re
import sys
from utility.job_manager import job_manager
from utility.date_time import date_time
import subprocess


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['refs']['cont'], config_data['refs']['obj'], config_data['params']['threads'],
            config_data['params']['ram'])


def list_bam(cont, obj, sample, threads):
    list_cmd = '. /home/ubuntu/.novarc;swift list ' + cont + ' --prefix ' + obj + '/' + sample
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    p = []

    for fn in re.findall('(.*)\n', flist):
        if re.match('.*.merged.final.ba', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + fn + ' --output '\
                     + os.path.basename(fn)
            p.append(dl_cmd)
    f = job_manager(p, threads)
    if f == 0:
        sys.stderr.write(date_time() + 'BAM download complete\n')
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)


def get_merged_bams(config_file, sample_list):
    fh = open(sample_list, 'r')
    (cont, obj, threads, ram) = parse_config(config_file)
    for sample in fh:
        sample = sample.rstrip('\n')
        list_bam(cont, obj, sample, threads)
    return 0


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
    get_merged_bams(config_file, sample_list)
