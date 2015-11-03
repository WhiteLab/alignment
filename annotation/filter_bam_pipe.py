#!/usr/bin/env python
__author__ = 'Miguel'
import json
import sys

sys.path.append('/home/ubuntu/TOOLS/utility')
from job_manager import job_manager

def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['threads'], config_data['refs']['cont'], config_data['refs']['obj'], config_data['tools']['mouse_filter']


def filter_bam_pipe(config_file, lane, ref_mnt):
    (th, cont, obj, mouse_filter) = parse_config(config_file)
    job_list = []
    src_cmd = ". /home/ubuntu/.novarc;"

    fh = open(lane,'r')
    for la in fh:
        la = la.rstrip('\n')
        info = la.split('\t')
        lanes = info[2].split(', ')
        for rg in lanes:
            fn = info[0] + '_' + rg + '.bam'
            stub = info[0] + '_' + rg
            swift_cmd = src_cmd + "swift download " + cont + " --skip-identical " + fn
            mf = mouse_filter + ' -b ' + fn + ' -o ' + stub
            cmd = swift_cmd + '; ' + mf + ';'
            job_list.append(cmd)
    fh.close()
    job_manager(cmd,th)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Filters unwanted cross-species hits and creates a new fastq file.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-l', '--lane', action='store', dest='lane', help='Lane list used for alignment')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Reference drive/folder mount location.  Example would be /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, lane, ref_mnt) = (inputs.config_file, inputs.lane, inputs.ref_mnt)
    filter_bam_pipe(config_file, lane, ref_mnt)
