#!/usr/bin/env python

import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from job_manager import job_manager
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['metalfox'],config_data['refs']['cont'], config_data['refs']['obj'],
            config_data['refs']['mapability'], config_data['params']['threads'], config_data['params']['ram'])


def metalfox_pipe(config_file, sample_pairs, ref_mnt):
    (metalfox_tool, cont, obj, map_ref, max_t, ram) = parse_config(config_file)
    map_ref = ref_mnt + '/' + map_ref
    src_cmd = '. ~/.novarc;'
    deproxy = 'unset http_proxy; unset https_proxy;'
    pairs = open(sample_pairs, 'r')
    job_list = []
    for sn in pairs:
        sn = sn.rstrip('\n')
        info = sn.split('\t')
        # stdout=PIPE
        sys.stderr.write('Getting bam file name for ' + info[1] + '\n')
        get_bam_name = 'swift list ' + cont + ' --prefix ' + obj + '/' + info[1] + '/BAM/' + info[1]\
                       + ' | grep .rmdup.srt.ba* '
        bam =subprocess.check_output(get_bam_name, shell=True).split('\n')
        dl_bam = 'swift download ' + cont + ' ' + bam[1] + ';swift download ' + cont + ' ' + bam[0] + ';'
        mut_out = 'ANALYSIS/' + info[0] + '/OUTPUT/' + info[0] + '.out.keep'
        dl_out = 'swift download ' + cont + ' ' + mut_out + ';'
        run_metal = metalfox_tool + ' -f1 ' + mut_out + ' -f3 ' + bam[1] + ' -m ' + map_ref + ' > ' + info[0] + \
                    '.foxog_scored_added.out;'
        cleanup = 'rm ' + ' '.join((bam[0],bam[1],mut_out)) + ';'
        job_list.append(deproxy + dl_bam + dl_out + run_metal)# + cleanup)
    pairs.close()
    sys.stderr.write('Queing jobs\n')
    job_manager(job_list, max_t)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Mini pipeline for adding FoxoG scores to calls for filtering.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt) = (inputs.config_file, inputs.sample_pairs, inputs.ref_mnt)
    metalfox_pipe(config_file, sample_pairs, ref_mnt)