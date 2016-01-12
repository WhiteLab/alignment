#!/usr/bin/env python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from job_manager import job_manager
import json
import pdb

def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['picard'], config_data['refs']['fa_ordered'],
            config_data['refs']['intervals'], config_data['refs']['cont'], config_data['refs']['obj'],
            config_data['params']['threads'], config_data['params']['ram'])

def oxog_check(config_file, lane_list, ref_mnt):
    (java, picard, fa_ordered, intervals, cont, obj, max_t, ram) = parse_config(config_file)
    src_cmd = ". /home/ubuntu/.novarc;"
    ram = str(int(ram)/int(max_t))
    job_list=[]
    lane_fh = open(lane_list, 'r')
    for sample in lane_fh:
        sample = sample.rstrip('\n')
        info = sample.split('\t')
        lanes = info[2].split(', ')
        bid = info[0]
        for lane in lanes:
            dl_bam = src_cmd + 'swift download ' + cont + ' --prefix ' + obj + '/' + bid + '/BAM/' + bid \
                 + '_' + lane + '.rmdup.srt.ba;'
            # pdb.set_trace()
            bam = obj + '/' + bid + '/BAM/' + bid + '_' + lane + '.rmdup.srt.bam'
            oxoG = java + ' -Xmx' + ram + 'g -jar ' + picard + ' CollectOxoGMetrics I=' + bam + ' O=' + bid + '_' + \
                    lane + '.oxo_summary.txt R=' + ref_mnt + '/' + fa_ordered + ' INTERVALS=' + ref_mnt + '/' + intervals + ' 2> ' + bid + '_' + lane + '.log;'
            del_bam = 'rm ' + obj + '/' + bid + '/BAM/' + bid + '_' + lane + '.rmdup.srt.bam;'
            job_list.append(dl_bam + oxoG + del_bam)
    lane_fh.close()
    job_manager(job_list, max_t)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Pipleine for checking custom capture samples for oxidative damage.'
                    '  Need BAM and bai files ahead of time.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-f', '--lane_list', action='store', dest='lane_list', help='Lane list used to run alignment - ')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, lane_list, ref_mnt) = (inputs.config_file, inputs.lane_list, inputs.ref_mnt)
    oxog_check(config_file, lane_list, ref_mnt)
