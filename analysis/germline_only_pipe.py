#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from analysis.variant_annot_pipe import *
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['analysis'], \
           config_data['refs']['annotation'], config_data['params']['user'], config_data['params']['group']


def temp_germline_pipe(config_file, sample, estep):
    (project_dir, project , analysis, annotation, user, group) = parse_config(config_file)
    src_env = '. /etc/environment'
    call(src_env, shell=True)
    if estep == 'start':
        germ_ana_dir = project_dir + project + '/' + analysis + '/' + sample
        if not os.path.isdir(germ_ana_dir):
            # might as well make log directory too
            mk_ana = 'mkdir -p ' + germ_ana_dir + '/LOGS'
            sys.stderr.write('Creating analysis output directories ' + mk_ana + '\n')
            call(mk_ana, shell=True)
        os.chdir(germ_ana_dir)
        check = platypus_germline(config_file, sample, 'LOGS/', 'n')
        if check != 0:
            sys.stderr.write(date_time() + 'platypus failed around ' + sample + '\n')
            exit(1)
        set_acls(germ_ana_dir, user, group)
    if estep == 'start' or estep == 'annot':
        germ_ann_dir = project_dir + project + '/' + annotation + '/' + sample
        if not os.path.isdir(germ_ann_dir):
            mk_ann = 'mkdir -p ' + germ_ann_dir + '/LOGS'
            sys.stderr.write('Creating analysis output directories ' + mk_ann + '\n')
            call(mk_ann, shell=True)
        os.chdir(germ_ann_dir)
        check = annot_platypus(config_file, sample, 'n')
        if check == 0:
            sys.stderr.write(date_time() + 'Germ line call complete\n')
        else:
            sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
            exit(1)
        reorg = 'mv *.log LOGS;'
        call(reorg, shell=True)
        sys.stderr.write('Reorganizing germline log files ' + reorg + '\n')
        set_acls(germ_ann_dir, user, group)
    sys.stderr.write(date_time() + 'Standalone germline pipe starting at ' + estep + ' for sample ' + sample
                     + ' complete. Check logs and outputs.\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mini temp pipe to run a part of variant_annot_pipe')
    parser.add_argument('-s', '--sample', action='store', dest='sample',
                        help='Sample ID')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-e', '--estep', action='store', dest='estep',
                        help='Step to start at, \'start\' or \'annot\'')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample, config_file, estep) = (
        inputs.sample, inputs.config_file, inputs.estep)
    temp_germline_pipe(config_file, sample, estep)
