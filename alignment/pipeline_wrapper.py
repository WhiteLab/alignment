#!/usr/bin/env python
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import os
import re
import argparse
import json
from utility.date_time import date_time
from utility.find_project_files import find_project_files
from subprocess import call

parser = argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially.')
parser.add_argument('-f', '--file', action='store', dest='fn',
                    help='File with bionimbus ID, seqtype and sample lane list')
parser.add_argument('-j', '--json', action='store', dest='config_file',
                    help='JSON config file with tools, references, and data storage locations')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

inputs = parser.parse_args()
fh = open(inputs.fn, 'r')


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project'], config_data['refs']['align'], config_data['refs']['config'], \
           config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['align_pipe'], \
           config_data['tools']['slurm_wrap']


(project, align_dir, pipe_cfg, cores, mem, align_pipe, slurm_wrap) = parse_config(inputs.config_file)
cwd = '/cephfs/PROJECTS/' + project
check_dir = os.path.isdir(cwd)
for line in fh:
    line = line.rstrip('\n')
    (bid, seqtype, lane_csv) = line.split('\t')
    try:
        os.chdir(cwd)
    except:
        sys.stderr.write(
            date_time() + 'Could not find working directory ' + cwd + '. Ensure correct project was set in config\n')
        exit(1)
    # All files for current bid to be stored in cwd

    cur_dir = cwd + '/RAW/' + bid
    # iterate through sample/lane pairs
    # dictionary to track status of success of pipelines for each sample and lane to help troubleshoot any failures
    for lane in lane_csv.split(', '):
        seq_dir = 'RAW/' + bid
        file_prefix = bid + '_' + lane
        (contents, seqfile, sf1, sf2) = ('', [], '', '')
        # attempt to find sequencing files
        try:
            sys.stderr.write(date_time() + 'Searching for sequencing files related to ' + lane + '\n')
            contents = find_project_files(seq_dir, file_prefix)
            # sequencing files found in pairs using simple iterator, as find gives files in alphanumeric order -
            # standard file naming should work with this
            seqfile = re.findall('(\S+[sequence|f*q]*\.gz)', contents)
            sf1 = seqfile[0]
            sf2 = seqfile[1]
        except:
            sys.stderr.write(date_time() + 'Getting sequencing files ' + sf1 + ' and ' + sf2 + ' failed.  Moving on\n')
            continue

        # Create sbatch script and submit
        #p = Pipeline(end1, end2, seqtype, pipe_cfg)
        # batch params $cores $mem $pipeline $f1 $f2 $t $j
        job_log = lane + '.log'
        batch = 'sbatch ' + slurm_wrap + ' --export=cores="' + cores + '",mem="' + mem + '",log="' + job_log \
                + '",pipeline="' + align_pipe + '",f1="' + sf1 + '",f2="' + sf2 + '",t="' + seqtype \
                + '",j="' + pipe_cfg + '"'
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + lane + ' failed! Check logs!\n')
        # change back to parent directory so that new sequencing files can be searched

sys.stderr.write(date_time() + "Jobs submitted.  Check logs for any errors\n")
