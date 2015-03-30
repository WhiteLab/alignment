#!/usr/bin/python
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
import os
import re
import argparse
import json
from date_time import date_time
import subprocess
from download_from_swift import download_from_swift
from pipeline import Pipeline
import pdb
from log import log

parser=argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially.')
parser.add_argument('-f','--file',action='store',dest='fn',help='File with bionimbus ID, seqtype and sample lane list')
parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tools, references, and data storage locations')
parser.add_argument('-m','--mount',action='store',dest='ref_mnt',help='Reference drive mount location.  Example would be /mnt/cinder/REFS_XXX')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputs=parser.parse_args()
fh=open(inputs.fn,'r')
src_cmd='. ~/.novarc;'
ref_mnt=inputs.ref_mnt
def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['refs']['cont'],config_data['refs']['obj'],config_data['refs']['config'])

(cont,obj,pipe_cfg)=parse_config(inputs.config_file)

for line in fh:
    line=line.rstrip('\n')
    (bid,seqtype,lane_csv)=line.split('\t')
    cwd=ref_mnt + '/SCRATCH'
    check_dir=os.path.isdir(cwd)
    if check_dir==False:
        subprocess.check_output('mkdir ' + cwd,shell=True)
    try:
        os.chdir(cwd)
    except:
        sys.stderr.write(date_time() + 'Creating directory for ' + bid + ' failed. Ensure correct machine being used for this sample set\n')
        exit(1)
    loc=cwd[:-7] + bid + '.run.log'
    log(loc,date_time() + 'Initializing scratch directory for ' + bid + '\n')
    # All files for current bid to be stored in cwd
    
    obj1='RAW/' + bid + '/' + bid + '_'
    cur_dir=cwd + '/RAW/' + bid
    # iterate through sample/lane pairs
    # dictionary to track status of success of pipelines for each sample and lane to help troubleshoot any failures
    lane_status={}
    for lane in lane_csv.split(', '):
        lane_status[lane]='Initializing'
        swift_cmd=src_cmd + 'swift list ' + cont + ' --prefix ' + obj1 + lane
        log(loc,date_time() + 'Getting sequence files for sample ' + lane + '\n' + swift_cmd + '\n')
        try:
            contents=subprocess.check_output(swift_cmd,shell=True)
        except:
            log(loc,date_time() + 'Can\'t find sequencing files for ' + lane + ' skipping!\n')
            continue
        end1=''
        end2=''
        sf1=''
        sf2=''

        lane_status[lane]='Running'
        # sequencing files downloaded in pairs using simple iterator, as swift gives files in alphanumeric order - standard file naming should work with this
        seqfile=re.findall('(\S+[sequence|f*q]*\.gz)',contents)
        sf1=seqfile[0]
        end1=os.path.basename(sf1)
        sf2=seqfile[1]
        end2=os.path.basename(sf2)
        lane_status[lane]='Downloading'
        prefix='RAW/' + bid + '/' + bid + '_' + lane
        # attempt to download sequencing files
        try:
            download_from_swift(cont,prefix)
        except:
            log(loc,date_time() + 'Getting sequencing files ' + sf1 + ' and ' + sf2 + ' failed.  Moving on\n')
            lane_status[lane]='Download failed'
            continue
            # pipeline needs to be run in same directory as sequencing files
        if os.path.isfile(cur_dir + '/' + end1) and os.path.isfile(cur_dir + '/' + end2):
            lane_status[lane]='Sequencing file download successful'
        else:
            lane_status[lane]='Sequencing file download failed'
            log(loc,lane + '\t' + lane_status[lane] + '\n')
            exit(3)
        try:
            os.chdir(cur_dir)
            l_dir=cur_dir + '/LOGS'
            l_check=os.path.isdir(l_dir)
            if l_check==False:
                subprocess.call('mkdir ' + l_dir, shell=True)
        except:
            log(loc,date_time() + 'Could not change to new directory ' + cur_dir + ' Skipping and removing sequencing files\n')
            rm_sf='rm ' + cur_dir + '/' + end1 + ' ' + cur_dir + '/' +  end2
            subprocess.call(rm_sf,shell=True)
            os.chdir(cwd)
            exit(3)
            # if pipeline fails, abandon process as a larger error might come up
        log(loc,date_time() + 'Running pipeline process for lane ' + lane + '\n')
        #check class status flag
        p=Pipeline(end1,end2,seqtype,pipe_cfg,ref_mnt)
        if p.status != 0:
            log(loc,date_time() + "Pipeline process for sample lane " + lane + " failed with status " + str(p.status) + " \n")
            lane_status[lane]='Pipeline return status failed'
            log(loc,lane + '\t' + lane_status[lane] + '\n')
            exit(3)
        # change back to parent directory so that new sequencing files can be downloaded in same place
        os.chdir(cwd)
        # clean out files for next run
        cleanup='rm -rf ' + cur_dir + '/BAM'
        subprocess.call(cleanup,shell=True)
        lane_status[lane]='Pipeline run and data uplaoded'
        log(loc,date_time() + lane + '\t' + lane_status[lane] + '\n')
    os.chdir(cur_dir)
sys.stderr.write(date_time() + "Process complete.  Check logs for any errors\n")
