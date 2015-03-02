#!/usr/bin/python
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
import os
import re
import argparse
from date_time import date_time
import subprocess
from download_from_swift import download_from_swift
from pipeline import Pipeline
import pdb
from log import log

parser=argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially.')
parser.add_argument('-f','--file',action='store',dest='fn',help='File with bionimbus ID, seqtype and sample list')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

inputs=parser.parse_args()
fh=open(inputs.fn,'r')
src_cmd='. ~/.novarc;'
obj='PANCAN'
cont='ALIGN_TEST'
# paired-end config file to use for pipeline
pipe_cfg='/home/ubuntu/TOOLS/Scripts/utility/hg19_pe_config.json'


for line in fh:
    line=line.rstrip('\n')
    (bid,seqtype,lane_csv)=line.split('\t')
    cwd='/mnt/cinder/REFS_' + bid + '/SCRATCH'
    loc=cwd[:-7] + bid + '.run.log'
    log(loc,date_time() + 'Initializing scratch directory for ' + bid + '\n')
    # All files for current bid to be stored in cwd
    check_dir=os.path.isdir(cwd)
    if check_dir==False:
        subprocess.check_output('mkdir ' + cwd,shell=True)
    try:
        os.chdir(cwd)
    except:
        log(loc,date_time() + 'Creating directory for ' + bid + ' failed. Ensure correct machine being used for this sample set\n')
        continue
    
    contain='RAW/' + bid + '/' + bid + '_'
    cur_dir=cwd + '/RAW/' + bid
    # iterate through sample/lane pairs
    # dictionary to track status of success of pipelines for each sample and lane to help troubleshoot any failures
    lane_status={}
    for lane in lane_csv.split(', '):
        lane_status[lane]='Initializing'
        swift_cmd=src_cmd + 'swift list ' + obj + ' --prefix ' + contain + lane
        log(loc,date_time() + 'Getting sequence files for lane ' + lane + '\n' + swift_cmd + '\n')
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
        seqfile=re.match('^(\S+)\n(\S+)\n$',contents)
        sf1=seqfile.group(1)
        end1=os.path.basename(sf1)
        sf2=seqfile.group(2)
        end2=os.path.basename(sf2)
        lane_status[lane]='Downloading'
        prefix='RAW/' + bid + '/' + bid + '_' + lane

        # attempt to download sequencing files
        try:
            download_from_swift(obj,prefix)
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
        p=Pipeline(end1,end2,seqtype,pipe_cfg)
        if str(p)!='0':
            log(loc,date_time() + "Pipeline process for sample lane " + lane + " failed with status " + str(p) + " \n")
            lane_status[lane]='Pipeline return status failed'
            log(loc,lane + '\t' + lane_status[lane] + '\n')
            exit(3)
        # change back to parent directory so that new sequencing files can be downloaded in same place
        os.chdir(cwd)
        lane_status[lane]='Success!'
        log(loc,lane + '\t' + lane_status[lane] + '\n')
    os.chdir(cur_dir)
    mv_gz='mv ../*.gz .'
    subprocess.call(mv_gz,shell=True)
    qc_cmd='/home/ubuntu/Scripts/parse_qc.pl 0'
    subprocess.call(qc_cmd)
    upload_to_swift(obj,cont)
sys.stderr.write(date_time() + "Process complete.  Check logs for any errors\n")
