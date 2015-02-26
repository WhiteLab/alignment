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
# dictionary to track status of success of pipelines for each sample and lane to help troubleshoot any failures
samp_status={}
for line in fh:
    line=line.rstrip('\n')
    (bid,seqtype,samp_csv)=line.split('\t')
    cwd='/mnt/cinder/REFS_' + bid + '/SCRATCH'
    sys.stderr.write(date_time() + 'Initializing scratch directory for ' + bid + '\n')
    # All files for current bid to be stored in cwd
    check_dir=os.path.isdir(cwd)
    if check_dir==False:
        subprocess.check_output('mkdir ' + cwd,shell=True)
    try:
        os.chdir(cwd)
    except:
        sys.stderr.write(date_time() + 'Creating directory for ' + bid + ' failed. Ensure correct machine being used for this sample set\n')
        continue
    
    contain='RAW/' + bid + '/' + bid + '_'
    cur_dir=cwd + '/RAW/' + bid
    # get list of paired end sequencing files for each pair
    for samp in samp_csv.split(', '):
        samp_status[samp]={}
        samp_status[samp]['s']='Initializing'
        samp_status[samp]['sf']={}
        swift_cmd=src_cmd + 'swift list ' + obj + ' --prefix ' + contain + samp
        sys.stderr.write(date_time() + 'Getting sequece files for sample ' + samp + '\n' + swift_cmd + '\n')
        try:
            contents=subprocess.check_output(swift_cmd,shell=True)
        except:
            sys.stderr.write(date_time() + 'Can\'t find sequencing files for ' + samp + ' skipping!\n')
            continue
        end1=''
        end2=''
        sf1=''
        sf2=''
        i=0
        samp_status[samp]['s']='Running'
        # sequencing files downloaded in pairs using simple iterator, as swift gives files in alphanumeric order - standard file naming should work with this
        for seqfile in re.findall('(.*)\n',contents):
            seqfile=seqfile.rstrip('\n')
            if i==0:
                sf1=seqfile
                end1=os.path.basename(sf1)
                i=i+1
            else:
                sf2=seqfile
                end2=os.path.basename(sf2)
                samp_status[samp]['sf'][seqfile]='Downloading'
                i=0
                m=re.match('^(\S+)_2_sequence\.txt\.gz$',sf2)
                lane=m.group(1)
                # attempt to download sequencing files
                sys.stderr.write(date_time() + 'Getting sequencing files\n')
                try:
                    download_from_swift(obj,lane)
                except:
                    sys.stderr.write(date_time() + 'Getting sequencing files ' + sf1 + ' and ' + sf2 + ' failed.  Moving on\n')
                    samp_status[samp]['sf'][seqfile]='Download failed'
                    continue
                # pipeline needs to be run in same directory as sequencing files
                try:
                    os.chdir(cur_dir)
                    l_dir=cur_dir + '/LOGS'
                    l_check=os.path.isdir(l_dir)
                    if l_check==False:
                        subprocess.call('mkdir ' + l_dir, shell=True)
                except:
                    sys.stderr.write(date_time() + 'Could not change to new directory ' + cur_dir + ' Skipping and removing sequencing files\n')
                    rm_sf='rm ' + cur_dir + '/' + sf1 + ' ' + cur_dir + '/' +  sf2
                    subprocess.call(rm_sf,shell=True)
                    os.chdir(cwd)
                    exit(3)
                # if pipeline fails, abandon process as a larger error might come up
                sys.stderr.write(date_time() + 'Running pipeline process \n')
                try:
                    p=Pipeline(end1,end2,seqtype,pipe_cfg)
                    p.pipeline()
                    if p==1:
                        sys.stderr.write("Pipeline process for sample " + seqfile + " failed!\n")
                        samp_status[samp]['sf'][seqfile]='Pipeline failed'
                        samp_status[samp]['s']='Failures detected'
                        for sf in samp_status[samp]['sf']:
                            sys.stderr.write(samp + '\t' + sf + '\t' + samp_status[samp]['sf'][sf] + '\n')
                            exit(3)
                except:
                    sys.stderr.write("Pipeline process for sample " + seqfile + " failed!\n")
                    samp_status[samp]['sf'][seqfile]='Pipeline failed'
                    samp_status[samp]['s']='Failures detected'
                    for sf in samp_status[samp]['sf']:
                        sys.stderr.write(samp + '\t' + sf + '\t' + samp_status[samp]['sf'][sf] + '\n')
                    exit(1)
                # change back to parent directory so that new sequencing files can be downloaded in same place
                os.chdir(cwd)
                samp_status[samp]['sf'][seqfile]='Success!'
        samp_status[samp]['s']='Succeeded'
        for sf in samp_status[samp]['sf']:
            sys.stderr.write(samp + '\t' + sf + '\t' + samp_status[samp]['sf'][sf] + '\n')
    os.chdir(cur_dir)
    mv_gz='mv ../*.gz .'
    subprocess.call(mv_gz,sgell=True)
    qc_cmd='/home/ubuntu/Scripts/parse_qc.pl 0'
    subprocess.call(qc_cmd)
    upload_to_swift(obj,cont)
sys.stderr.write(date_time() + "Process complete.  Check logs for any errors\n")
