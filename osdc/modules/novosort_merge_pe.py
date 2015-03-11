#!/usr/bin/python
import sys
import os
import re
import json
from date_time import date_time
from subprocess import call
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from upload_to_swift import upload_to_swift
from download_from_swift import download_from_swift
import subprocess

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['novosort'],config_data['refs']['obj'],config_data['refs']['cont'])

def list_bam(obj,sample,wait,cont):
    ct=0
    list_cmd='. /home/ubuntu/.novarc;swift list ' + obj + ' --prefix ' +  cont + '/' + sample
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist=subprocess.check_output(list_cmd,shell=True)
    # Use to check on download status
    p=[]

    bam_list=''
    bai_list=''
    for fn in re.findall('(.*)\n',flist):
        if re.match('^\S+_\d+\.rmdup.srt.ba[m|i]$', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd='. /home/ubuntu/.novarc;swift download ' + obj + ' --skip-identical ' + fn
            p.append(subprocess.Popen(dl_cmd,shell=True))
            if fn[-3:] == 'bam':
                bam_list=bam_list + fn + ' '
                ct=ct+1
            else:
                bai_list=bai_list + fn + ' '
    bam_list=bam_list.rstrip(' ')
    bai_list=bai_list.rstrip(' ')
    n=0
    f=0
    x=len(p)
    while(n < wait):
        sys.stderr.write(date_time() + 'Checking status of download processes. ' + str(n) + ' seconds have passed\n')
        s=0
        for cur in p:
            check=cur.poll()
            if str(check)!= 'None':
                s=s+1
        if s==x:
            f=1
            break
        sys.stderr.write(date_time() + str(s) + ' of ' + str(x) + ' downloads have been completed\n')
        n=n+30
        sleep_cmd='sleep 30s;'
        subprocess.call(sleep_cmd,shell=True)
    if f==1:
        sys.stderr.write(date_time() + 'BAM download complete\n')
        return (bam_list,bai_list,ct)
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)

def novosort_merge_pe(config_file,sample_list,wait):
    fh=open(sample_list,'r')
    (novosort,obj,cont)=parse_config(config_file)
    for sample in fh:
        sample=sample.rstrip('\n')
        (bam_list,bai_list,n)=list_bam(obj,sample,wait,cont)
        if n > 1:
            novosort_merge_pe_cmd=novosort + " --threads 8 --ram 28G --assumesorted --output " + sample + '.merged.bam --index --tmpdir ./TMP ' + bam_list
            sys.stderr.write(date_time() + novosort_merge_pe_cmd + "\n")
            try:
                subprocess.check_output(novosort_merge_pe_cmd,shell=True) 
            except:
                sys.stderr.write(date_time() + 'novosort failed for sample ' + sample + '\n')
                exit(1)
        else:
            rename_bam='cp ' + bam_list + ' ' + sample + '.merged.final.bam;cp ' + bai_list + ' ' + sample + '.merged.final.bai'
            sys.stderr.write(date_time() + rename_bam + ' Only one associated bam file, renaming\n')
            subprocess.call(rename_bam,shell=True)
    sys.stderr.write(date_time() + 'Merge process complete\n')
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl','--sample_list',action='store',dest='sample_list',help='Sample/project name prefix list')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and ref locations')
    parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time to download bam files.  900 (seconds) recommended')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (sample_list,config_file,wait)=(inputs.sample_list,inputs.config_file,inputs.wait)
    novosort_merge_pe(config_file,sample_list,wait)
