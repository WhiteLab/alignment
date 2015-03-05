#!/usr/bin/python
import sys
import os
import re
import json
from date_time import date_time
from subprocess import call
from upload_to_swift import upload_to_swift
from download_from_swift import download_from_swift
import subprocess
sys.path.append('/home/ubuntu/TOOLS/Scripts/modules')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['novosort'])

def list_bam(obj,sample,wait):
    sys.stderr.write(date_time() + 'Getting BAM list\n')
    list_cmd='. /home/ubuntu/.novarc;swift list ' + obj + ' --prefix ' +  sample
    flist=subprocess.check_output(list_cmd,shell=True)
    # Use to check on download status
    p=[]

    bam_list=''
    for fn in re.findall('(.*)\n',flist):
        if re.match('^' + sample + '_\d+\.rmdup.srt.ba[m|i]$', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd='. /home/ubuntu/.novarc;swift download ' + obj + ' --skip-identical ' + fn
            p.append(subprocess.Popen(dl_cmd,shell=True))
            if fn[-3:] == 'bam':
                bam_list=bam_list + fn + ' '
    bam_list=bam_list.rstrip(' ')
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
        return bam_list
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)

def novosort_merge_pe(config_file,obj,sample,wait):
    bam_list=list_bam(obj,sample,wait)
    (novosort)=parse_config(config_file)
    novosort_merge_pe_cmd=novosort + " --threads 8 --ram 28G --assumesorted --output " + sample + '.merged.bam --index --tmpdir ./TMP ' + bam_list
    log(log_dir + sample + ".novosort.merge.pe.log",date_time() + novosort_merge_pe_cmd + "\n")
    try:
        subprocess.check_output(novosort_merge_pe_cmd,shell=True) 
    except:
        sys.stderr.write(date_time() + 'novosort failed for sample ' + sample + '\n')
        exit(1)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-o','--obj',action='store',dest='obj',help='Swift object name')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and ref locations')
    parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time to download bam files.  900 (seconds) recommended')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (sample,obj,config_file,wait)=(inputs.sample,inputs.obj,inputs.config_file,inputs.wait)
    novosort_merge_pe(config_file,obj,sample,wait)
