#!/usr/bin/python
import sys
import os
import re
import json
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
import subprocess

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['refs']['cont'],config_data['refs']['obj'],config_data['params']['threads'],config_data['params']['ram'])

def list_bam(cont,obj,sample,wait):
    ct=0
    list_cmd='. /home/ubuntu/.novarc;swift list ' + cont + ' --prefix ' +  obj + '/' + sample
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist=subprocess.check_output(list_cmd,shell=True)
    # Use to check on download status
    p=[]
    
    for fn in re.findall('(.*)\n',flist):
        if re.match('.*.merged.final.ba', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd='. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + fn + ' --output ' + os.path.basename(fn)
            p.append(subprocess.Popen(dl_cmd,shell=True))
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
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)
        
def get_merged_bams(config_file,sample_list,wait):
    fh=open(sample_list,'r')
    (cont,obj,threads,ram)=parse_config(config_file)
    for sample in fh:
        sample=sample.rstrip('\n')
        list_bam(cont,obj,sample,wait)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl','--sample_list',action='store',dest='sample_list',help='Sample/project prefix list')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and ref locations')
    parser.add_argument('-w','--wait',action='store',dest='wait',help='Wait time to download bam files.  900 (seconds) recommended')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (sample_list,config_file,wait)=(inputs.sample_list,inputs.config_file,inputs.wait)
    get_merged_bams(config_file,sample_list,wait)
