#!/usr/bin/python
import sys
import json
import subprocess
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'],config_data['tools']['picard'],config_data['refs']['fa_ordered'])

def ksort(config_file,bam_list,flag,ref_mnt):
    (java,picard,fa)=parse_config(config_file)
    k_proc=[]
    c=0
    s=0
    fh=open(bam_list,'r')
    blist=[]
    for bam in fh:
        bam=bam.rstrip('\n')
        blist.append(bam)
        ksort_cmd=java + ' -Xmx4g -jar ' + picard + ' ReorderSam INPUT=' + bam + '  REFERENCE=' + ref_mnt + '/' + fa + ' OUTPUT=' + bam + '.reordered CREATE_INDEX=True'
        sys.stderr.write(date_time() + 'Started reorder process for ' + bam + '\n' + ksort_cmd + '\n')
        k_proc.append(subprocess.Popen(ksort_cmd,shell=True))
    n=len(blist)
    while(c < n):
        sys.stderr.write(date_time() + 'Checking status of reorder processes. ' + str(s) + ' seconds have passed\n')
        for cur in k_proc:
            check=cur.poll()
            sys.stderr.write(str(check) + '\n')
            if str(check) == '0':
                c=c+1
                k_proc.remove(cur)
            if str(check) == '1':
                sys.stderr.write(date_time() + 'Picard command failed.  Check references and bam file names and try again!\n')
                exit(1)
        sys.stderr.write(date_time() + str(c) + ' of ' + str(n) + ' processes have been completed\n')
        s=s+20
        sleep_cmd='sleep 20s;'
        subprocess.call(sleep_cmd,shell=True)
    if flag == 'y':
        for bam in blist:
            bai=bam[:-1] + 'i'
            sys.stderr.write('Renaming ' + bam + '.reordered to ' + bam + ' and ' + bam + '.reordered.bai to ' + bai)
            rn='mv ' + bam + '.reordered ' + bam + ';mv ' + bam + '.reordered.bai ' + bai 
            subprocess.call(rn,shell=True)
    sys.stderr.write(date_time() + 'Reorder complete\n')
    return 0
        

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Picard tool to reorder BAM file by karyotypic order, necessary for running muTect.')
    parser.add_argument('-b','--bam_list',action='store',dest='bam_list',help='BAM file list')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and ref locations')
    parser.add_argument('-o','--overwrite',action='store',dest='flag',help='Enter \'y\' or \'n.\' Flag to overwrite original after reordering or not.')
    parser.add_argument('-r','--reference',action='store',dest='ref_mnt',help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (bam_list,config_file,flag,ref_mnt)=(inputs.bam_list,inputs.config_file,inputs.flag,inputs.ref_mnt)
    ksort(config_file,bam_list,flag,ref_mnt)
