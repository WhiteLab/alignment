#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['snpEff'],config_data['tools']['report'])

def job_manage(cmd_list,max_t):
    # num commands
    x = len(cmd_list)
    # cur position in command list
    cur=0
    #running
    r=0
    #completed
    comp=0
    # initialize process list
    p={}
    sys.stderr.write(date_time() + 'Initializing run\n')
    n=max_t
    if n > x:
        n=x
    for i in xrange(0,n,1):
        p[i]={}
        p[i]['snp']=subprocess.Popen(cmd_list[i],shell=True)
        p[i]['cmd']=cmd_list[i]
        sys.stderr.write(cmd_list[i]+ '\n')
        cur+=1
    s=0
    j=30
    m=30
    while comp <= x:
        if s % m == 0:
            sys.stderr.write(date_time() + 'Checking job statuses. ' + str(comp) + ' of ' + str(x) + ' completed. ' + str(s) + ' seconds have passed\n')
        for i in xrange(0,n,1):
            check=p[i]['snp'].poll()
            if str(check) == '1':
                sys.stderr.write(date_time() + 'Job returned an error while running ' + p[i]['cmd'] + '  aborting!\n')
                exit(1)
            if str(check) == '0':
                comp+=1
                if comp <= (x-n):
                    p[i]['snp']=subprocess.Popen(cmd_list[cur],shell=True)
                    p[i]['cmd']=cmd_list[cur]
                    cur+=1
        s+=j
        sleep_cmd='sleep ' + str(j) + 's'
        subprocess.call(sleep_cmd,shell=True)
def snpeff_pipe(config_file,sample_pairs):
    max_t=8
    (snpeff,report)=parse_config(config_file)
    fh=open(sample_pairs)
    mk_log_dir='mkdir logs'
    subprocess.call(mk_log_dir,shell=True)
    cmd_list=[]
    for line in fh:
        #array will store commands to run, next def will take care of job management using popen
        line=line.rstrip('\n')
        (sample,tumor_id,normal_id)=line.split('\t')
        run_snp=snpeff + ' ' + sample + '.out.keep  2> logs/' + sample + '.snpeff.log;' + report + ' -i ' + sample + '.out.keep.eff.vcf > ' + sample + '.vcf.keep.eff.xls'  
        cmd_list.append(run_snp)
    job_manage(cmd_list,max_t)
    sys.stderr.write(date_time() + 'SNP annotation  completed!\n')

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='muTect pipleine for variant calling.  Need BAM and bai files ahead of time.')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and reference locations')
    parser.add_argument('-sp','--sample_pairs',action='store',dest='sample_pairs',help='Sample tumor/normal pairs')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (config_file,sample_pairs)=(inputs.config_file,inputs.sample_pairs)
    snpeff_pipe(config_file,sample_pairs)
