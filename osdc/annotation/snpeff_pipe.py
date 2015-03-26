#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'],config_data['tools']['snpEff'],config_data['tools']['snpsift'],config_data['tools']['report'],config_data['refs']['dbsnp'])

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
def snpeff_pipe(config_file,sample_pairs,ref_mnt):
    max_t=8
    (java,snpeff,snpsift,report,dbsnp)=parse_config(config_file)
    dbsnp=ref_mnt + '/' + dbsnp
    fh=open(sample_pairs)
    mk_log_dir='mkdir LOGS'
    subprocess.call(mk_log_dir,shell=True)
    cmd_list=[]
    run_snpsift=java + ' -jar ' + snpsift + ' annotate ' + dbsnp
    run_snpeff= java + ' -jar ' + snpeff + ' eff -t hg19 '
    for line in fh:
        #array will store commands to run, next def will take care of job management using popen
        line=line.rstrip('\n')
        (sample,tumor_id,normal_id)=line.split('\t')
        #run snpsift first, then snpeff
        run_snp= run_snpsift  + ' ' + sample + '.out.keep > ' + sample + '.out.keep.sift.vcf 2> LOGS/' + sample + '.snpeff.log;'+ run_snpeff + ' ' + sample + '.out.keep.sift.vcf -v > ' + sample + '.out.keep.eff.vcf  2>> LOGS/' + sample + '.snpeff.log;' + report + ' -i ' + sample + '.out.keep.eff.vcf > ' + sample + '.out.keep.eff.xls'  
        cmd_list.append(run_snp)
    job_manage(cmd_list,max_t)
    sys.stderr.write(date_time() + 'SNP annotation  completed!\n')
    return 0
if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and reference locations')
    parser.add_argument('-sp','--sample_pairs',action='store',dest='sample_pairs',help='Sample tumor/normal pairs')
    parser.add_argument('-r','--ref_mnt',action='store',dest='ref_mnt',help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (config_file,sample_pairs,ref_mnt)=(inputs.config_file,inputs.sample_pairs,inputs.ref_mnt)
    snpeff_pipe(config_file,sample_pairs,ref_mnt)
