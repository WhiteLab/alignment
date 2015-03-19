#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'],config_data['tools']['mutect'],config_data['refs']['intervals'],config_data['refs']['fa_ordered'])

def job_manage(cmd_list,out,max_t):
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
        p[i]['class']=subprocess.Popen(cmd_list[i],shell=True)
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
            check=p[i]['class'].poll()
            if str(check) == '1':
                sys.stderr.write(date_time() + 'Job returned an error while running ' + p[i]['cmd'] + '  aborting!\n')
                exit(1)
            if str(check) == '0':
                comp+=1
                if comp <= (x-n):
                    p[i]['class']=subprocess.Popen(cmd_list[cur],shell=True)
                    p[i]['cmd']=cmd_list[cur]
                    cur+=1
        s+=j
        sleep_cmd='sleep ' + str(j) + 's'
        subprocess.call(sleep_cmd,shell=True)
    sys.stderr.write(date_time() + 'Jobs completed for ' + out + '\n')
def mutect_pipe(config_file,sample_pairs,ref_mnt):
    max_t=8
    (java,mutect,intervals,fa_ordered)=parse_config(config_file)
    intervals=ref_mnt + '/' + intervals
    #break up intervals into max threads junks to runn all in parallel
    int_fh=open(intervals,'r')
    int_dict={}
    i=0
    # create temp directory
    tmp_cmd='mkdir temp'
    subprocess.call(tmp_cmd,shell=True)
    # create sub-interval files - split by chromosome
    for interval in int_fh:
        (chrom,intvl)=interval.split(':')
        try:
            int_dict[chrom]['fh'].write(interval)
        except:
            int_dict[chrom]={}
            int_dict[chrom]['fn']='intervals_' + chrom + '.list'
            int_dict[chrom]['fh']=open(int_dict[chrom]['fn'],'w')
            int_dict[chrom]['fh'].write(interval)
        i+=1
    fa_ordered=ref_mnt+ '/' + fa_ordered
    fh=open(sample_pairs)
    run_mut=java + ' -Djava.io.tmpdir=./temp -Xmx2g -jar ' + mutect
    mk_log_dir='mkdir LOGS'
    subprocess.call(mk_log_dir,shell=True)
    for line in fh:
        #array will store commands to run, next def will take care of job management using popen
        cmd_list=[]
        line=line.rstrip('\n')
        (sample,tumor_id,normal_id)=line.split('\t')
        tumor_bam=tumor_id + '.merged.final.bam'
        normal_bam=normal_id + '.merged.final.bam'
        sys.stderr.write(date_time() + 'Processing pair T: ' + tumor_bam + ' N: ' + normal_bam + '\n' )
        out=tumor_id + '_' + normal_id
        # make result directory for current pair
        mk_res='mkdir ' + out
        subprocess.call(mk_res,shell=True)
        i=1
        for intvl in sorted(int_dict):
            int_dict[intvl]['fh'].close()
            cur=run_mut
            #            sys.stderr.write('interval: ->' + interval + '<-\n')
            output_file=out + '.' + intvl + '.out'
            vcf_file=out + '.' + intvl +  '.vcf'
            log_file='LOGS/' + out + '.mut.' + intvl +  '.log'
            cur=cur+ ' -T MuTect -fixMisencodedQuals -R ' + fa_ordered + ' --intervals ' + int_dict[intvl]['fn'] + '  --input_file:normal ' + normal_bam + '  --input_file:tumor ' + tumor_bam + ' --out ' + out + '/' + output_file + ' -vcf ' + out + '/' + vcf_file + ' --enable_extended_output --strand_artifact_power_threshold 0 -log ' + log_file + ' >> ' + log_file + ' 2>> ' + log_file + '; cat ' + out + '/' + output_file + ' | grep -v REJECT > ' + out + '/' + output_file + '.keep; cat ' + out + '/' + vcf_file + ' | grep -v REJECT > ' + out + '/' + vcf_file + '.keep '
            cmd_list.append(cur)
            i=i+1
        job_manage(cmd_list,out,max_t)
    sys.stderr.write(date_time() + 'Variant calling completed!\n')
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='muTect pipleine for variant calling.  Need BAM and bai files ahead of time.')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file with tool and reference locations')
    parser.add_argument('-sp','--sample_pairs',action='store',dest='sample_pairs',help='Sample tumor/normal pairs')
    parser.add_argument('-r','--ref_mnt',action='store',dest='ref_mnt',help='Reference drive path - i.e. /mnt/cinder/REFS_XXXX')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (config_file,sample_pairs,ref_mnt)=(inputs.config_file,inputs.sample_pairs,inputs.ref_mnt)
    mutect_pipe(config_file,sample_pairs,ref_mnt)