#!/usr/bin/python
import sys
from date_time import date_time
import subprocess
import json

def parse_config(config_file):
    config_data=json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'],config_data['tools']['mutect'],config_data['refs']['intervals'],config_data['refs']['fa_ordered'])

def job_manage(cmd_list):
    
    for cmd in cmd_list:
        sys.stderr.write(cmd + '\n')

def muTect_pipe(config_file,sample_pairs,ref_mnt):
    (java,mutect,intervals,fa_ordered)=parse_config(config_file)
    intervals=ref_mnt + '/' + intervals
    fa_ordered=ref_mnt+ '/' + fa_ordered
    fh=open(sample_pairs)
    run_mut=java + ' -Djava.io.tmpdir=' + ref_mnt + '/temp -Xmx2g -jar ' + mutect
    mk_log_dir='mkdir logs'
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
        int_fh=open(intervals,'r')
        for interval in int_fh:
            cur=run_mut
#            sys.stderr.write('interval: ->' + interval + '<-\n')
            interval=interval.rstrip('\n')
            (chrom,intvl)=interval.split(':')
            output_file=out + '.' + chrom + '_' + intvl + '.out'
            vcf_file=out + '.' + chrom  + '_' + intvl + '.vcf'
            log_file=out + '.' + chrom + '_' + intvl + '.log'
            cur=cur+ ' -T MuTect -R ' + fa_ordered + ' --intervals ' + interval + '  --input_file:normal ' + normal_bam + '  --input_file:tumor ' + tumor_bam + ' --out ' + output_file + ' -vcf ' + vcf_file + ' --enable_extended_output --strand_artifact_power_threshold 0 -log ' + log_file + '; cat ' + output_file + ' | grep -v REJECT > ' + output_file + '.keep; cat ' + vcf_file + ' | grep -v REJECT > ' + vcf_file + '.keep'
            cmd_list.append(cur)
        job_manage(cmd_list)

            
        
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
    muTect_pipe(config_file,sample_pairs,ref_mnt)
