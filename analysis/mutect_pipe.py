#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.job_manager import job_manager
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['java'], config_data['tools']['mutect'], config_data['refs']['genome'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads'], config_data['params']['ram'], \
           config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir']


def mutect_pipe(config_file, tumor_id, normal_id):
    (java, mutect, intervals, fa_ordered, max_t, ram, project_dir, project, align) = parse_config(config_file)

    # break up intervals into max threads junks to run all in parallel
    int_fh = open(intervals, 'r')
    int_dict = {}
    i = 0
    # create temp directory
    tmp_cmd = 'mkdir temp'
    subprocess.call(tmp_cmd, shell=True)
    # create sub-interval files - split by chromosome
    mk_dir_bed = 'mkdir bed'
    subprocess.call(mk_dir_bed, shell=True)
    for interval in int_fh:
        (chrom, start, end) = interval.split('\t')
        try:
            int_dict[chrom]['fh'].write(interval)
        except:
            int_dict[chrom] = {}
            int_dict[chrom]['fn'] = 'bed/intervals_' + chrom + '.bed'
            int_dict[chrom]['fh'] = open(int_dict[chrom]['fn'], 'w')
            int_dict[chrom]['fh'].write(interval)
        i += 1
    job_ram = int(int(ram) / int(max_t))
    run_mut = java + ' -Djava.io.tmpdir=./temp -Xmx' + str(job_ram) + 'g -jar ' + mutect
    # array will store commands to run, next def will take care of job management using popen
    cmd_list = []
    bam_dir = project_dir + project + '/' + align
    tumor_bam = bam_dir + '/' + tumor_id + '/BAM/' + tumor_id + '.merged.final.bam'
    normal_bam = bam_dir + '/' + normal_id + '/BAM/' + normal_id + '.merged.final.bam'
    sys.stderr.write(date_time() + 'Processing pair T: ' + tumor_bam + ' N: ' + normal_bam + '\n')
    out = tumor_id + '_' + normal_id
    # make result directory for current pair
    i = 1
    for intvl in sorted(int_dict):
        int_dict[intvl]['fh'].close()
        cur = run_mut
        output_file = out + '.' + intvl + '.out'
        vcf_file = out + '.' + intvl + '.vcf'
        log_file = 'LOGS/' + out + '.mut.' + intvl + '.log'
        cur = cur + ' -T MuTect -fixMisencodedQuals -R ' + fa_ordered + ' --intervals ' + int_dict[intvl][
            'fn'] + '  --input_file:normal ' + normal_bam + '  --input_file:tumor ' + tumor_bam + \
              ' --max_alt_alleles_in_normal_count 1000 --max_alt_alleles_in_normal_qscore_sum 37 ' \
              '--max_alt_allele_in_normal_fraction 0.05 --out ' + output_file + ' -vcf ' + vcf_file \
              + ' --enable_extended_output --strand_artifact_power_threshold 0 -log ' + log_file \
              + ' >> ' + log_file + ' 2>> ' + log_file + '; cat ' + output_file \
              + ' | grep -v REJECT > ' + output_file + '.keep; cat ' + vcf_file \
              + ' | grep -v REJECT > ' + vcf_file + '.keep '
        cmd_list.append(cur)
        i += 1
    # fix encode flag won't work if already phred 33, if a job fails try without
    try:
        job_manager(cmd_list, max_t)
    except:
        for i in range(0, len(cmd_list), 1):
            cmd_list[i] = cmd_list[i].replace('-fixMisencodedQuals ', '')
        job_manager(cmd_list, max_t)
    cleanup_temp_dirs = 'rm -rf temp bed'
    sys.stderr.write('Cleaning up temp dirs ' + cleanup_temp_dirs + '\n')
    subprocess.call(cleanup_temp_dirs, shell=True)
    sys.stderr.write(date_time() + 'SNV calling completed for ' + out + '\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='muTect pipleine for variant calling.  Need BAM and bai files ahead of time.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-t', '--tumor', action='store', dest='tumor', help='Tumor sample id')
    parser.add_argument('-n', '--normal', action='store', dest='normal', help='Normal sample id')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, tumor_id, normal_id) = (inputs.config_file, inputs.tumor, inputs.normal)
    mutect_pipe(config_file, tumor_id, normal_id)
