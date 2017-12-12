#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from annotation.germline_report import gen_report
from utility.date_time import date_time
import subprocess
import json
import psutil
import os
import signal


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['fa_ordered'],\
           config_data['params']['threads'], config_data['tools']['snpsift'], config_data['tools']['java'], \
           config_data['refs']['cadd'], config_data['refs']['tx_index']


def pass_filter(sample):
    in_fn = sample + '.germline_calls.vcf'
    out_fn = sample + '.germline_pass.vcf'
    out = open(out_fn, 'w')
    infile = open(in_fn, 'r')
    for line in infile:
        if line[0] == '#':
            out.write(line)
        else:
            fields = line.split('\t')
            if fields[6] == 'PASS':
                out.write(line)
    infile.close()
    out.close()


def run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size):
    cmd = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf + ' --symbol --vcf_info_field ANN ' \
        '--canonical --html --variant_class --sift both --offline --maf_exac --no_whole_genome --buffer_size ' \
          + buffer_size + ' --fork ' + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache \
          + ' --plugin CADD,' + cadd + ' 2>> LOGS/' + sample + '.vep.log >> LOGS/' + sample + '.vep.log;'
    return cmd


def watch_mem(proc_obj, sample):
    from time import sleep
    while proc_obj.poll() is None:
        mem_pct = psutil.virtual_memory().percent
        sys.stderr.write(date_time() + 'Current memory usage at ' + str(mem_pct) + '% processing sample ' + sample
            + ' from platypus ' + '\n')
        if mem_pct >= 99:
            sys.stderr.write(date_time() + 'Memory exceeded while running VEP.')
            return 1
        sleep(30)

    return proc_obj.poll()


def annot_platypus(config_file, samp_list, ref_mnt):
    (vep_tool, vep_cache, fasta, threads, snpsift, java, cadd, tx_index) = parse_config(config_file)
    fasta = ref_mnt + '/' + fasta
    cadd = ref_mnt + '/' + cadd
    vep_cache = ref_mnt + '/' + vep_cache
    tx_index = ref_mnt + '/' + tx_index
    # parse sample file, use only last if pairs
    samp_fh = open(samp_list, 'r')
    # track to prevent repeat annotation if same sample used as comparison
    samp_flag = {}
    for line in samp_fh:
        info = line.rstrip('\n').split('\t')
        sample = info[0]
        if len(info) > 1:
            sample = info[2]
        if sample not in samp_flag:
            pass_filter(sample)
            mk_log_dir = 'mkdir LOGS'
            subprocess.call(mk_log_dir, shell=True)
            in_vcf = sample + '.germline_pass.vcf'
            out_vcf = sample + '.germline_pass.vep.vcf'
            buffer_size = '2000'
            run_cmd = run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size)
            sys.stderr.write(date_time() + 'Annotating sample ' + in_vcf + '\n')
            # from stack overflow to allow killing of spawned processes in main process fails for cleaner restart
            check = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
            check_run = watch_mem(check, sample)
            if check_run != 0:

                buffer_size = str(int(buffer_size) / 2)
                clean_up = 'rm ' + out_vcf + '*'
                sys.stderr.write(date_time() + 'VEP failed. Status of run was ' + str(check_run)
                                 + ' Trying smaller buffer size of ' + buffer_size + '\n' + clean_up + '\n')
                os.killpg(os.getpgid(check.pid), signal.SIGINT)

                subprocess.call(clean_up, shell=True)
                run_cmd = run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size)
                sys.stderr.write(date_time() + 'Annotating sample ' + sample + in_vcf + '\n')
                check = subprocess.call(run_cmd, shell=True)
                if check != 0:
                    sys.stderr.write(date_time() + 'VEP failed for sample ' + sample + '\n')
                    exit(1)
            else:
                sys.stderr.write(date_time() + 'VEP annotation of ' + in_vcf + ' successful!\n')

            check = gen_report(out_vcf, sample, tx_index)
            if check == 0:
                sys.stderr.write(date_time() + 'Summary table of germline calls completed!\n')
                samp_flag[samp_fh] = 1
            else:
                sys.stderr.write(date_time() + 'Summary table for ' + out_vcf + ' FAILED!\n')
                return 1

    return 0
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs,'
                                                                                           ' or sample list')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt) = (
        inputs.config_file, inputs.sample_pairs, inputs.ref_mnt)
    annot_platypus(config_file, sample_pairs, ref_mnt)
