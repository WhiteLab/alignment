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
    return config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['params']['plugin_dir'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads'], config_data['tools']['java'], \
           config_data['refs']['cadd'], config_data['refs']['tx_index'], config_data['refs']['project_dir'], \
           config_data['refs']['project'], config_data['refs']['analysis']


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


def run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size, plugin_dir):
    cmd = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf + ' --symbol --vcf_info_field ANN ' \
        '--canonical --html --variant_class --sift both --offline --maf_exac --no_whole_genome --buffer_size ' \
          + buffer_size + ' --fork ' + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache + ' --dir_plugins '\
          + plugin_dir + ' --plugin CADD,' + cadd + ' 2>> ' + sample + '.vep.log >> ' + sample + '.vep.log;'
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


def annot_platypus(config_file, sample):
    (vep_tool, vep_cache, plugin_dir, fasta, threads, java, cadd, tx_index, project_dir, project, analysis) \
        = parse_config(config_file)
    ana_dir = project_dir + project + '/' + analysis + '/' + sample
    pass_filter(ana_dir + '/' + sample)
    in_vcf = ana_dir + '/' + sample + '.germline_pass.vcf'
    out_vcf = sample + '.germline_pass.vep.vcf'
    buffer_size = '2000'
    if int(threads) > 1:
        threads = str(int(threads) - 1)
    run_cmd = run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size, plugin_dir)
    sys.stderr.write(date_time() + 'Annotating sample ' + in_vcf + ' ' + run_cmd + '\n')
    # from stack overflow to allow killing of spawned processes in main process fails for cleaner restart
    check = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    check_run = watch_mem(check, sample)
    if check_run != 0:

        buffer_size = str(int(buffer_size) / 2)
        clean_up = 'rm \'' + out_vcf + '*\''
        sys.stderr.write(date_time() + 'VEP failed. Status of run was ' + str(check_run)
                         + ' Trying smaller buffer size of ' + buffer_size + '\n' + clean_up + '\n')
        try:
            os.killpg(os.getpgid(check.pid), signal.SIGINT)
        except:
            sys.stderr.write(date_time() + 'Killing process failed.  Might have already died for other reasons...\n')

        subprocess.call(clean_up, shell=True)
        run_cmd = run_vep(vep_tool, in_vcf, out_vcf, threads, fasta, vep_cache, cadd, sample, buffer_size, plugin_dir)
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
    else:
        sys.stderr.write(date_time() + 'Summary table for ' + out_vcf + ' FAILED!\n')
        return 1

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-s', '--sample', action='store', dest='sample', help='Normal sample ID')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample) = (inputs.config_file, inputs.sample)
    annot_platypus(config_file, sample)
