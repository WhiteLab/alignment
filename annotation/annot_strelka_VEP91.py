#!/usr/bin/env python3

import sys
import os
import signal
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from annotation.VEP91_strelka_report import gen_report as gen_strelka_report
from utility.date_time import date_time
from utility.log import log
import subprocess
import json
import psutil


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['fa_ordered'], \
           config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['params']['vep_cache_version'], \
           config_data['params']['plugin_dir'], config_data['params']['threads'], config_data['refs']['intervals'], \
           config_data['params']['wg_flag'], config_data['refs']['tx_index'],\
           config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['annotation']


def run_vep(vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc, plugin_dir):
    run_cmd = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf \
              + ' --symbol --sift b --vcf_info_field ANN --canonical --variant_class --buffer_size ' + buffer_size \
              + ' --offline --af_gnomad --hgvs --hgvsg --fork ' + threads + ' --fasta ' + fasta + ' --dir_cache ' \
              + vep_cache + ' --cache_version ' + vcache + ' --dir_plugins ' + plugin_dir + '  2>> ' + loc + ' >> ' \
              + loc

    return run_cmd


def watch_mem(proc_obj, source, sample, loc):
    from time import sleep
    while proc_obj.poll() is None:
        mem_pct = psutil.virtual_memory().percent
        log(loc, date_time() + 'Current memory usage at ' + str(mem_pct) + '% processing sample ' + sample
            + ' from source ' + source + '\n')
        if mem_pct >= 99:
            log(loc, date_time() + 'Memory exceeded while running VEP.')
            return 1
        sleep(30)

    return proc_obj.poll()


def annot_vcf_vep_pipe(config_file, sample_pair, in_suffix, out_suffix):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, plugin_dir, threads, intvl, wg_flag, tx_index,
     project_dir, project, annotation) = parse_config(config_file)
    # scale back on the forking a bit

    if int(threads) > 2:
        # threads = str(int(threads)/2 - 1)
        threads = str(int(threads) - 1)
    # track to prevent repeat annotation if same sample used as comparison

    ana_dir = project_dir + project + '/' + annotation + '/' + sample_pair + '/OUTPUT'
    loc = '../LOGS/' + sample_pair + '.vep91_anno.log'
    cwd = os.getcwd()
    os.chdir(ana_dir)
    log(loc, date_time() + 'Changed to annotation directory ' + ana_dir + '\n')
    in_vcf = ana_dir + '/' + sample_pair + in_suffix
    out_vcf = sample_pair + out_suffix
    # run_vep = ''
    buffer_size = '5000'
    run_cmd = run_vep(vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc,
                      plugin_dir)
    log(loc, date_time() + 'Annotating sample ' + sample_pair + in_suffix + ' ' + run_cmd + '\n')
    # from stack overflow to allow killing of spawned processes in main process fails for cleaner restart
    check = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    check_run = watch_mem(check, source, sample_pair, loc)
    if check_run != 0:

        buffer_size = str(int(buffer_size)/2)
        clean_up = 'rm ' + out_vcf + '*'
        log(loc, date_time() + 'VEP failed. Status of run was ' + str(check_run) + ' Trying smaller buffer size of '
            + buffer_size + '\n' + clean_up + '\n')
        os.killpg(os.getpgid(check.pid), signal.SIGINT)

        subprocess.call(clean_up, shell=True)
        run_cmd = run_vep(vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc,
                          plugin_dir)
        log(loc, date_time() + 'Annotating sample ' + sample_pair + in_suffix + '\n')
        check = subprocess.call(run_cmd, shell=True)
        if check != 0:
            log(loc, date_time() + 'VEP failed for sample ' + sample_pair + '\n')
            exit(1)
    else:
        log(loc, date_time() + 'VEP annotation of ' + sample_pair + in_suffix + ' successful!\n')

    check = gen_strelka_report(out_vcf, intvl, tx_index)
    if check != 0:
        log(loc, date_time() + 'Report generation for ' + out_vcf + ' failed\n')
        exit(1)
    log(loc, date_time() + 'Changing back to working directory ' + cwd + '\n')
    os.chdir(cwd)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pair', action='store', dest='sample_pair', help='Sample tumor/normal pair')
    parser.add_argument('-is', '--in_suffix', action='store', dest='in_suffix', help='Suffix of input files')
    parser.add_argument('-os', '--out_suffix', action='store', dest='out_suffix', help='Suffix of output vcf files')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pair, in_suffix, out_suffix, in_mutect, source) = (inputs.config_file,
            inputs.sample_pair, inputs.in_suffix, inputs.out_suffix, inputs.in_mutect, inputs.source)
    annot_vcf_vep_pipe(config_file, sample_pair, in_suffix, out_suffix)
