#!/usr/bin/env python3

import sys
import os
import signal
from annotation.vep_subsitution_report import gen_report as gen_snv_report
from annotation.vep_indel_report import gen_report as gen_indel_report
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
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
           config_data['params']['dustmask_flag'], config_data['params']['wg_flag'], config_data['refs']['tx_index'],\
           config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['analysis']


def pass_filter(ana_dir, sample, in_suffix, dustmask_flag):
    in_fn = ana_dir + '/somatic' + in_suffix
    if dustmask_flag == 'Y':
        in_fn = ana_dir + '/' + sample + '.somatic.indel.dustmasked.vcf'
    pass_val = 'PASS'
    out_suffix = '.somatic_indel.PASS.vcf'
    out_fn = sample + out_suffix
    out = open(out_fn, 'w')
    infile = open(in_fn, 'r')
    for line in infile:
        if line[0] == '#':
            out.write(line)
        else:
            fields = line.split('\t')
            if fields[6] == pass_val:
                out.write(line)
    infile.close()
    out.close()


def run_vep(wg_flag, vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc, plugin_dir):
    if wg_flag == 'n':
        run_cmd = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf \
                  + ' --symbol --vcf_info_field ANN --canonical --variant_class --buffer_size ' + buffer_size \
                  + ' --offline --maf_exac --no_whole_genome --fork ' + threads + ' --fasta ' \
                  + fasta + ' --dir_cache ' + vep_cache + ' --cache_version ' + vcache + ' --dir_plugins ' + plugin_dir\
                  + '  2>> ' + loc + ' >> ' \
                  + loc
    else:
        run_cmd = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf \
                  + ' --symbol --vcf_info_field ANN --canonical --variant_class --buffer_size ' + buffer_size \
                  + ' --offline --maf_exac --fork ' + threads + ' --fasta ' + fasta + \
                  ' --dir_cache ' + vep_cache + ' --cache_version ' + vcache + ' --dir_plugins ' + plugin_dir \
                  + ' 2>> ' + loc + ' >> ' + loc
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


def annot_vcf_vep_pipe(config_file, sample_pair, in_suffix, out_suffix, in_mutect, source):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, plugin_dir, threads, intvl, dustmask_flag, wg_flag, tx_index,
     project_dir, project, analysis) = parse_config(config_file)
    # scale back on the forking a bit

    if int(threads) > 2:
        # threads = str(int(threads)/2 - 1)
        threads = str(int(threads) - 1)
    # track to prevent repeat annotation if same sample used as comparison
    loc = 'LOGS/' + sample_pair + '.vep_anno.log'
    ana_dir = project_dir + project + '/' + analysis + '/' + sample_pair + '/OUTPUT'
    in_vcf = ana_dir + '/' + sample_pair + in_suffix
    out_vcf = sample_pair + out_suffix
    if source == 'scalpel':
        pass_filter(ana_dir, sample_pair, in_suffix, dustmask_flag)
        in_vcf = sample_pair + '.somatic_indel.PASS.vcf'
    # run_vep = ''
    buffer_size = '2000'
    run_cmd = run_vep(wg_flag, vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc,
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
        run_cmd = run_vep(wg_flag, vep_tool, in_vcf, out_vcf, buffer_size, threads, fasta, vep_cache, vcache, loc,
                          plugin_dir)
        log(loc, date_time() + 'Annotating sample ' + sample_pair + in_suffix + '\n')
        check = subprocess.call(run_cmd, shell=True)
        if check != 0:
            log(loc, date_time() + 'VEP failed for sample ' + sample_pair + '\n')
            exit(1)
    else:
        log(loc, date_time() + 'VEP annotation of ' + sample_pair + in_suffix + ' successful!\n')
    if source == 'mutect':
        if wg_flag == 'y':
            intvl = 'n'
        check = gen_snv_report(out_vcf, sample_pair + in_mutect, intvl, tx_index)
        if check != 0:
            log(loc, date_time() + 'Report generation for ' + out_vcf + ' failed\n')
            exit(1)
    else:
        check = gen_indel_report(out_vcf, tx_index)
        if check != 0:
            log(loc, date_time() + 'Report generation for ' + out_vcf + ' failed\n')
            exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pair', action='store', dest='sample_pair', help='Sample tumor/normal pair')
    parser.add_argument('-so', '--source', action='store', dest='source', help='Variant call annot source - mutect or'
                                                                               ' scalpel')
    parser.add_argument('-is', '--in_suffix', action='store', dest='in_suffix', help='Suffix of input files')
    parser.add_argument('-os', '--out_suffix', action='store', dest='out_suffix', help='Suffix of output vcf files')
    parser.add_argument('-im', '--in_mutect', action='store', dest='in_mutect', help='Suffix of input mutect files, '
                                                                                     'can be \'NA\' if not from mutect')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pair, in_suffix, out_suffix, in_mutect, source) = (inputs.config_file,
            inputs.sample_pair, inputs.in_suffix, inputs.out_suffix, inputs.in_mutect, inputs.source)
    annot_vcf_vep_pipe(config_file, sample_pair, in_suffix, out_suffix, in_mutect, source)
