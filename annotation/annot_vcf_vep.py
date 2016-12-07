#!/usr/bin/env python

import sys
from vep_subsitution_report import gen_report as gen_snv_report
from vep_indel_report import gen_report as gen_indel_report
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.log import log
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['fa_ordered'],
            config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['params']['vep_cache_version'],
            config_data['params']['threads'], config_data['refs']['intervals'], config_data['params']['dustmask_flag'])


def pass_filter(sample, in_suffix, dustmask_flag):
    in_fn = sample + '/somatic' + in_suffix
    if dustmask_flag == 'Y':
        in_fn = sample + '/' + sample + '.somatic.indel.dustmasked.vcf'
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


def annot_vcf_vep_pipe(config_file, sample_pairs, ref_mnt, in_suffix, out_suffix, source):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, threads, intvl, dustmask_flag) = parse_config(config_file)
    fasta = ref_mnt + '/' + fasta
    vep_cache = ref_mnt + '/' + vep_cache
    intvl = ref_mnt + '/' + intvl
    # vep is a beast, so leaving it at 2 for now
    threads = '2'
    # scale back on the forking a bit
    if int(threads) > 2:
        threads = str(int(threads)/2)
    # parse sample file, use only last if pairs
    samp_fh = open(sample_pairs, 'r')
    # track to prevent repeat annotation if same sample used as comparison
    for line in samp_fh:
        info = line.rstrip('\n').split('\t')
        sample = info[0]
        mk_log_dir = 'mkdir LOGS'
        subprocess.call(mk_log_dir, shell=True)
        loc = 'LOGS/' + sample + '.vep_anno.log'
        in_vcf = sample + in_suffix
        out_vcf = sample + out_suffix
        if source == 'scalpel':
            pass_filter(sample, in_suffix, dustmask_flag)
            in_vcf = sample + '.somatic_indel.PASS.vcf'
        run_vep = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf + ' --symbol --vcf_info_field' \
                ' ANN --canonical --variant_class --no_whole_genome --offline --maf_exac --no_whole_genome ' \
                '--fork ' + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache + ' --cache_version ' + vcache \
                + ' 2>> ' + loc + ' >> ' + loc
        log(loc, date_time() + 'Annotating sample ' + sample + in_suffix + '\n')
        check = subprocess.call(run_vep, shell=True)
        if check != 0:
            log(loc, date_time() + 'VEP annotation for ' + sample + in_suffix + ' failed\n')
            exit(1)
        else:
            log(loc, date_time() + 'VEP annotation ' + sample + in_suffix + ' successful!\n')
        if source == 'mutect':
            check = gen_snv_report(out_vcf, sample + '.out.keep', intvl)
            if check != 0:
                log(loc, date_time() + 'Report generation for ' + out_vcf + ' failed\n')
                exit(1)
        else:
            check = gen_indel_report(out_vcf)
            if check != 0:
                log(loc, date_time() + 'Report generation for ' + out_vcf + ' failed\n')
                exit(1)
    return 0

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-so', '--source', action='store', dest='source', help='Variant call annot source - mutect or'
                                                                               ' scalpel')
    parser.add_argument('-is', '--in_suffix', action='store', dest='in_suffix', help='Suffix of input files')
    parser.add_argument('-os', '--out_suffix', action='store', dest='out_suffix', help='Suffix of output files')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt, in_suffix, out_suffix, source) = (inputs.config_file, inputs.sample_pairs,
                                                    inputs.ref_mnt, inputs.in_suffix, inputs.out_suffix, inputs.source)
    annot_vcf_vep_pipe(config_file, sample_pairs, ref_mnt, in_suffix, out_suffix, source)
