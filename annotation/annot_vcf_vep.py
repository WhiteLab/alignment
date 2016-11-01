#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.job_manager import job_manager
from utility.log import log
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['fa_ordered'],
            config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['params']['vep_cache_version'],
            config_data['params']['threads'], config_data['params']['pass_flag'])


def pass_filter(sample, in_suffix, pflag):
    in_fn = sample + in_suffix
    pass_val = 'PASS'
    out_suffix = '.somatic_indel.PASS.vcf'
    if pflag == 'sub':
        out_suffix = '.keep.vcf'
        pass_val = 'KEEP'
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


def annot_vcf_vep_pipe(config_file, sample_pairs, ref_mnt, in_suffix, out_suffix):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, threads, pflag) = parse_config(config_file)
    fasta = ref_mnt + '/' + fasta
    vep_cache = ref_mnt + '/' + vep_cache
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
        if pflag != 'N':
            pass_filter(sample, in_suffix, pflag)
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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-is', '--in_suffix', action='store', dest='in_suffix', help='Suffix of input files')
    parser.add_argument('-os', '--out_suffix', action='store', dest='out_suffix', help='Suffix of output files')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt, in_suffix, out_suffix) = (inputs.config_file, inputs.sample_pairs,
                                                                  inputs.ref_mnt, inputs.in_suffix, inputs.out_suffix)
    annot_vcf_vep_pipe(config_file, sample_pairs, ref_mnt, in_suffix, out_suffix)
