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
            config_file['params']['threads'])


def mutect_annot_vep_pipe(config_file, sample_pairs, ref_mnt):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, threads) = parse_config(config_file)
    fasta = ref_mnt + '/' + fasta
    dbsnp = ref_mnt + '/' + dbsnp
    vcache = ref_mnt + '/' + vcache
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
        in_vcf = sample + '.out.keep.vcf'
        out_vcf = sample + '.out.keep.vep.vcf'
        run_vep = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf + ' --symbol --vcf_info_field' \
                ' ANN --canonical --html --variant_class --no_whole_genome --offline --maf_exac --no_whole_genome ' \
                '--fork ' + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache + '--cache_version ' + vcache \
                + ' --plugin dbNSFP,' + dbsnp + ' 2>> LOGS/' + sample + '.vep.log >> LOGS/' + sample + '.vep.log;'
        log(loc, date_time() + 'Annotating sample ' + sample + '\n')
        check = subprocess.call(run_vep, shell=True)
        if check != 0:
            log(loc, date_time() + 'VEP annotation for ' + sample + ' failed\n')
            exit(1)
        else:
            log(loc, date_time() + 'VEP annotation ' + sample + ' successful!\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt) = (inputs.config_file, inputs.sample_pairs, inputs.ref_mnt)
    mutect_annot_vep_pipe(config_file, sample_pairs, ref_mnt)
