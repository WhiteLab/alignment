#!/usr/bin/env python

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.job_manager import job_manager
from utility.source_novarc import source_novarc
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['snpEff'], config_data['tools']['snpsift'],
            config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['refs']['capture'],
            config_data['refs']['cont'], config_data['refs']['obj'], config_data['params']['threads'])


def snpeff_pipe(config_file, sample_list, ref_mnt, novarc):
    (java, snpeff, snpsift, report, dbsnp, bed, cont, obj, max_t) = parse_config(config_file)
    dbsnp = ref_mnt + '/' + dbsnp
    bed = ref_mnt + '/' + bed
    fh = open(sample_list)
    mk_log_dir = 'mkdir LOGS'
    subprocess.call(mk_log_dir, shell=True)
    cmd_list = []
    run_snpsift = java + ' -jar ' + snpsift + ' annotate ' + dbsnp
    run_snpeff = java + ' -jar ' + snpeff + ' eff -t hg19 -interval ' + bed
    source_novarc(novarc)
    for line in fh:
        line = line.rstrip('\n')
        in_vcf = obj + '/' + line + '.merged.final.bam.germline_calls.vcf'
        sift_vcf = obj + '/' + line + '.snpSift.vcf'
        final_vcf = obj + '/' + line + '.snpSift.snpEff.vcf'
        dl_vcf = 'swift download ' + cont + ' ' + in_vcf + ';'
        log = 'LOGS/' + line + '.snpeff.log'
        run_cmd = dl_vcf + run_snpsift + ' ' + in_vcf + ' > ' + sift_vcf + ' 2> ' + log + ';' + run_snpeff + ' ' \
                  + sift_vcf + ' -v > ' + final_vcf + ' 2>> ' + log
        cmd_list.append(run_cmd)
    job_manager(cmd_list, max_t)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Normal sample list')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-n', '--novarc', action='store', dest='novarc', help='novarc file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_list, ref_mnt, novarc) = (
        inputs.config_file, inputs.sample_list, inputs.ref_mnt, inputs.novarc)
    snpeff_pipe(config_file, sample_list, ref_mnt, novarc)
