#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json
from job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['snpEff'], config_data['tools']['snpsift'],
            config_data['refs']['dbsnp'], config_data['refs']['intervals'], config_data['params']['threads'])

def pass_filter(sample):
    in_fn = sample + '/somatic.indel.vcf'
    out_fn = sample + '/' + sample + '.somatic_indel.PASS.vcf'
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


def annot_scalpel(config_file, sample_pairs, ref_mnt):
    (java, snpeff, snpsift, dbsnp, intervals, th) = parse_config(config_file)
    job_list = []
    dbsnp = ref_mnt + '/' + dbsnp
    for line in open(sample_pairs, 'r'):
        sample = line.rstrip('\n')
        pass_filter(sample)
        out_fn = sample + '/' + sample + '.somatic_indel.PASS.vcf'
        out_fn1 = sample + '/' + sample + '.somatic_indel.PASS.sift.vcf'
        out_fn2 = sample + '/' + sample + '.somatic_indel.PASS.eff.vcf'
        mk_log_dir = 'mkdir LOGS'
        subprocess.call(mk_log_dir, shell=True)
        run_snpsift = java + ' -jar ' + snpsift + ' annotate ' + dbsnp
        run_snpeff = java + ' -jar ' + snpeff + ' eff -t hg19 '
        run_snp = run_snpsift + ' ' + out_fn + ' > ' + out_fn1 + ' 2> LOGS/'\
                  + sample + '.snpeff.log;' + run_snpeff + ' ' + out_fn1 + ' -v > ' + out_fn2 \
                  + ' 2>> LOGS/' + sample + '.snpeff.log;'
        job_list.append(run_snp)
    job_manager(job_list, th)
    return 0

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
    (config_file, sample_pairs, ref_mnt) = (
        inputs.config_file, inputs.sample_pairs, inputs.ref_mnt)
    annot_scalpel(config_file, sample_pairs, ref_mnt)
