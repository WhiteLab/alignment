#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['fa_ordered'],\
           config_data['params']['threads']


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


def annot_platypus(config_file, sample, ref_mnt):
    (vep_tool, vep_cache, fasta, threads) = parse_config(config_file)
    pass_filter(sample)
    vep_cache = ref_mnt + '/' + vep_cache
    mk_log_dir = 'mkdir LOGS'
    subprocess.call(mk_log_dir, shell=True)
    run_vep = 'perl ' + vep_tool + ' --cache -i ' + sample + '.germline_pass.vcf --vcf -o ' + sample\
              + '.germline_pass.vep.vcf --symbol --html --variant_class --offline --maf_exac --no_whole_genome --fork '\
              + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache + ' 2>> LOGS/' + sample + '.vep.log;'
    check = subprocess.call(run_vep, shell=True)
    if check == 0:
        sys.stderr.write(date_time() + 'SNP annotation of germline calls completed!\n')
    else:
        sys.stderr.write(date_time() + 'SNP annotation of germline calls for ' + sample + ' FAILED!\n')
        return 1

    table_cmd = java + ' -jar ' + snpsift + ' extractFields ' + sample + '.germline_pass.eff.vcf CHROM POS ID REF ALT '\
                '"EFF[0].EFFECT" "EFF[0].CODON" "EFF[0].AA" "EFF[0].AA_LEN" "EFF[0].GENE" ' \
                '"EFF[0].BIOTYPE" "EFF[0].CODING" > ' + sample + '.germline_pass.xls'
    check = subprocess.call(table_cmd, shell=True)
    if check == 0:
        sys.stderr.write(date_time() + 'Germline table created!\n')
        return 0
    else:
        sys.stderr.write(date_time() + 'Germline table failed!\n')

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
    annot_platypus(config_file, sample_pairs, ref_mnt)
