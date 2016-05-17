#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['snpEff'], config_data['tools']['snpsift'],
            config_data['tools']['gatk'], config_data['refs']['dbsnp'], config_data['refs']['intervals'], 
            config_data['refs']['fa_ordered'])

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
    # edit to grab from config max thread count
    (java, snpeff, snpsift, gatk, dbsnp, intervals, fasta) = parse_config(config_file)
    # pass_filter(sample)
    dbsnp = ref_mnt + '/' + dbsnp
    mk_log_dir = 'mkdir LOGS'
    '''
    subprocess.call(mk_log_dir, shell=True)
    run_snpsift = java + ' -jar ' + snpsift + ' annotate ' + dbsnp
    run_snpeff = java + ' -jar ' + snpeff + ' eff -t hg19 '
    run_snp = run_snpsift + ' ' + sample + '.germline_pass.vcf > ' + sample + '.germline_pass.sift.vcf 2> LOGS/'\
              + sample + '.snpeff.log;' + run_snpeff + ' ' + sample + '.germline_pass.sift.vcf -v > ' + sample \
              + '.germline_pass.eff.vcf  2>> LOGS/' + sample + '.snpeff.log;'
    check = subprocess.call(run_snp, shell=True)
    if check == 0:
        sys.stderr.write(date_time() + 'SNP annotation of germline calls completed!\n')
    else:
        sys.stderr.write(date_time() + 'SNP annotation of germline calls for ' + sample + ' FAILED!\n')
        return 1
    '''
    # use GATK tool to convert vcf to table
    fasta = ref_mnt + '/' + fasta
    table_cmd = java + ' -jar ' + gatk + ' -T VariantsToTable -SMA -V ' + sample + '.germline_pass.eff.vcf -o ' + sample + \
                '.germline_pass.xls -F CHROM -F POS -F ID -F REF -F ALT -GF INFO.EFF.Effect -GF Codon_Change -GF Amino_Acid_Change -GF Amino_Acid_length -GF Gene_Name -GF Transcript_BioType -GF Gene_Coding -R ' + fasta
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
