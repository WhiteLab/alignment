#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from analysis.variant_annot_pipe import *


def get_out_vcf(cont, obj, pairs):
    src_cmd = '. /home/ubuntu/.novarc;'
    for pair in open(pairs):
        pair = pair.rstrip('\n')
        dirname = obj + '/' + pair + '/OUTPUT/'
        out = pair + '.out.keep'
        vcf = pair + '.vcf.keep'
        # output to cur dir to make vep work
        cmd = src_cmd + 'swift download ' + cont + ' ' + dirname + out + ' --skip-identical  --output ' + out \
              + '; swift download ' + cont + ' ' + dirname + vcf + ' --skip-identical --output ' + vcf
        sys.stderr.write(date_time() + cmd + '\n')
        check = call(cmd, shell=True)
        if check != 0:
            return 1, pair
    return 0, 'OK'


def temp_annot_pipe(config_file, sample_pairs, ref_mnt):
    (novosort, obj, cont, analysis, annotation, germ_flag, indel_flag, annot_used) = parse_config(config_file)
    sys.stderr.write(date_time() + 'Downloading mutect out files and vcf files\n')
    check = get_out_vcf(cont, analysis, sample_pairs)
    if check[0] != 0:
        sys.stderr.write(date_time() + 'Could not get vcfs for ' + check[1] + '\n')
        exit(1)
    if annot_used == 'vep':
        vep(config_file, sample_pairs, ref_mnt, '.vcf.keep', '.snv.vep.vcf', 'mutect')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mini temp pipe to run a part of variant_annot_pipe')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-k', '--karyo', action='store', dest='kflag',
                        help='Flag to perform karyotypic reordering of BAM files.  Only need if original reference used'
                             ' wasn\'t sorted in the manner. \'y\' to do so')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, kflag, ref_mnt) = (
        inputs.sample_pairs, inputs.config_file, inputs.kflag, inputs.ref_mnt)
    temp_annot_pipe(config_file, sample_pairs, ref_mnt)