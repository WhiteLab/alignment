#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import os
from utility.date_time import date_time
from subprocess import call
from analysis.mutect_pipe import mutect_pipe
from analysis.mutect_merge_sort import mutect_merge_sort
from utility.upload_variants_to_swift import upload_variants_to_swift
from annotation.annot_platypus_VEP import annot_platypus
from analysis.platypus_germline import platypus_germline
from analysis.scalpel_indel import scalpel_indel


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align'], \
           config_data['refs']['analysis'], config_data['refs']['annotation'], config_data['params']['germflag'], \
           config_data['params']['indelflag'], config_data['params']['annotator'], config_data['params']['wg_flag']


def vep(config_file, sample_pairs, in_suffix, out_suffix, in_mutect, source):
    from annotation.annot_vcf_vep import annot_vcf_vep_pipe
    check = annot_vcf_vep_pipe(config_file, sample_pairs, in_suffix, out_suffix, in_mutect, source)
    if check == 0:
        sys.stderr.write(date_time() + 'vep annotation of ' + source + ' output successful.\n')
    else:
        sys.stderr.write(date_time() + 'vep annotation of ' + source + ' output failed.\n')
        exit(1)
    return 0


def variant_annot_pipe(tumor_id, normal_id, config_file):
    (project_dir, project, align, analysis, annotation, germ_flag, indel_flag, annot_used, wg) \
        = parse_config(config_file)
    sample_pair = tumor_id + '_' + normal_id
    ana_dir = project_dir + project + '/' + analysis + '/' + sample_pair
    if not os.path.isdir(ana_dir):
        mk_ana = 'mkdir -p ' + ana_dir + ' ' + ana_dir + '/LOGS ' + ana_dir + '/OUTPUT'
        sys.stderr.write('Creating anaylsis output directories ' + mk_ana + '\n')
        call(mk_ana, shell=True)
    os.chdir(ana_dir)
    check = mutect_pipe(config_file, tumor_id, normal_id)
    if check == 0:
        sys.stderr.write(date_time() + 'Mutect variant calls successful\n')
    else:
        sys.stderr.write(date_time() + 'Mutect variant calls failed.\n')
        exit(1)
    check = mutect_merge_sort(config_file, sample_pair)
    if check == 0:
        sys.stderr.write(date_time() + 'Mutect file merge successful.\n')
    else:
        sys.stderr.write(date_time() + 'Mutect file merge failed.\n')
    # create def to do vep or snpeff mode for annotation
    check = scalpel_indel(tumor_id, normal_id, 'LOGS/', config_file)
    if check == 0:
        sys.stderr.write(date_time() + 'scalpel successful.\n')
    else:
        sys.stderr.write(date_time() + 'scalpel failed.\n')
        exit(1)

    if annot_used == 'vep':
        vep(config_file, sample_pair, '.vcf.keep', '.snv.vep.vcf', '.out.keep', 'mutect')
        vep(config_file, sample_pair, '.indel.vcf', '.somatic.indel.vep.vcf', 'NA', 'scalpel')

    if germ_flag == 'Y':
        sys.stderr.write(date_time() + 'Germ line call flag indicated\n')
        check = platypus_germline(config_file, sample_pair, 'LOGS/', wg, ref_mnt)
        check += annot_platypus(config_file, sample_pair, ref_mnt)
        if check == 0:
            sys.stderr.write(date_time() + 'Germ line call complete\n')
        else:
            sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
            exit(1)

    # relocate stuff, then upload
    mv_cmds = 'rm -rf outdir; mv *.bai *.bam BAM; mv *.xls *eff* *sift* *vep* ANNOTATION; mv *out* *vcf* ANALYSIS;'
    call(mv_cmds, shell=True)
    check = upload_variants_to_swift(cont, obj, sample_list, sample_pairs, analysis, annotation, annot_used)
    if check == 0:
        sys.stderr.write(date_time() + 'Uploading data to swift successful!\n')
    else:
        sys.stderr.write(date_time() + 'Uploading data to swift failed!\n')
        exit(1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-t', '--tumor', action='store', dest='tumor',
                        help='Tumor id')
    parser.add_argument('-n', '--normal', action='store', dest='normal',
                        help='Normal id')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (tumor, normal, config_file) = (inputs.tumor, inputs.normal, inputs.config_file)
    variant_annot_pipe(tumor, normal, config_file)
