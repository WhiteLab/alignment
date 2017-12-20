#!/usr/bin/env python3

import json
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import os
from utility.date_time import date_time
from subprocess import call
from analysis.mutect_pipe import mutect_pipe
from analysis.mutect_merge_sort import mutect_merge_sort
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
    # Working directory is sample analysis directory
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
    # organize analysis files
    reorg = 'mv *.log LOGS; rmdir outdir; find . -maxdepth 1 -type f -exec mv {} OUTPUT \;'
    sys.stderr.write('Reorganizing analysis files ' + reorg + '\n')
    call(reorg, shell=True)
    # Working directory now annotation directory
    ann_dir = project_dir + project + '/' + annotation + '/' + sample_pair
    if not os.path.isdir(ann_dir):
        mk_ann = 'mkdir -p ' + ann_dir + ' ' + ann_dir + '/LOGS ' + ann_dir + '/OUTPUT'
        sys.stderr.write('Creating annotation output directories ' + mk_ann + '\n')
        call(mk_ann, shell=True)

    os.chdir(ann_dir)
    if annot_used == 'vep':
        check = vep(config_file, sample_pair, '.vcf.keep', '.snv.vep.vcf', '.out.keep', 'mutect')
        check += vep(config_file, sample_pair, '.indel.vcf', '.somatic.indel.vep.vcf', 'NA', 'scalpel')

    if check == 0:
        sys.stderr.write(date_time() + 'File annotation successful, reorganizing files\n')
    else:
        sys.stderr.write(date_time() + 'File annotation failed! Check logs!\n')
    # organize annotation files
    reorg = 'mv *.log LOGS; find . -maxdepth 1 -type f -exec mv {} OUTPUT \;'
    sys.stderr.write(reorg + '\n')
    call(reorg, shell=True)

    if germ_flag == 'Y':
        sys.stderr.write(date_time() + 'Germ line call flag indicated\n')
        # check for germline analysis dir
        germ_ana_dir = project_dir + project + '/' + analysis + '/' + normal_id
        if not os.path.isdir(germ_ana_dir):
            mk_ana = 'mkdir -p ' + germ_ana_dir
            sys.stderr.write('Creating analysis output directories ' + mk_ana + '\n')
            call(mk_ana, shell=True)
        os.chdir(germ_ana_dir)
        check = platypus_germline(config_file, normal_id, ana_dir + '/LOGS/', wg)
        reorg = 'mv *.log ' + ana_dir + '/LOGS;'
        sys.stderr.write('Reorganizing germline analysis files ' + reorg + '\n')
        call(reorg, shell=True)
        # check for germline annotation dir
        germ_ann_dir = project_dir + project + '/' + annotation + '/' + normal_id
        if not os.path.isdir(germ_ann_dir):
            mk_ann = 'mkdir -p ' + germ_ann_dir
            sys.stderr.write('Creating analysis output directories ' + mk_ann + '\n')
            call(mk_ann, shell=True)
        os.chdir(germ_ann_dir)
        check += annot_platypus(config_file, normal_id)
        reorg = 'mv *.log ' + ann_dir + '/LOGS;'
        sys.stderr.write('Reorganizing germline analysis files ' + reorg + '\n')
        call(reorg, shell=True)
        if check == 0:
            sys.stderr.write(date_time() + 'Germ line call complete\n')
        else:
            sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
            exit(1)

    if check == 0:
        sys.stderr.write(date_time() + 'Variant calling pipeline complete.  Check outputs\n')
    else:
        sys.stderr.write(date_time() + 'Variant calling pipeline FAILED.  Check outputs!\n')
        exit(1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-t', '--tumor', action='store', dest='tumor',
                        help='Tumor id')
    parser.add_argument('-n', '--normal', action='store', dest='normal',
                        help='Normal id')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations - use FULL PATH!')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (tumor, normal, config_file) = (inputs.tumor, inputs.normal, inputs.config_file)
    variant_annot_pipe(tumor, normal, config_file)
