#!/usr/bin/env python

import json
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
from utility.date_time import date_time
from subprocess import call
from alignment.novosort_merge_pe import novosort_merge_pe
from alignment.picard_ksort import ksort
from alignment.get_merged_bams import get_merged_bams
from mutect_pipe import mutect_pipe
from mutect_merge_sort import mutect_merge_sort
from annotation.snpeff_pipe import snpeff_pipe
from utility.upload_variants_to_swift import upload_variants_to_swift
from annotation.annot_scalpel_snpeff import annot_scalpel
from annotation.snpeff_scalpel_vcf2table import convert_vcf
from annotation.annot_platypus_VEP import annot_platypus
from platypus_germline import platypus_germline
from scalpel_indel import scalpel_indel


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['novosort'], config_data['refs']['obj'], config_data['refs']['cont'],
            config_data['refs']['analysis'], config_data['refs']['annotation'], config_data['params']['germflag'],
            config_data['params']['indelflag'])


def check_existing_bams(sample_list):
    temp_list = []
    test = open(sample_list, 'r')
    for fn in test:
        fn = fn.rstrip('\n')
        if not os.path.isfile(fn + '.merged.final.bam'):
            temp_list.append(fn)
    test.close()
    return temp_list


def variant_annot_pipe(config_file, sample_pairs, kflag, ref_mnt, wg, sm):
    # create eventual output location directories

    mk_dir = 'mkdir BAM LOGS ANALYSIS ANNOTATION'
    call(mk_dir, shell=True)
    (novosort, obj, cont, analysis, annotation, germ_flag, indel_flag) = parse_config(config_file)
    # create sample list
    sample_list = 'sample_list.txt'
    fh = open(sample_pairs, 'r')
    sl = open(sample_list, 'w')
    temp = {}
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        if cur[1] not in temp:
            sl.write(cur[1] + '\n')
            temp[cur[1]] = 1
        if cur[2] not in temp:
            sl.write(cur[2] + '\n')
            temp[cur[2]] = 1
    sl.close()
    fh .close()
    del temp
    # download and merge (if necessary) bam files
    if sm == 'n':
        check = novosort_merge_pe(config_file, sample_list)
        if check == 0:
            sys.stderr.write(date_time() + 'File download and merge complete!\n')
            # rm unmerged bams, no longer needed
            rm_bam = 'rm -rf ' + obj
            call(rm_bam, shell=True)
        else:
            sys.stderr.write(date_time() + 'File download and merge failed.\n')
            exit(1)
        if kflag == 'y':
            # create bam list for ksort
            bam_list = 'bam_list.txt'
            blist_cmd = 'ls *.merged.final.bam > ' + bam_list
            call(blist_cmd, shell=True)
            check = ksort(config_file, bam_list, kflag, ref_mnt)
            if check == 0:
                sys.stderr.write(date_time() + 'Karyotypic reorder of BAM files completed\n')
            else:
                sys.stderr.write(date_time() + 'Karyotypic reorder of BAM files failed.\n')
                exit(1)
    else:
        # quick check to see if just need to restart pipleine from mutect, or actually get merged bams
        sys.stderr.write(date_time() + 'Skip merge indicated, checking for already merged bam files\n')
        temp_list = check_existing_bams(sample_list)
        if len(temp_list) > 0:
            sys.stderr.write(date_time() + 'Missing files detected, downloading merged bam files\n')
            temp_fn = 'temp_samp_list.txt'
            temp_fh = open(temp_fn, 'w')
            temp_fh.write('\n'.join(temp_list))
            temp_fh.close()
            check = get_merged_bams(config_file, temp_fn)
            if check == 0:
                sys.stderr.write(date_time() + 'Merged bam files successfully download\n')
            else:
                sys.stderr.write(date_time() + 'Merged bam file download failed.  Make sure container, object correctly'
                                               ' configured in config file\n')
                exit(1)
        else:
            sys.stderr.write(date_time() + 'All bams found.  Moving on\n')
    check = mutect_pipe(config_file, sample_pairs, ref_mnt)
    if check == 0:
        sys.stderr.write(date_time() + 'Mutect variant calls successful\n')
    else:
        sys.stderr.write(date_time() + 'Mutect variant calls failed.\n')
        exit(1)
    check = mutect_merge_sort(config_file, sample_pairs, ref_mnt)
    if check == 0:
        sys.stderr.write(date_time() + 'Mutect file merge successful.\n')
    else:
        sys.stderr.write(date_time() + 'Mutect file merge failed.\n')

    check = snpeff_pipe(config_file, sample_pairs, ref_mnt, wg)
    if check == 0:
        sys.stderr.write(date_time() + 'snpEff successful.\n')
    else:
        sys.stderr.write(date_time() + 'snpEff failed.\n')
        exit(1)
    check = scalpel_indel(sample_pairs, 'LOGS/', config_file, ref_mnt)
    if check == 0:
        sys.stderr.write(date_time() + 'scalpel successful.\n')
    else:
        sys.stderr.write(date_time() + 'scalpel failed.\n')
        exit(1)
    check = annot_scalpel(config_file, sample_pairs, ref_mnt)
    if check == 0:
        sys.stderr.write(date_time() + 'annot scalpel successful.\n')
    else:
        sys.stderr.write(date_time() + 'annot scalpel failed.\n')
        exit(1)
    indel_vcf_suffix = '.somatic_indel.PASS.eff.vcf'
    check = convert_vcf(config_file, sample_pairs, indel_vcf_suffix)
    if check == 0:
        sys.stderr.write(date_time() + 'scalpel vcf2table successful.\n')
    else:
        sys.stderr.write(date_time() + 'scalpel vcf2table failed.\n')
        exit(1)

    if germ_flag == 'Y':
        sys.stderr.write(date_time() + 'Germ line call flag indicated\n')
        check = platypus_germline(config_file, sample_pairs, 'LOGS/', wg, ref_mnt)
        check += annot_platypus(config_file, sample_pairs, ref_mnt)
        if check == 0:
            sys.stderr.write(date_time() + 'Germ line call complete\n')
        else:
            sys.stderr.write(date_time() + 'Error during germline calls.  Check output\n')
            exit(1)

    # relocate stuff, then upload
    mv_cmds = 'rm -rf outdir; mv *.bai *.bam BAM; mv *.xls *eff* *sift* *vep* ANNOTATION; mv *out* *vcf* ANALYSIS;'
    call(mv_cmds, shell=True)
    check = upload_variants_to_swift(cont, obj, sample_list, sample_pairs, analysis, annotation)
    if check == 0:
        sys.stderr.write(date_time() + 'Uploading data to swift successful!\n')
    else:
        sys.stderr.write(date_time() + 'Uploading data to swift failed!\n')
        exit(1)

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-k', '--karyo', action='store', dest='kflag',
                        help='Flag to perform karyotypic reordering of BAM files.  Only need if original reference used'
                             ' wasn\'t sorted in the manner. \'y\' to do so')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-wg', '--whole-genome', action='store', dest='wg',
                        help='\'y\' or \'n\' flag if whole genome or not.  will determine whether to flag for on/off '
                             'target')
    parser.add_argument('-sm', '--skip-merge', action='store', dest='sm',
                        help='\'y\' or \'n\' flag to skip merge files.  Useful for repeating variant calls when BAMs'
                             ' were already merged, sorted, etc the first time')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, kflag, ref_mnt, wg, sm) = (
        inputs.sample_pairs, inputs.config_file, inputs.kflag, inputs.ref_mnt, inputs.wg, inputs.sm)
    variant_annot_pipe(config_file, sample_pairs, kflag, ref_mnt, wg, sm)
