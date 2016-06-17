#!/usr/bin/python
import json
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
sys.path.append('/home/ubuntu/TOOLS/Scripts/annotation')
from date_time import date_time
from subprocess import call
from novosort_merge_pe import novosort_merge_pe
from picard_ksort import ksort
from get_merged_bams import get_merged_bams
from mutect_pipe import mutect_pipe
from mutect_merge_sort import mutect_merge_sort
from snpeff_pipe import snpeff_pipe
from upload_variants_to_swift import upload_variants_to_swift


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['novosort'], config_data['refs']['obj'], config_data['refs']['cont'],
            config_data['refs']['analysis'], config_data['refs']['annotation'])


def variant_annot_pipe(config_file, sample_pairs, wait, kflag, ref_mnt, wg, sm):
    # create eventual output location directories

    mk_dir = 'mkdir BAM LOGS ANALYSIS ANNOTATION'
    call(mk_dir, shell=True)
    (novosort, obj, cont, analysis, annotation) = parse_config(config_file)
    # create sample list
    samp_cmd = 'cut -f 2 ' + sample_pairs + ' > sample_list.txt;' + 'cut -f 3 ' + sample_pairs + ' >> sample_list.txt'
    call(samp_cmd, shell=True)
    sample_list = 'sample_list.txt'
    # download and merge (if necessary) bam files
    if sm == 'n':
        check = novosort_merge_pe(config_file, sample_list, wait)
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
        sys.stderr.write(date_time() + 'Skip merge indicated, downloading merged bam files\n')
        check = get_merged_bams(config_file, sample_list, wait)
        if check == 0:
            sys.stderr.write(date_time() + 'Merged bam files successfully download\n')
        else:
            sys.stderr.write(
                date_time() + 'Merged bam file download failed.  Make sure container, object correctly configured in config file\n')
            exit(1)
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
    # relocate stuff, then upload
    mv_cmds = 'mv *.bai *.bam BAM;mv *eff* *sift* ANNOTATION; mv *out* *vcf* ANALYSIS'
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
    parser.add_argument('-w', '--wait', action='store', dest='wait',
                        help='Wait time to download bam files.  900 (seconds) recommended')
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
    (sample_pairs, config_file, wait, kflag, ref_mnt, wg, sm) = (
        inputs.sample_pairs, inputs.config_file, inputs.wait, inputs.kflag, inputs.ref_mnt, inputs.wg, inputs.sm)
    variant_annot_pipe(config_file, sample_pairs, wait, kflag, ref_mnt, wg, sm)
