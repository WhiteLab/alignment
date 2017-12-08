#!/usr/bin/env python

import glob
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import os
from subprocess import call
from date_time import date_time


def upload_variants_to_swift(cont, obj, sample_list, sample_pairs, analysis, annotation, annot_used):
    src_cmd = '. ~/.novarc;'
    ONE_GB = 1073741824
    fh = open(sample_list, 'r')
    for sample in fh:
        sample = sample.rstrip('\n')
        swift_cmd = src_cmd + 'swift upload ' + cont + ' BAM/' + sample + '.merged.final.bam -S ' + str(ONE_GB)\
                    + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample \
                    + '.merged.final.bam >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
        check = call(swift_cmd, shell=True)
        # depending on software used to index, .bai extension may follow bam
        bai = 'BAM/' + sample + '.merged.final.bai'
        if not os.path.isfile(bai):
            bai = 'BAM/' + sample + '.merged.final.bam.bai'
        swift_cmd = src_cmd + 'swift upload ' + cont + ' ' + bai + '  -S ' + str(ONE_GB)\
                    + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample\
                    + '.merged.final.bai >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
        try:
            call(swift_cmd, shell=True)
        except:
            swift_cmd = src_cmd + 'swift upload ' + cont + ' BAM/' + sample + '.merged.final.bam.bai -S ' + str(ONE_GB)\
                        + ' --skip-identical --object-name ' + obj + '/' + sample + '/BAM/' + sample\
                        + '.merged.final.bai >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
            check += call(swift_cmd, shell=True)
        if check == 0:
            sys.stderr.write(date_time() + 'Uploading final BAMs for ' + sample + ' successful!\n')
        else:
            sys.stderr.write(date_time() + 'Uploading final BAMs for ' + sample + ' failed\n')
            exit(1)
        # check for germline call of file, skip if not present
        germ = 'ANALYSIS/' + sample + '.germline_calls.vcf'
        if os.path.isfile(germ):
            suffix_list = ['.germline_calls.vcf', '.germline_pass.vcf']

            for suffix in suffix_list:
                swift_cmd = src_cmd + 'swift upload ' + cont + ' ANALYSIS/' + sample + suffix + ' -S ' + str(ONE_GB)\
                            + ' --skip-identical --object-name ' + analysis + '/' + sample + '/' + sample\
                            + suffix + ' >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
                check += call(swift_cmd, shell=True)
            suffix_list = ['.germline_pass.vep.vcf', '.germline_pass.xls',
                           '.germline_pass.vep.vcf.html', '.germline_pass.vep.vcf_summary.html']
            for suffix in suffix_list:
                swift_cmd = src_cmd + 'swift upload ' + cont + ' ANNOTATION/' + sample + suffix + ' -S ' + str(ONE_GB)\
                            + ' --skip-identical --object-name ' + annotation + '/' + sample + '/' + sample\
                            + suffix + ' >> LOGS/' + sample + '.upload.log 2>> LOGS/' + sample + '.upload.log'
                check += call(swift_cmd, shell=True)
            if check == 0:
                sys.stderr.write(date_time() + 'Uploading germline calls for ' + sample + ' successful!\n')
            else:
                sys.stderr.write(date_time() + 'Uploading germline calls for ' + sample + ' failed\n')
                exit(1)
    fh.close()

    # upload analysis and annotation files
    suffix_list1 = ['.out', '.out.keep', '.vcf', '.vcf.keep']
    suffix_list2 = ['.out.keep.eff.vcf', '.vcf.keep.eff.xls']
    if annot_used == 'vep':
        suffix_list2 = ['.snv.vep.vcf', '.subsitutions.vep.prioritized_impact.report.xls', '.snv.vep.vcf_summary.html']
    fh = open(sample_pairs, 'r')
    for line in fh:
        pair = line.rstrip('\n').split('\t')[0]
        # upload analysis files
        for suffix in suffix_list1:
            swift_cmd = src_cmd + 'swift upload ' + cont + ' ANALYSIS/' + pair + suffix + ' -S ' + str(
                ONE_GB) + ' --skip-identical --object-name ' + analysis + '/' + pair + '/OUTPUT/' + pair + suffix\
                        + ' >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
            check = call(swift_cmd, shell=True)
            if check == 0:
                sys.stderr.write(date_time() + 'Uploading analysis file ' + pair + suffix + ' successful!\n')
            else:
                sys.stderr.write(date_time() + 'Uploading analysis file ' + pair + suffix + ' failed!\n')
                exit(1)
        # upload annotation files
        for suffix in suffix_list2:
            swift_cmd = src_cmd + 'swift upload ' + cont + ' ANNOTATION/' + pair + suffix + ' -S ' + str(
                ONE_GB) + ' --skip-identical --object-name ' + annotation + '/' + pair + '/OUTPUT/' + pair + suffix\
                        + ' >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
            check = call(swift_cmd, shell=True)
            if check == 0:
                sys.stderr.write(date_time() + 'Uploading annotation file ' + pair + suffix + ' successful!\n')
            else:
                sys.stderr.write(date_time() + 'Uploading annotation file ' + pair + suffix + ' failed!\n')
                exit(1)
        # upload log files
        mut_list = glob.glob('LOGS/' + pair + '.mut*')
        mut_list.extend(glob.glob('LOGS/' + pair + '.scalpel*'))
        for mut in mut_list:
            swift_cmd = src_cmd + 'swift upload ' + cont + ' ' + mut + ' -S ' + str(
                ONE_GB) + ' --skip-identical --object-name ' + analysis + '/' + pair + '/' + mut + ' >> LOGS/'\
                        + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
            check = call(swift_cmd, shell=True)
            if check == 0:
                sys.stderr.write(date_time() + 'Uploading analysis log file ' + mut + ' successful!\n')
            else:
                sys.stderr.write(date_time() + 'Uploading analysis log file ' + mut + ' failed!\n')
                exit(1)
        swift_cmd = src_cmd + 'swift upload ' + cont + ' LOGS/' + pair + '.snpeff.log -S ' + str(
            ONE_GB) + ' --skip-identical --object-name ' + annotation + '/' + pair + '/LOGS/' + pair\
                    + '.snpeff.log >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
        if annot_used == 'vep':
            swift_cmd = src_cmd + 'swift upload ' + cont + ' LOGS/' + pair + '.vep_anno.log -S ' + str(
            ONE_GB) + ' --skip-identical --object-name ' + annotation + '/' + pair + '/LOGS/' + pair\
                    + '.vep_anno.log >> LOGS/' + pair + '.upload.log 2>> LOGS/' + pair + '.upload.log'
        check = call(swift_cmd, shell=True)
        if check == 0:
            sys.stderr.write(date_time() + 'Uploading annotation log file ' + pair + ' using ' + annot_used
                             + ' successful!\n')
        else:
            sys.stderr.write(date_time() + 'Uploading annotation log file ' + pair + ' using ' + annot_used
                             + ' failed!\n')
            exit(1)

        # check for indel call files, upload in present
        indel_vcf = pair + '/' + pair + '.somatic_indel.PASS.eff.vcf'
        if annot_used == 'vep':
            indel_vcf  = pair + '/somatic.indel.vcf'
        if os.path.isfile(indel_vcf):
            if annot_used == 'vep':
                ana_list = glob.glob(pair + '/*indel*')
                # need to rename variats subfile
                ana_list.append(indel_vcf)
                ann_list = glob.glob('ANNOTATION/' + pair + '*indel*')
            else:
                ana_list = glob.glob(pair + '/*PASS.sift*')
                ana_list.extend(glob.glob(pair + '/*PASS.vcf'))
                ana_list.extend(glob.glob(pair + '/*.indel.vcf'))
                ann_list = (pair + '/' + pair + '.indels.xls', indel_vcf)
            for ana in ana_list:

                fn = os.path.basename(ana)
                swift_cmd = src_cmd + 'swift upload ' + cont + ' ' + ana + ' --object-name ' + analysis + '/' + pair \
                            + '/OUTPUT/' + fn
                check = call(swift_cmd, shell=True)
                if check == 0:
                    sys.stderr.write(date_time() + 'Uploading analysis vcf file ' + fn + ' successful!\n')
                else:
                    sys.stderr.write(date_time() + 'Uploading analysis vcf file ' + fn + ' failed!\n')
                    exit(1)
            for ann in ann_list:
                fn = os.path.basename(ann)
                swift_cmd = src_cmd + 'swift upload ' + cont + ' ' + ann + ' --object-name ' + annotation + '/' + pair \
                            + '/OUTPUT/' + fn
                check = call(swift_cmd, shell=True)
                if check == 0:
                    sys.stderr.write(date_time() + 'Uploading annotation vcf file ' + fn + ' successful!\n')
                else:
                    sys.stderr.write(date_time() + 'Uploading annotation vcf file ' + fn + ' failed!\n')
                    exit(1)

    fh.close()

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uploads current directory contents to specified object and container')
    parser.add_argument('-c', '--container', action='store', dest='cont',
                        help='Swift container name to upload to.  i.e. PANCAN')
    parser.add_argument('-o', '--object', action='store', dest='obj',
                        help='Swift object name root to use for aligned merged bam files.  i.e. ALIGN/2015-1234')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample list, one per line')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs',
                        help='Sample tumor/normal pairs, tsv file with bid pair, sample1, sample2')
    parser.add_argument('-ana', '--analysis', action='store', dest='analysis',
                        help='Analysis output object root, i.e. ANALYSIS')
    parser.add_argument('-ann', '--annotation', action='store', dest='annotation',
                        help='Analysis output object root, i.e. ANNOTATION')
    parser.add_argument('-anu', '--annot_used', action='store', dest='annot_used',
                        help='Software used for annotation, likely snpeff or vep')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (cont, obj, sample_list, sample_pairs, analysis, annotation, annot_used) = (
        inputs.cont, inputs.obj, inputs.sample_list, inputs.sample_pairs, inputs.analysis, inputs.annotation,
        inputs.annot_used)
    upload_variants_to_swift(cont, obj, sample_list, sample_pairs, analysis, annotation, annot_used)
