#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import json
from utility.date_time import date_time
from subprocess import call
from subprocess import check_output
from utility.log import log
from annotation.annot_strelka_VEP91 import annot_vcf_vep_pipe
from utility.set_acls import set_acls


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['manta_cfg'], config_data['tools']['strelka_cfg'], config_data['refs']['capture_bgzip'],\
           config_data['refs']['fa_ordered'], config_data['params']['threads'], config_data['params']['wg_flag'],\
           config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'],\
           config_data['refs']['analysis'], config_data['refs']['annotation'],config_data['params']['user'],\
           config_data['params']['group']


# def wg_mode(scalpel, tumor_bam, normal_bam, fasta, cpus, pair, config_file):
#     config_data = json.loads(open(config_file, 'r').read())
#     exome = config_data['refs']['exome']
#     loc = 'LOGS/' + pair + '_' + pair + '.genome_as_exome.strelka.log'
#     cmd = scalpel + ' --somatic --logs --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' \
#           + normal_bam + ' --window 600 --two-pass --bed ' + exome + ' --ref ' + fasta + ' 2> ' + loc
#     log(loc, date_time() + cmd + '\n')
#     check = call(cmd, shell=True)
#     if check != 0:
#         return 1, pair
#     return 0, pair


def run_strelka(tumor_id, normal_id, config_file):
    (manta_cfg, strelka_cfg, bed, fasta, cpus, wg, project_dir, project, align, analysis_dir, annotation_dir, user,
     group) = parse_config(config_file)

    sample_pair = tumor_id + '_' + normal_id
    run_dir_prefix = project_dir + project + '/' + analysis_dir + '/' + sample_pair
    anno_complete_path = project_dir + project + '/' + annotation_dir + '/' + sample_pair
    manta_dir = run_dir_prefix + '/manta_out'
    manta_run = manta_dir + '/runWorkflow.py'
    strelka_dir = run_dir_prefix + '/strelka_out'
    strelka_run = strelka_dir + '/runWorkflow.py'
    loc = run_dir_prefix + '/LOGS/' + sample_pair + '.strelka.log'
    bam_dir = project_dir + project + '/' + align
    tumor_bam = bam_dir + '/' + tumor_id + '/BAM/' + tumor_id + '.merged.final.bam'
    normal_bam = bam_dir + '/' + normal_id + '/BAM/' + normal_id + '.merged.final.bam'
    if wg == 'n':
        manta_setup_cmd = manta_cfg + ' --tumorBam ' + tumor_bam + ' --normalBam ' + normal_bam + ' --callRegions ' \
                          + bed + ' --exome --referenceFasta ' + fasta + ' --runDir ' + manta_dir + ' 2>> ' + loc
        sys.stderr.write(date_time() + 'Starting indel calls for ' + sample_pair + '\n')
        log(loc, date_time() + 'Starting indel calls for ' + sample_pair + ' in capture mode with command:\n'
            + manta_setup_cmd + '\n')
        check = call(manta_setup_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'Indel calling setup failed for pair ' + sample_pair + ' with command:\n' +
                             manta_setup_cmd + '\n')
            exit(1)
        manta_run_cmd = manta_run + ' -m local -j ' + cpus
        check = call(manta_run_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'Indel calling run failed for pair ' + sample_pair + ' with command:\n' +
                     manta_setup_cmd + '\n')
            exit(1)
            strelka_setup_cmd = strelka_cfg + ' --tumorBam ' + tumor_bam + ' --normalBam ' + normal_bam \
                                + ' --callRegions ' + bed + ' --exome --referenceFasta ' + fasta + ' --runDir ' \
                                + strelka_dir + ' --indelCandidates ' + manta_dir \
                                + '/results/variants/candidateSmallIndels.vcf.gz 2>> ' + loc
            sys.stderr.write(date_time() + 'Starting snv calls for ' + sample_pair + '\n')
            log(loc, date_time() + 'Starting indel calls for ' + sample_pair + ' in capture mode with command:\n'
                + strelka_setup_cmd + '\n')
            check = call(strelka_setup_cmd, shell=True)
            if check != 0:
                sys.stderr.write(date_time() + 'SNV calling setup failed for pair ' + sample_pair + ' with command:\n' +
                                 strelka_setup_cmd + '\n')
                exit(1)
            strelka_run_cmd = strelka_run + ' -m local -j ' + cpus
            check = call(strelka_run_cmd, shell=True)
            if check != 0:
                sys.stderr.write(date_time() + 'SNV calling run failed for pair ' + sample_pair + ' with command:\n' +
                                 strelka_setup_cmd + '\n')
                exit(1)
    # else:
    #     check = wg_mode(scalpel, tumor_bam, normal_bam, fasta, cpus, sample_pair, config_file)
    #     if check[0] != 0:
    #         sys.stderr.write('Scalpel failed for ' + normal_id + ' at ' + tumor_id + '\n')
    #         exit(1)
    log(loc, date_time() + 'Variant calling complete for pair ' + sample_pair + ' filtering output files\n')
    strelka_snv_vcf = strelka_dir + '/results/variants/somatic.snvs.vcf.gz'
    strelka_snv_pass = run_dir_prefix + '/' + sample_pair + '.strelka.snv_PASS.vcf'
    strelka_indel_vcf = strelka_dir + '/results/variants/somatic.indels.vcf.gz'
    strelka_indel_pass = run_dir_prefix + '/' + sample_pair + '.strelka.indel_PASS.vcf'
    filter_vcf_cmd = 'zcat ' + strelka_snv_vcf + ' | grep -E "^#|PASS" > ' + strelka_snv_pass + ';'
    filter_vcf_cmd += 'zcat ' + strelka_indel_vcf + ' | grep -E "^#|PASS" > ' + strelka_indel_pass
    log(loc, date_time() + filter_vcf_cmd + '\n')
    call(filter_vcf_cmd, shell=True)
    sys.stderr.write(date_time() + 'Completed variant calls for ' + sample_pair + '\n')
    # check pass vcf before calling vep
    check_snv = 'cat ' + strelka_snv_pass + ' grep PASS | wc -l'
    check_snv = check_output(check_snv, shell=True)
    if check_snv.decode().rstrip('\n') == '0':
        sys.stderr.write(date_time() + 'Starting vep strelka snv annotation for ' + sample_pair + '\n')
        annot_vcf_vep_pipe(config_file, sample_pair, '.strelka.snv_PASS.vcf', '.strelka.snv.VEP91.vcf')
    else:
        sys.stderr.write(date_time() + 'No PASS calls for ' + strelka_snv_pass + ', skipping annotation!\n')
    check_indel = 'cat ' + strelka_indel_pass + ' grep PASS | wc -l'
    check_indel = check_output(check_indel, shell=True)
    if check_indel.decode().rstrip('\n') == '0':
        sys.stderr.write(date_time() + 'Starting vep strelka indel annotation for ' + sample_pair + '\n')
        annot_vcf_vep_pipe(config_file, sample_pair, '.strelka.indel_PASS.vcf', '.strelka.indel.VEP91.vcf')
    else:
        sys.stderr.write(date_time() + 'No PASS calls for ' + strelka_indel_pass + ', skipping annotation!\n')
    # cleanup and reorg
    rm_manta = 'rm -rf ' + manta_dir

    mv_vcf = 'mv -t ' + run_dir_prefix + '/OUTPUT ' + ' '.join((strelka_snv_pass, strelka_indel_pass,
                                                                strelka_indel_vcf, strelka_snv_vcf))
    rm_strelka = 'rm -rf ' + strelka_dir
    sys.stderr.write(date_time() + 'Reorganizing files ' + ';'.join((rm_manta, mv_vcf, rm_strelka)))
    call(rm_manta, shell=True)
    call(mv_vcf, shell=True)
    call(rm_strelka, shell=True)
    set_acls(run_dir_prefix, user, group)
    set_acls(anno_complete_path, user, group)
    sys.stderr.write(date_time() + 'Strelka somatic variant calling for ' + sample_pair + ' complete\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Strelka somatic caller wrapper script.  Samples must have been aligned and bams merged '
                    'ahead of time')
    parser.add_argument('-t', '--tumor', action='store', dest='tumor', help='Tumor sample id')
    parser.add_argument('-n', '--normal', action='store', dest='normal', help='Normal sample id')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (tumor_id, normal_id, config_file) = (inputs.tumor, inputs.normal, inputs.config_file)
    run_strelka(tumor_id, normal_id, config_file)
