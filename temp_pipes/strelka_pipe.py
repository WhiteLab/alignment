#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
from utility.date_time import date_time
import subprocess
from annotation.annot_vcf_vep import parse_config
from annotation.annot_vcf_vep import run_vep
from batch_coverage import get_bam_name
from annotation.vep_strelka_report import gen_report


def run_strelka(strelka_tools, norm_bam, tum_bam, pair, threads, fasta):
    cwd = os.getcwd()
    strelka_config = strelka_tools + ' --normal=' + norm_bam + ' --tumor=' + tum_bam + ' --ref=' + fasta\
                     + ' --config=' + cwd + '/config.ini --output-dir=' + cwd + '/' + pair
    sys.stderr.write(date_time() + strelka_config + '\n')
    check = subprocess.call(strelka_config, shell=True)
    run_cmd = 'cd ' + pair + ' && make -j ' + threads
    sys.stderr.write(date_time() + run_cmd + '\n')
    check += subprocess.call(run_cmd)
    return check


def annot_strelka_pipe(strelka_tools, pairs, config, ref_mnt):
    (vep_tool, vep_cache, fasta, report, dbsnp, vcache, threads, intvl, dustmask_flag, wg_flag, tx_index) \
        = parse_config(config)
    fasta = ref_mnt + '/' + fasta
    vep_cache = ref_mnt + '/' + vep_cache
    intvl = ref_mnt + '/' + intvl
    tx_index = ref_mnt + '/' + tx_index
    src_cmd = '. /home/ubuntu/.novarc;'
    obj = 'ALIGN'
    cwd = os.getcwd()
    # scale back on the forking a bit for VEP
    vep_threads = threads
    if int(threads) > 2:
        # threads = str(int(threads)/2 - 1)
        vep_threads = str(int(threads) - 1)

    for pair in open(pairs):
        pair = pair.rstrip('\n')
        loc = pair + '.vep.log'
        sys.stderr.write(date_time() + 'Processing ' + pair + '\n')
        bnids = pair.split('_')
        # get bams
        (dl_cmd, tum_bam, tum_bai) = get_bam_name(bnids[0], src_cmd, 'PDX', obj)
        if len(tum_bam) < 1:
            sys.stderr.write('Did not find bam for ' + bnids[0] + ' in PDX container trying PANCAN\n')
            (dl_cmd, tum_bam, tum_bai) = get_bam_name(bnids[0], src_cmd, 'PANCAN', obj)
            if len(tum_bam) < 1:
                sys.stderr.write(date_time() + 'Could not find tumor bam ' + bnids[0] + '! Aborting\n')
                exit(1)
        sys.stderr.write(date_time() + dl_cmd + '\n')
        subprocess.call(dl_cmd, shell=True)

        (dl_cmd, norm_bam, norm_bai) = get_bam_name(bnids[1], src_cmd, 'PANCAN', obj)
        if len(norm_bam) < 1:
            sys.stderr.write('Did not find bam for ' + bnids[1] + ' in PANCAN container trying PDX\n')
            (dl_cmd, norm_bam, norm_bai) = get_bam_name(bnids[1], src_cmd, 'PDX', obj)
            if len(norm_bam) < 1:
                sys.stderr.write(date_time() + 'Could not find normal bam ' + bnids[1] + '! Aborting\n')
        sys.stderr.write(date_time() + dl_cmd + '\n')
        subprocess.call(dl_cmd, shell=True)

        check = run_strelka(strelka_tools, norm_bam, tum_bam, pair, threads, fasta)
        if check != 0:
            sys.stderr.write(date_time() + 'Running strelka for ' + pair + ' failed!\n')
            exit(1)
        cleanup = ' '.join(('rm', norm_bam, norm_bai, tum_bam, tum_bai))
        sys.stderr.write(date_time() + cleanup)
        subprocess.call(cleanup, shell=True)
        buffer_size = '2000'

        in_vcf = cwd + '/' + pair + '/results/passed.somatic.snvs.vcf'
        out_vcf = cwd + '/' + pair + '/results/' + pair + '.somatic.snv.strelka.vep.vcf'
        run_cmd = run_vep(wg_flag, vep_tool, in_vcf, out_vcf, buffer_size, vep_threads, fasta, vep_cache, vcache, loc)
        check = subprocess.call(run_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'vep failed for ' + in_vcf + ' for pair ' + pair + '\n')
            exit(1)
        gen_report(out_vcf, intvl, tx_index)

        in_vcf = cwd + '/' + pair + '/results/passed.somatic.indels.vcf'
        out_vcf = cwd + '/' + pair + '/results/' + pair + '.somatic.indel.strelka.vep.vcf'
        run_cmd = run_vep(wg_flag, vep_tool, in_vcf, out_vcf, buffer_size, vep_threads, fasta, vep_cache, vcache, loc)
        check = subprocess.call(run_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'vep failed for ' + in_vcf + ' for pair ' + pair + '\n')
            exit(1)
        gen_report(out_vcf, intvl, tx_index)

    sys.stderr.write(date_time() + 'FIN.\n')






if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Quick n dirty annotation with strelka. Container locations hardcoded')
    parser.add_argument('-s', '--strelka', action='store', dest='strelka', help='Location of strelka tools.')
    parser.add_argument('-j', '--json', action='store', dest='config',
                        help='config file with reference locations ')
    parser.add_argument('-p', '--pairs', action='store', dest='pairs', help='bnid pairs in \\n separated file, '
                                                                            'bnids seprated by _')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    strelka_tools = inputs.strelka
    pairs = inputs.pairs
    config = inputs.json
    ref_mnt = inputs.ref_mnt

    annot_strelka_pipe(strelka_tools, pairs, config, ref_mnt)
