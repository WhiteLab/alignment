#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import json
from utility.date_time import date_time
from subprocess import call
from utility.log import log
from analysis.dustmask_filter import filter_indel


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['scalpel'], config_data['tools']['bedtools'], config_data['refs']['capture'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads'], config_data['params']['dustmask_flag'],\
           config_data['refs']['dustmask'], config_data['params']['wg_flag'], config_data['refs']['project_dir'], \
           config_data['refs']['project'], config_data['refs']['align_dir']


def wg_mode(scalpel, tumor_bam, normal_bam, fasta, cpus, pair, config_file):
    config_data = json.loads(open(config_file, 'r').read())
    exome = config_data['refs']['exome']
    loc = 'LOGS/' + pair + '_' + pair + '.genome_as_exome.scalpel.log'
    cmd = scalpel + ' --somatic --logs --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' \
    + normal_bam + ' --window 600 --two-pass --bed ' + exome + ' --ref ' + fasta + ' 2> ' + loc
    log(loc, date_time() + cmd + '\n')
    check = call(cmd, shell=True)
    if check != 0:
        return 1, pair
    return 0, pair


def scalpel_indel(tumor_id, normal_id, log_dir, config_file):
    (scalpel, bedtools, bed, fasta, cpus, dustmask_flag, dustmask_bed, wg, project_dir, project, align) \
        = parse_config(config_file)
    
    sample_pair = tumor_id + '_' + normal_id
    loc = log_dir + sample_pair + '.scalpel.log'
    bam_dir = project_dir + project + '/' + align
    tumor_bam = bam_dir + '/' + tumor_id + '/BAM/' + tumor_id + '.merged.final.bam'
    normal_bam = bam_dir + '/' + normal_id + '/BAM/' + normal_id + '.merged.final.bam'
    if wg == 'n':
        scalpel_cmd = scalpel + ' --somatic --logs --numprocs ' + cpus + ' --tumor ' + tumor_bam + ' --normal ' \
                      + normal_bam + ' --bed ' + bed + ' --ref ' + fasta + ' 2>> ' + loc
        sys.stderr.write(date_time() + 'Starting indel calls for ' + sample_pair + '\n')
        log(loc, date_time() + 'Starting indel calls for ' + sample_pair + ' in capture mode with command:\n'
            + scalpel_cmd + '\n')
        check = call(scalpel_cmd, shell=True)
        if check != 0:
            sys.stderr.write(date_time() + 'Indel calling failed for pair ' + sample_pair + ' with command:\n' +
                             scalpel_cmd + '\n')
            exit(1)
    else:
        check = wg_mode(scalpel, tumor_bam, normal_bam, fasta, cpus, sample_pair, config_file)
        if check[0] != 0:
            sys.stderr.write('Scalpel failed for ' + normal_id + ' at ' + tumor_id + '\n')
            exit(1)
    log(loc, date_time() + 'Indel calling complete for pair ' + sample_pair + ' moving output files\n')
    mv_cmd = 'mv outdir/main/* .; rmdir outdir/main;'
    log(loc, date_time() + mv_cmd + '\n')
    call(mv_cmd, shell=True)
    sys.stderr.write(date_time() + 'Completed indel calls for ' + sample_pair + '\n')
    if dustmask_flag == 'Y':
        log(loc, date_time() + 'Filter dustmask flag given\n')
        check = filter_indel(bedtools, dustmask_bed, sample_pair, loc)
        if check != 0:
            sys.stderr.write(date_time() + 'Dustmask failed for ' + sample_pair + '\n')
            exit(1)
        else:
            log(loc, date_time() + 'Dustmask complete for ' + sample_pair + '\n')
    sys.stderr.write(date_time() + 'Indel call completed\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Scalpel indel caller wrapper script.  Samples must have been aligned and bams merged '
                    'ahead of time')
    parser.add_argument('-t', '--tumor', action='store', dest='tumor', help='Tumor sample id')
    parser.add_argument('-n', '--normal', action='store', dest='normal', help='Normal sample id')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (tumor_id, normal_id, log_dir, config_file) = (inputs.tumor, inputs.normal, inputs.log_dir, inputs.config_file)
    scalpel_indel(tumor_id, normal_id, log_dir, config_file)
