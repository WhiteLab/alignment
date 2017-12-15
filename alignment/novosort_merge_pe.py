#!/usr/bin/env python3

import json
import re
import sys
import os
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.log import log
import subprocess
from utility.job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['novosort'], config_data['tools']['java'], config_data['tools']['picard'], \
           config_data['refs']['project'], config_data['refs']['align'], config_data['params']['threads'], \
           config_data['params']['ram'], config_data['params']['novaflag'], \
           config_data['tools']['novo_merge_rmdup_slurm']


def list_unsorted_bam(bam_list, cont, threads):
    p = []
    for i in range(0, len(bam_list), 1):
        bam_list[i] = bam_list[i].replace('.rmdup.srt.bam', '.bam')
        dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' ' + bam_list[i]
        p.append(dl_cmd)
    f = job_manager(p, threads)
    if f == 0:
        sys.stderr.write(date_time() + 'BAM download complete\n')
        return bam_list

    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)


def list_bam(project, align, sample,rmdup):
    bam_dir = '/cephfs/PROJECTS/' + project + '/' + align + '/' + sample + '/BAM/'
    find_bam_cmd = 'find ' + bam_dir + '*.bam'
    sys.stderr.write(date_time() + find_bam_cmd + '\nGetting BAM list\n')
    bam_find = subprocess.check_output(find_bam_cmd, shell=True).decode()

    bam_list = bam_find.split('\n')
    find_bai_cmd = 'find ' + bam_dir + '*.bai'
    sys.stderr.write(date_time() + find_bai_cmd + '\nGetting bai list\n')
    bai_find = subprocess.check_output(find_bai_cmd, shell=True).decode()
    bai_list = bai_find.split('\n')
    ct = len(bam_list)
    if ct >= 1:
        sys.stderr.write(date_time() + 'BAM files\n')
        if rmdup == 'Y':
            return bam_list, ct
        else:
            return bam_list, bai_list, ct
    else:
        sys.stderr.write(date_time() + 'No bams found for ' + sample + '\n')
        exit(1)


def novosort_merge_pe(config_file, sample_list):
    fh = open(sample_list, 'r')
    (novosort, java_tool, picard_tool, project, align, threads, ram, rmdup, novo_merge_rmdup_slurm) \
        = parse_config(config_file)

    for sample in fh:
        sample = sample.rstrip('\n')
        loc = 'LOGS/' + sample + '.novosort_merge.log'
        if rmdup == 'Y':
            (bam_list, n) = list_bam(project, align, sample, rmdup)
        else:
            (bam_list, bai_list, n) = list_bam(project, align, sample, rmdup)
        bam_string = " ".join(bam_list)
        cur_dir = '/cephfs/PROJECTS/' + project + '/' + align + '/' + sample + '/BAM/'
        os.chdir(cur_dir)
        tmp_dir = 'mkdir TMP'
        job_log = sample + '.novosort_merge.log'
        subprocess.call(tmp_dir, shell=True)
        if n > 1:
            if rmdup == 'Y':
                # $novosort -c $threads -m $ram --rd -o $sample .merged.final.bam --index --tmpdir ./TMP $bam_string
                # 2>> $loc
                batch = 'sbatch -c ' + threads + ' --mem ' + ram + ' -o ' + job_log \
                + ' --export=novosort="' + novosort + '",threads="' + threads + '",ram="' + ram + 'G",$bam_string="' \
                        + bam_string + '",loc="' + loc + '"' + ' ' + novo_merge_rmdup_slurm
                log(loc, date_time() + 'Submitting merge bam job for sample ' + batch + "\n")

                subprocess.call(batch, shell=True)

            else:
                novosort_merge_cmd = novosort + " --threads " + threads + " --ram " + ram + "G --assumesorted -o "\
                                        + sample + '.merged.bam --index --tmpdir ./TMP ' + bam_string
                log(loc, date_time() + novosort_merge_cmd + "\n")
                try:
                    subprocess.check_output(novosort_merge_cmd, shell=True).decode()
                    # delete old bams to free up space
                    rm_bam = 'rm ' + bam_string
                    log(loc, date_time() + 'Removing bams that were already merged\n')
                    subprocess.call(rm_bam, shell=True)
                    # rm dups
                    picard_tmp = 'picard_tmp'
                    # setting max records in ram to half of ram
                    recs = (int(ram) / 2) * (1000000000 / 200)
                    picard_rmdup_cmd = java_tool + " -Xmx" + ram + "g -jar " + picard_tool + " MarkDuplicates " \
                                       "CREATE_INDEX=true TMP_DIR=" + picard_tmp + " REMOVE_DUPLICATES=true" \
                                       " ASSUME_SORTED=true MAX_RECORDS_IN_RAM=" + str(recs) + \
                                       " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 INPUT=" + sample + ".merged.bam OUTPUT="\
                                       + sample + ".merged.final.bam METRICS_FILE=" + sample +\
                                       ".rmdup.srt.metrics VALIDATION_STRINGENCY=LENIENT"
                    sys.stderr.write(date_time() + 'Removing dups\n' + picard_rmdup_cmd + '\n')
                    subprocess.call(picard_rmdup_cmd, shell=True)
                    # delete merged bam after removing dups
                    rm_merged_bam = 'rm ' + sample + '.merged.bam ' + sample + '.merged.bai'
                    subprocess.call(rm_merged_bam, shell=True)
                except:
                    sys.stderr.write(date_time() + 'novosort and picard merge failed for sample ' + sample + '\n')
                    exit(1)
        else:
            try:
                # first need to get bai file before renaming
                if rmdup == 'N':
                    dl_cmd = '. /home/ubuntu/.novarc;swift download ' + project + ' --skip-identical ' + bai_list[0]
                    check = subprocess.call(dl_cmd, shell=True)
                    if check != 0:
                        log(loc, 'Could not find bai file for ' + bam_list[0])
                        exit(1)
                    mv_bam = 'mv ' + bam_list[0] + ' ' + sample + '.merged.final.bam; mv ' + bai_list[0]\
                            + ' ' + sample + '.merged.final.bam.bai'
                    subprocess.call(mv_bam, shell=True)
                else:
                    # get name of rmdup bai file
                    root = bam_list[0].replace('.bam', '')
                    get_name = '. /home/ubuntu/.novarc; swift list ' + project + ' --prefix ' + root + ' | grep bai$'
                    bai = subprocess.check_output(get_name, shell=True).decode()
                    bai = bai.rstrip('\n')
                    dl_cmd = '. /home/ubuntu/.novarc;swift download ' + project + ' --skip-identical ' + bai
                    check = subprocess.call(dl_cmd, shell=True)
                    if check != 0:
                        log(loc, 'Could not find bai file for ' + bai)
                        exit(1)
                    mv_bam = 'mv ' + bam_list[0] + ' ' + sample + '.merged.final.bam; mv ' + bai + ' ' + sample \
                             + '.merged.final.bam.bai'
                    log(loc, date_time() + mv_bam + '\n')
                    subprocess.call(mv_bam, shell=True)
            except:
                log(loc, 'Rename for single file failed.  Command was ' + mv_bam + '\n')
                exit(1)
            sys.stderr.write(date_time() + mv_bam + ' Only one associated bam file, renaming\n')
            subprocess.call(mv_bam, shell=True)
    sys.stderr.write(date_time() + 'Merge process complete\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample/project prefix list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file) = (inputs.sample_list, inputs.config_file)
    novosort_merge_pe(config_file, sample_list)
