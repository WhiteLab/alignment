#!/usr/bin/python
import json
import re
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from log import log
import subprocess
from job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['novosort'], config_data['tools']['java'], config_data['tools']['picard'], \
           config_data['refs']['cont'], config_data['refs']['obj'], config_data['params']['threads'], \
           config_data['params']['ram'], config_data['params']['novaflag']


def list_bam(cont, obj, sample, threads, rmdup):
    ct = 0
    # added trailing slash since we're dealing with bids - otherwise will end up pulling more samples than intended
    list_cmd = '. /home/ubuntu/.novarc;swift list ' + cont + ' --prefix ' + obj + '/' + sample + '/BAM/'
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    p = []

    bam_list = []
    bai_list = []
    for fn in re.findall('(.*)\n', flist):
        # depending on software used to index, .bai extension may follow bam
        if rmdup == 'N':
            test = re.match('^\S+_\w*\d+\.rmdup.srt.ba[m|i]$', fn) or re.match('^\S+_\w*\d+\.rmdup.srt.bam.bai$', fn)
        else:
            test = re.match('^\S+_\w*\d+\.bam$', fn)
        if test:
            dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + fn
            p.append(dl_cmd)
            if fn[-3:] == 'bam':
                bam_list.append(fn)
                ct += 1
            else:
                bai_list.append(fn)
    f = 0
    if len(p) > 1:
        f = job_manager(p, threads)
    # need to replace with rmdup version if nothing there to merge
    elif rmdup == 'Y':
        cur = bam_list[0]
        new = cur.replace('.bam', '.rmdup.srt.bam')
        new_dl = '. /home/ubuntu/.novarc; swift download ' + cont + ' ' + new
        subprocess.call(new_dl, shell=True)
        bam_list[0] = new
    if f == 0:
        sys.stderr.write(date_time() + 'BAM download complete\n')
        if rmdup == 'Y':
            return bam_list, ct
        else:
            return bam_list, bai_list, ct
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)


def novosort_merge_pe(config_file, sample_list, wait):
    fh = open(sample_list, 'r')
    (novosort, java_tool, picard_tool, cont, obj, threads, ram, rmdup) = parse_config(config_file)
    tmp_dir = 'mkdir TMP'
    subprocess.call(tmp_dir, shell=True)
    for sample in fh:
        sample = sample.rstrip('\n')
        loc = 'LOGS/' + sample + '.novosort_merge.log'
        if rmdup == 'Y':
            (bam_list, n) = list_bam(cont, obj, sample, threads, rmdup)
        else:
            (bam_list, bai_list, n) = list_bam(cont, obj, sample, threads, rmdup)
        bam_string = " ".join(bam_list)
        if n > 1:
            if rmdup == 'Y':
                novosort_merge_cmd = novosort + " -c " + threads + " -m " + ram + "G -f --rd --output " + sample \
                                     + '.merged.final.bam --index --tmpdir ./TMP ' + bam_string + ' 2>> ' + loc
                log(loc, date_time() + novosort_merge_cmd + "\n")
                try:
                    subprocess.check_output(novosort_merge_cmd, shell=True)
                    # delete old bams to free up space
                    rm_bam = 'rm ' + bam_string
                    log(loc, date_time() + 'Removing bams that were already merged\n')
                    subprocess.call(rm_bam, shell=True)
                except:
                    log(loc, date_time() + 'novosort sort and merge failed for sample ' + sample + '\n')
                    exit(1)

            else:
                novosort_merge_cmd = novosort + " --threads " + threads + " --ram " + ram + "G --assumesorted --output "\
                                        + sample + '.merged.bam --index --tmpdir ./TMP ' + bam_string
                log(loc, date_time() + novosort_merge_cmd + "\n")
                try:
                    subprocess.check_output(novosort_merge_cmd, shell=True)
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
                    dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + bai_list[0]
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
                    get_name = '. /home/ubuntu/.novarc; swift list ' + cont + ' --prefix ' + root + ' | grep bai$'
                    bai = subprocess.check_output(get_name, shell=True)
                    bai = bai.rstrip('\n')
                    dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + bai
                    check = subprocess.call(dl_cmd, shell=True)
                    if check != 0:
                        log(loc, 'Could not find bai file for ' + bai)
                        exit(1)
                    mv_bam = 'mv ' + bam_list[0] + ' ' + sample + '.merged.final.bam; mv ' + bai + ' ' + sample \
                             + '.merged.final.bam.bai'
                    log(loc, date_time() + mv_bam)
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
    parser.add_argument('-w', '--wait', action='store', dest='wait',
                        help='Wait time to download bam files.  900 (seconds) recommended')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file, wait) = (inputs.sample_list, inputs.config_file, inputs.wait)
    novosort_merge_pe(config_file, sample_list, wait)
