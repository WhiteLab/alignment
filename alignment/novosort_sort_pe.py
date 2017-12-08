#!/usr/bin/env python

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import subprocess
from utility.log import log


def novosort_sort_pe(novosort, sample, log_dir, threads, ram, rmdup):
    if rmdup == 'Y':
        logfile = sample + ".novosort.rmdup.sort.pe.log"
        novosort_sort_cmd = 'mkdir novosort_tmp;' + novosort + " -c " + threads + " -m " + ram \
                            + "G --tmpdir novosort_tmp --rd --kt -o " + sample + ".rmdup.srt.bam --index  "\
                               + sample + ".bam > " + log_dir + logfile + " 2>&1"
        log(log_dir + logfile, date_time() + novosort_sort_cmd + "\n")
    else:
        logfile = sample + ".novosort.sort.pe.log"
        novosort_sort_cmd = 'mkdir novosort_tmp;' + novosort + " --threads " + threads + " --ram " \
                               + ram + "G --tmpdir novosort_tmp -o " + sample + ".srt.bam --index  " \
                               + sample + ".bam > " + log_dir + logfile + " 2>&1"
        log(log_dir + logfile, date_time() + novosort_sort_cmd + "\n")
    f = 0
    try:
        f = subprocess.call(novosort_sort_cmd, shell=True)
        rm_tmp = 'rm -rf novosort_tmp'
        subprocess.call(rm_tmp, shell=True)
    except:
        log(log_dir + logfile, 'novosort sort failed for sample ' + sample + '\n')
        exit(1)
    return f


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to sort BAM module.')
    parser.add_argument('-n', '--novosort', action='store', dest='novosort', help='novosort binary location')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-t', '--threads', action='store', dest='threads',
                        help='Number of threads to use. Assuming standard vm, 8 recommonded ')
    parser.add_argument('-r', '--ram', action='store', dest='ram',
                        help='RAM to use in GB.  24 recommended for standard vm')
    parser.add_argument('-f', '--flag', action='store', dest='rmdup',
                        help='To also remove duplicates while sorting, give Y')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (novosort, sample, log_dir, threads, ram, rmdup) = (
        inputs.novosort, inputs.sample, inputs.log_dir, inputs.threads, inputs.ram, inputs.rmdup)
    novosort_sort_pe(novosort, sample, log_dir, threads, ram, rmdup)
