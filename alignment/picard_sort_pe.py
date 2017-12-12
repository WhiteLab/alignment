#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
import subprocess
from utility.log import log


def picard_sort_pe(java_tool, picard_tool, picard_tmp, sample, log_dir):
    picard_sort_pe_cmd = java_tool + " -Xmx8g -jar " + picard_tool + " SortSam CREATE_INDEX=true TMP_DIR=" \
                         + picard_tmp + " INPUT=" + sample + ".bam OUTPUT=" + sample + ".srt.bam SORT_ORDER=" \
                        "coordinate VALIDATION_STRINGENCY=LENIENT > " + log_dir + sample + ".picard.sort.pe.log 2>&1"
    log(log_dir + sample + ".picard.sort.pe.log", date_time() + picard_sort_pe_cmd + "\n")
    try:
        subprocess.check_output(picard_sort_pe_cmd, shell=True)
    except:
        log(log_dir + sample + ".picard.sort.pe.log",
            'Picard sort failed for sample ' + sample + '.  Check for borg!\n')
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Picard tools to sort BAM module.')
    parser.add_argument('-j', '--java', action='store', dest='java_tool',
                        help='Java location directory, version jdk1.7.0_45 preferred')
    parser.add_argument('-p', '--picard', action='store', dest='picard_tool', help='Picard jar file location')
    parser.add_argument('-pt', '--picard_temp', action='store', dest='picard_tmp',
                        help='Picard temp folder location to create')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (java_tool, picard_tool, picard_tmp, sample, log_dir) = (
        inputs.java_tool, inputs.picard_tool, inputs.picard_tmp, inputs.sample, inputs.log_dir)
    picard_sort_pe(java_tool, picard_tool, picard_tmp, sample, log_dir)
