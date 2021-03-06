#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.log import log
import subprocess


def filter_wrap(mmu_filter, bwa_tool, RGRP, bwa_ref, end1, end2, samtools_tool, samtools_ref, sample, log_dir, threads):
    loc = log_dir + sample + ".mmu.bwa.pe.log"
    bwa_cmd = "(" + bwa_tool + " mem -O 60 -L 0 -E 10 -t " + threads + " -R \"" + RGRP + "\" -v 2 " + bwa_ref + " "\
              + end1 + " " + end2 + " 2>> " + loc + " | " + samtools_tool + " view -bT " \
              + samtools_ref + " - 2>> " + loc + " | tee " + sample + ".mmu.bam | python " \
              + mmu_filter + " -s " + sample + " -n 0 -t DNA | gzip -4 -c - > " \
              + sample + "_1.filtered.fq.gz;) 2>&1 | gzip -4 -c - > " + sample + "_2.filtered.fq.gz"

    log(loc, date_time() + bwa_cmd + "\n")
    try:
        subprocess.check_output(bwa_cmd, shell=True).decode()
        log(loc, date_time() + 'Filtering completed, replacing fastq file\n')
        rn_fq = 'mv ' + sample + '_1.filtered.fq.gz ' + end1 + '; mv ' + sample + '_2.filtered.fq.gz ' + end2
        subprocess.call(rn_fq, shell=True)
    except:
        sys.stderr.write('Filtering failed\n.')
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mouse read filter alignment module. Fairly PDX-specific process.')
    parser.add_argument('-mmu', '--mmu_filter', action='store', dest='mmu_filter',
                        help='Location of bam filter tool.')
    parser.add_argument('-b', '--bwa', action='store', dest='bwa_tool',
                        help='Location of bwa alignment tool.  Version 0.7.8 preferred.')
    parser.add_argument('-rg', '--RGRP', action='store', dest='RGRP', help='SAM header read group string')
    parser.add_argument('-br', '--bwa_reference', action='store', dest='bwa_ref', help='Location of bwa reference file')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-s', '--samtools', action='store', dest='samtools_tool',
                        help='Location of samtools tool.  Version 1.19 preferred.')
    parser.add_argument('-sr', '--samtools_reference', action='store', dest='samtools_ref',
                        help='Location of samtools reference')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-t', '--threads', action='store', dest='threads',
                        help='Number of threads to use.  8 recommended for standard vm')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (mmu_filter, bwa_tool, RGRP, bwa_ref, end1, end2, samtools_tool, samtools_ref, sample, log_dir, threads) = (
        inputs.mmu_filter, inputs.bwa_tool, inputs.RGRP, inputs.bwa_ref, inputs.end1, inputs.end2, inputs.samtools_tool,
        inputs.samtools_ref, inputs.sample, inputs.log_dir, inputs.threads)
    filter_wrap(mmu_filter, bwa_tool, RGRP, bwa_ref, end1, end2, samtools_tool, samtools_ref, sample, log_dir, threads)
