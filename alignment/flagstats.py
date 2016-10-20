#!/usr/bin/python
import sys
import os
from utility.date_time import date_time
from subprocess import Popen

def flagstats(samtools_tool, sample):

    # test for sorted bam, otherwise use unsorted bam
    raw_bam = sample + ".srt.bam"
    res_file = sample + ".srt.bam.flagstats"
    if os.path.isfile(raw_bam) == False:
        raw_bam = sample + ".bam"
        res_file = sample + ".bam.flagstats"

    flagstats_cmd = samtools_tool + " flagstat " + raw_bam + " > " + res_file
    sys.stderr.write(date_time() + flagstats_cmd + "\n")
    Popen(flagstats_cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)


    flagstats_cmd = samtools_tool + " flagstat " + sample + ".rmdup.srt.bam > " + sample + ".rmdup.srt.bam.flagstats"
    sys.stderr.write(date_time() + flagstats_cmd + "\n")
    Popen(flagstats_cmd, shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Flag stats from samtools module.  Assumes bwa alignent and picard tools or novosort have been run '
                    'to sort bam and remove duplicates.',
        add_help=True)
    parser.add_argument('-s', '--samtools', action='store', dest='samtools_tool',
                        help='Location of samtools tool.  Version 1.19 preferred.')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (samtools_tool, sample) = (inputs.samtools_tool, inputs.sample)
    flagstats(samtools_tool, sample)
