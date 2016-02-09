#!/usr/bin/env python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import os
import re
from date_time import date_time
from cutadapter import cutadapter
from filter_wrap import filter_wrap
from fastx import fastx
from bwa_mem_pe import bwa_mem_pe
from novosort_sort_pe import novosort_sort_pe
from picard_rmdup import picard_rmdup
from picard_insert_size import picard_insert_size
from flagstats import flagstats
import coverage
from subprocess import call
import subprocess
import json
from log import log
from parse_qc import parse_qc


class Pipeline():
    def __init__(self, end1, end2, seqtype, json_config, ref_mnt):
        self.json_config = json_config
        self.end1 = end1
        self.end2 = end2
        self.seqtype = seqtype
        self.status = 0
        self.ref_mnt = ref_mnt
        self.parse_config()

    def parse_config(self):
        self.config_data = json.loads(open(self.json_config, 'r').read())
        s = re.match('^(\S+)_1_sequence\.txt\.gz$', self.end1)
        if s:
            self.sample = s.group(1)
        else:
            s = re.match('(^\S+)_\D*\d\.f\w*q\.gz$', self.end1)
            self.sample = s.group(1)
        self.loc = 'LOGS/' + self.sample + '.pipe.log'
        HGACID = self.sample.split("_")
        self.bid = HGACID[0]
        self.cutadapter = self.config_data['tools']['cutadapt']
        self.mmu_filter = self.config_data['tools']['mouse_filter']
        self.fastx_tool = self.config_data['tools']['fastx']
        self.bwa_tool = self.config_data['tools']['bwa']
        self.mmu_bwa_ref = self.ref_mnt + '/' + self.config_data['refs']['mmu_bwa']
        self.hsa_bwa_ref = self.ref_mnt + '/' + self.config_data['refs']['hsa_bwa']
        self.samtools_tool = self.config_data['tools']['samtools']
        self.mmu_samtools_ref = self.ref_mnt + '/' + self.config_data['refs']['mmu_samtools']
        self.hsa_samtools_ref = self.ref_mnt + '/' + self.config_data['refs']['hsa_samtools']
        self.java_tool = self.config_data['tools']['java']
        self.picard_tool = self.config_data['tools']['picard']
        self.novosort = self.config_data['tools']['novosort']
        self.picard_tmp = 'picard_tmp'
        self.bedtools2_tool = self.config_data['tools']['bedtools']
        self.bed_ref = self.ref_mnt + '/' + self.config_data['refs'][self.seqtype]
        self.obj = self.config_data['refs']['obj']
        self.cont = self.config_data['refs']['cont']
        self.qc_stats = self.config_data['tools']['qc_stats']
        self.threads = self.config_data['params']['threads']
        self.run_cut_flag = self.config_data['params']['cutflag']
        self.ram = self.config_data['params']['ram']
        self.cflag = 'y'
        if self.seqtype == 'capture':
            self.cflag = 'n'

        self.pipeline()

    def pipeline(self):
        log_dir = 'LOGS/'
        if not os.path.isdir(log_dir):
            mk_log_dir = 'mkdir ' + log_dir
            call(mk_log_dir, shell=True)
            log(self.loc, date_time() + 'Made log directory ' + log_dir + "\n")
        # create BAM and QC directories if they don't exist already
        bam_dir = 'BAM/'
        qc_dir = 'QC/'
        if not os.path.isdir(bam_dir):
            mk_bam_dir = 'mkdir ' + bam_dir
            call(mk_bam_dir, shell=True)
            log(self.loc, date_time() + 'Made bam directory ' + bam_dir + "\n")
        if not os.path.isdir(qc_dir):
            mk_qc_dir = 'mkdir ' + qc_dir
            call(mk_qc_dir, shell=True)
            log(self.loc, date_time() + 'Made qc directory ' + qc_dir + "\n")
        log(self.loc,
            date_time() + "Starting alignment qc for paired end sample files " + self.end1 + " and " + self.end2 + "\n")
        # inputs

        SAMPLES = {}
        SAMPLES[self.sample] = {}
        SAMPLES[self.sample]['f1'] = self.end1
        SAMPLES[self.sample]['f2'] = self.end2
        RGRP = "@RG\\tID:" + self.sample + "\\tLB:" + self.bid + "\\tSM:" + self.bid + "\\tPL:illumina"
        if self.run_cut_flag == 'Y':
            check = cutadapter(self.sample, self.end1, self.end2, self.json_config)
            if check != 0:
                log(self.loc, date_time() + 'cutadapt failure for ' + self.sample + '\n')
                exit(1)
        log(self.loc, date_time() + 'Aligning and filtering reads for mouse contamination')
        check = filter_wrap(self.mmu_filter, self.bwa_tool, RGRP, self.mmu_bwa_ref, self.end1, self.end2,
                            self.samtools_tool, self.mmu_samtools_ref, self.sample, log_dir, self.threads)
        if check != 0:
            log(self.loc, date_time() + 'Read filter failure for ' + self.sample + '\n')
            exit(1)
        log(self.loc, date_time() + 'Starting BWA align\n')
        wait_flag = 1
        # check certain key processes
        # skip aligning if bam already exists
        if not os.path.isfile(self.sample + '.bam'):
            check = bwa_mem_pe(self.bwa_tool, RGRP, self.hsa_bwa_ref, self.end1, self.end2, self.samtools_tool,
                               self.hsa_samtools_ref, self.sample, log_dir, self.threads)
            # rest won't run until completed
            if check != 0:
                log(self.loc, date_time() + 'BWA failure for ' + self.sample + '\n')
                exit(1)
            log(self.loc, date_time() + 'Getting fastq quality score stats\n')
            fastx(self.fastx_tool, self.sample, self.end1, self.end2)  # will run independently of rest of output
        else:
            log(self.loc, date_time() + 'bam file already exists, skipping alignment as well as fastx!\n')
        # skip sort if sorted file exists already
        log(self.loc, date_time() + 'Sorting BAM file\n')
        if not os.path.isfile(self.sample + '.srt.bam'):
            check = novosort_sort_pe(self.novosort, self.sample, log_dir, self.threads,
                                     self.ram)  # rest won't run until completed
            if check != 0:
                log(self.loc, date_time() + 'novosort sort failure for ' + self.sample + '\n')
                exit(1)
        else:
            log(self.loc, date_time() + 'Sorted bam file already exists, skipping\n')
        # skip next steps in insert size already calculated
        if not os.path.isfile(self.sample + '.insert_metrics.hist'):
            log(self.loc, date_time() + 'Removing PCR duplicates\n')
            picard_rmdup(self.java_tool, self.picard_tool, self.picard_tmp, self.sample, log_dir,
                         self.ram)  # rest won't run until emopleted
            log(self.loc, date_time() + 'Gathering SAM flag stats\n')
            flagstats(self.samtools_tool,
                      self.sample)  # flag determines whether to run independently or hold up the rest of the pipe until completion
            log(self.loc, date_time() + 'Calculating insert sizes\n')
            picard_insert_size(self.java_tool, self.picard_tool, self.sample, log_dir,
                               self.ram)  # get insert size metrics.
        else:
            log(self.loc,
                date_time() + 'Insert size file detected, skipping remove duplicates, flagstats, and remove duplicates steps')
        # figure out which coverage method to call using seqtype
        log(self.loc, date_time() + 'Calculating coverage for ' + self.seqtype + '\n')
        method = getattr(coverage, (self.seqtype + '_coverage'))

        method(self.bedtools2_tool, self.sample, self.bed_ref,
               wait_flag)  # run last since this step slowest of the last
        log(self.loc, date_time() + 'Checking outputs and uploading results\n')
        # check to see if last expected file was generated search for seqtype + .hist suffix

        flist = os.listdir('./')
        if self.seqtype == 'capture':
            suffix = self.seqtype + '_t2.hist'
        else:
            suffix = self.seqtype + '.hist'

        for fn in flist:
            if fn == (self.sample + '.' + suffix):
                self.status = 1
                break
        if self.status == 1:
            p_tmp_rm = "rm -rf picard_tmp"
            call(p_tmp_rm, shell=True)
            # move files into appropriate place and run qc_stats
            log(self.loc, date_time() + 'Calculating qc stats and prepping files for upload\n')
            mv_bam = 'mv *.bam *.bai BAM/'
            subprocess.call(mv_bam, shell=True)
            rm_sf = 'rm ' + self.end1 + ' ' + self.end2
            subprocess.call(rm_sf, shell=True)
            parse_qc(self.json_config, self.sample, self.cflag)
            mv_rest = 'find . -maxdepth 1 -type f -exec mv {} QC \;'
            subprocess.call(mv_rest, shell=True)
            mv_config = ' cp ' + self.json_config + ' QC/'
            subprocess.call(mv_config, shell=True)
            from upload_to_swift import upload_to_swift
            obj = self.obj + "/" + self.bid

            check = upload_to_swift(self.cont, obj)
            if check == 0:
                log(self.loc, date_time() + 'Couchdb successfully updated\n')
                self.status = 0
                log(self.loc,
                    date_time() + "Pipeline complete, files successfully uploaded.  Files may be safely removed\n")
            else:
                log(self.loc, date_time() + "All but file upload succeeded\n")
                self.status = 1

        else:
            (self.loc,
             date_time() + "File with suffix " + suffix + " is missing!  If intentional, ignore this message.  Otherwise, check logs for potential failures\n")
            self.status = 1


def main():
    import argparse
    parser = argparse.ArgumentParser(description='DNA alignment paired-end QC pipeline')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second fastq file')
    parser.add_argument('-t', '--seqtype', action='store', dest='seqtype',
                        help='Type of sequencing peformed.  Likely choices are genome, exome, and capture')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    end1 = inputs.end1
    end2 = inputs.end2
    seqtype = inputs.seqtype
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt
    Pipeline(end1, end2, seqtype, config_file, ref_mnt)


if __name__ == "__main__":
    main()
