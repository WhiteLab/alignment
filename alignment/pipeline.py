#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
import re
from utility.date_time import date_time
from cutadapter import cutadapter
from fastx import fastx
from bwa_mem_pe import bwa_mem_pe
from novosort_sort_pe import novosort_sort_pe
from filter_wrap import filter_wrap
from picard_rmdup import picard_rmdup
from picard_insert_size import picard_insert_size
from flagstats import flagstats
import coverage
from subprocess import call
import subprocess
import json
from utility.log import log
from parse_qc import parse_qc


class Pipeline:
    def __init__(self, end1, end2, seqtype, json_config, ref_mnt):
        self.config_data = json.loads(open(json_config, 'r').read())
        self.end1 = end1
        self.end2 = end2
        s = re.match('^(\S+)_1_sequence\.txt\.gz$', self.end1)
        if s:
            self.sample = s.group(1)
        else:
            s = re.match('(^\S+)_\D*\d\.f\w*q\.gz$', self.end1)
            self.sample = s.group(1)
        hgac_ID = self.sample.split("_")
        self.ref_mnt = ref_mnt
        self.seqtype = seqtype
        self.cflag = 'y'
        if self.seqtype == 'capture':
            self.cflag = 'n'

        # flag for whether to run adapter/quality trimming
        self.run_cut_flag = self.config_data['params']['cutflag']
        self.pdxflag = self.config_data['params']['pdxflag']
        if self.pdxflag == 'Y':
            self.mmu_filter = self.config_data['tools']['mouse_filter']
            self.mmu_bwa_ref = self.ref_mnt + '/' + self.config_data['refs']['mmu_bwa']
            self.hsa_bwa_ref = self.ref_mnt + '/' + self.config_data['refs']['hsa_bwa']
            self.mmu_samtools_ref = self.ref_mnt + '/' + self.config_data['refs']['mmu_samtools']
            self.hsa_samtools_ref = self.ref_mnt + '/' + self.config_data['refs']['hsa_samtools']
        else:
            self.samtools_ref = self.ref_mnt + '/' + self.config_data['refs']['samtools']
            self.bwa_ref = self.ref_mnt + '/' + self.config_data['refs']['bwa']

        # flag for whether to use novosort rmdup capabilities
        self.use_nova_flag = self.config_data['params']['novaflag']
        self.ram = self.config_data['params']['ram']
        self.threads = self.config_data['params']['threads']
        self.qc_stats = self.config_data['tools']['qc_stats']
        self.cont = self.config_data['refs']['cont']
        self.obj = self.config_data['refs']['obj']
        self.bed_ref = self.ref_mnt + '/' + self.config_data['refs'][self.seqtype]
        self.bedtools2_tool = self.config_data['tools']['bedtools']
        self.picard_tmp = 'picard_tmp'
        self.novosort = self.config_data['tools']['novosort']
        self.picard_tool = self.config_data['tools']['picard']
        self.java_tool = self.config_data['tools']['java']
        self.samtools_tool = self.config_data['tools']['samtools']
        self.bwa_tool = self.config_data['tools']['bwa']
        self.fastx_tool = self.config_data['tools']['fastx']
        self.cutadapter = self.config_data['tools']['cutadapt']
        self.bid = hgac_ID[0]
        self.sample = s.group(1)
        self.loc = 'LOGS/' + self.sample + '.pipe.log'
        self.json_config = json_config

        self.status = 0
        self.pipeline()

    def check_outputs(self):
        # check for outputs and ensure file size greater than 0, if not, wait
        if self.seqtype == 'capture':
            suffix = self.seqtype + '_t2.hist'
        else:
            suffix = self.seqtype + '.hist'
        hist = self.sample + '.' + suffix
        ins = self.sample + '.insert_metrics.hist'
        f1 = self.sample + '_1.qs'
        f2 = self.sample + '_2.qs'
        fstat1 = self.sample + '.srt.bam.flagstats'
        if self.use_nova_flag == 'Y':
            fstat1 = self.sample + '.bam.flagstats'
        fstat2 = self.sample + '.rmdup.srt.bam.flagstats'
        status = {hist: 0, ins: 0, f1: 0, f2: 0, fstat1: 0, fstat2: 0}
        # currently giving 20 minutes to complete before giving up
        wait = 1200
        intvl = 30
        # if whole genome, increase wait time substantially
        if self.seqtype == 'genome':
            wait = 14400
            intvl = 300
        cur = 0
        while cur < wait:
            comp_flag = 1
            for fn in status:
                if os.path.isfile(fn) and (os.path.getsize(fn) > 0):
                    status[fn] = 1
                else:
                    sys.stderr.write(date_time() + 'Still waiting on file ' + fn + ' to be created.\n')
                    comp_flag = 0
            if comp_flag == 1:
                self.status = 1
                return 0
            else:
                cur += intvl
                sleep_cmd = 'sleep ' + str(intvl) + 's'
                call(sleep_cmd, shell=True)
        failed = []
        for fn in status:
            if status[fn] == 0:
                failed.append(fn)
        sys.stderr.write('Outputs ' + ', '.join(failed) + ' failed to complete!\n')
        exit(1)

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

        # initialize fail return values
        check = 1
        if self.run_cut_flag == 'Y':
            check = cutadapter(self.sample, self.end1, self.end2, self.json_config)
            if check != 0:
                log(self.loc, date_time() + 'cutadapt failure for ' + self.sample + '\n')
                exit(1)
        if self.pdxflag == 'Y':
            log(self.loc, date_time() + 'Aligning and filtering reads for mouse contamination')
            check = filter_wrap(self.mmu_filter, self.bwa_tool, RGRP, self.mmu_bwa_ref, self.end1, self.end2,
                            self.samtools_tool, self.mmu_samtools_ref, self.sample, log_dir, self.threads)
            if check != 0:
                log(self.loc, date_time() + 'Read filter failure for ' + self.sample + '\n')
                exit(1)
            if not os.path.isfile(self.sample + '.bam'):
                check = bwa_mem_pe(self.bwa_tool, RGRP, self.hsa_bwa_ref, self.end1, self.end2, self.samtools_tool,
                               self.hsa_samtools_ref, self.sample, log_dir, self.threads)
        else:
            log(self.loc, date_time() + 'Starting BWA align\n')
            # check certain key processes
            # skip aligning if bam already exists
            if not os.path.isfile(self.sample + '.bam'):
                check = bwa_mem_pe(self.bwa_tool, RGRP, self.bwa_ref, self.end1, self.end2, self.samtools_tool,
                                   self.samtools_ref, self.sample, log_dir, self.threads)  # rest won't run until completed

            else:
                log(self.loc, date_time() + 'bam file already exists, skipping alignment as well as fastx!\n')

        if check != 0:
            log(self.loc, date_time() + 'BWA failure for ' + self.sample + '\n')
            exit(1)
        # skip sort if sorted file exists already
        log(self.loc, date_time() + 'Sorting BAM file\n')
        if not os.path.isfile(self.sample + '.srt.bam'):
            check = novosort_sort_pe(self.novosort, self.sample, log_dir, self.threads,
                                     self.ram, self.use_nova_flag)  # rest won't run until completed
            if check != 0:
                log(self.loc, date_time() + 'novosort sort failure for ' + self.sample + '\n')
                exit(1)
        else:
            log(self.loc, date_time() + 'Sorted bam file already exists, skipping\n')
        # skip next steps in insert size already calculated
        log(self.loc, date_time() + 'Getting fastq quality score stats\n')
        fastx(self.fastx_tool, self.sample, self.end1, self.end2)  # will run independently of rest of output
        if not os.path.isfile(self.sample + '.insert_metrics.hist') and self.use_nova_flag != 'Y':
            log(self.loc, date_time() + 'Removing PCR duplicates\n')
            picard_rmdup(self.java_tool, self.picard_tool, self.picard_tmp, self.sample, log_dir,
                         self.ram)  # rest won't run until completed
        if self.use_nova_flag == 'Y':
            log(self.loc, date_time() + 'Duplicates removed using novosort.\n')

        log(self.loc, date_time() + 'Gathering SAM flag stats\n')
        flagstats(self.samtools_tool, self.sample)
        log(self.loc, date_time() + 'Calculating insert sizes\n')
        picard_insert_size(self.java_tool, self.picard_tool, self.sample, log_dir, self.ram)  # get insert size metrics.

        # figure out which coverage method to call using seqtype
        log(self.loc, date_time() + 'Calculating coverage for ' + self.seqtype + '\n')
        method = getattr(coverage, (self.seqtype + '_coverage'))
        # run last since this step slowest of the last
        wait_flag = 1
        method(self.bedtools2_tool, self.sample, self.bed_ref, wait_flag)
        log(self.loc, date_time() + 'Checking outputs and uploading results\n')
        # check to see if last expected files have been generated suffix
        self.check_outputs()

        if self.use_nova_flag != 'Y':
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
        obj = self.obj + "/" + self.bid + "/"

        check = upload_to_swift(self.cont, obj)
        if check == 0:
            log(self.loc, date_time() + 'Couchdb successfully updated\n')
            self.status = 0
            log(self.loc,
                date_time() + "Pipeline complete, files successfully uploaded.  Files may be safely removed\n")
        else:
            log(self.loc, date_time() + "All but file upload succeeded\n")
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
