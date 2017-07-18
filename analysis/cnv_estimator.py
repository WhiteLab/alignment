#!/usr/bin/python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
import re
import subprocess
from utility.job_manager import job_manager
import json


def get_genes(bed):
    temp = {}
    glist = []
    for line in open(bed):
        info = line.rstrip('\n').split('\t')
        m = re.match('(\S+)_chr.*', info[-1])
        if m.group(1) not in temp:
            temp[m.group(1)] = 1
            glist.append(m.group(1))
    return glist


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['cont'], config_data['refs']['obj'], config_data['tools']['samtools'], \
           config_data['refs']['analysis'], config_data['refs']['annotation'], config_data['refs']['capture']


def get_bam_name(bnid, src_cmd, cont, obj):
    list_cmd = src_cmd + 'swift list ' + cont + ' --prefix ' + obj + '/' + bnid + '/BAM/'
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    bam = ''
    bai = ''
    dl_cmd = ''
    for fn in re.findall('(.*)\n', flist):
        # depending on software used to index, .bai extension may follow bam
        # test = re.match('^\S+_\w*\d+\.rmdup.srt.ba[m|i]$', fn) or re.match('^\S+_\w*\d+\.rmdup.srt.bam.bai$', fn)
        test = re.search(bnid + '.merged.final.ba[m|i]$', fn) or re.search(bnid + '.merged.final.bam.bai$', fn)
        if test:
            dl_cmd += src_cmd + 'swift download ' + cont + ' ' + fn + ';'
            if fn[-3:] == 'bam':
                bam = fn
            else:
                bai = fn
    return dl_cmd, bam, bai


def cnv_pipe(config_file, sample_pairs, ref_mnt):
    src_cmd = '. /home/ubuntu/.novarc;'
    job_list = []
    (cont, obj, samtools, ana, ann, bed) = parse_config(config_file)
    bed = ref_mnt + '/' + bed
    bed_t1 = bed.replace('.bed', '_t1.bed')
    bed_t2 = bed.replace('.bed', '_t2.bed')
    pair_list = []
    t1_genes = get_genes(bed_t1)
    t2_genes = get_genes(bed_t2)
    # calc coverage for all gene capture regions
    for pairs in open(sample_pairs):
        pair_set = pairs.rstrip('\n').split('\t')
        pair_list.append('\t'.join((pair_set[1], pair_set[2])))
        (tum_dl_cmd, tum_bam, tum_bai) = get_bam_name(pair_set[1], src_cmd, cont, obj)
        (norm_dl_cmd, norm_bam, norm_bai) = get_bam_name(pair_set[2], src_cmd, cont, obj)
        job_list.append(tum_dl_cmd)
        job_list.append(norm_dl_cmd)
        job_manager(job_list, '8')
        job_list = []
        job_list.append(samtools + ' bedcov ' + bed_t1 + ' ' + tum_bam + ' > ' + pair_set[1] + '.t1.sam.bedcov.txt;')
        job_list.append(samtools + ' bedcov ' + bed_t2 + ' ' + tum_bam + ' > ' + pair_set[1] + '.t2.sam.bedcov.txt;')
        job_list.append(samtools + ' bedcov ' + bed_t1 + ' ' + norm_bam + ' > ' + pair_set[2] + '.t1.sam.bedcov.txt;')
        job_list.append(samtools + ' bedcov ' + bed_t2 + ' ' + norm_bam + ' > ' + pair_set[2] + '.t2.sam.bedcov.txt;')
        job_manager(job_list, '8')
        cleanup = 'rm ' + tum_bam + ' ' + tum_bai + ' ' + norm_bam + ' ' + norm_bai + ';'
        subprocess.call(cleanup, shell=True)
    # process coverage files, assess cnv


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pipeline for variant calls and annotation using mutect and snpEff')
    parser.add_argument('-sp', '--sample-pairs', action='store', dest='sample_pairs',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-r', '--reference', action='store', dest='ref_mnt',
                        help='Directory references are mounted, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_pairs, config_file, ref_mnt) = (inputs.sample_pairs, inputs.config_file, inputs.ref_mnt)
    cnv_pipe(config_file, sample_pairs, ref_mnt)