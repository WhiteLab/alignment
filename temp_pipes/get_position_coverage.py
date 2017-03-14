#!/usr/bin/env python

import sys
import os
import subprocess
import re
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.job_manager import job_manager
from utility.date_time import date_time


def create_bed(line, bed):
    b_dict = {}
    var = line.rstrip('\n').split('\t')
    m = re.search('\S+-(chr\w+)_(\d+)_\w+->\w+', var[0])
    (chrom, pos) = (m.group(1), m.group(2))
    if chrom not in b_dict:
        b_dict[chrom] = {}
    if pos not in b_dict[chrom]:
        bed.write(chrom + '\t' + str(int(pos)-1) + '\t' + pos + '\n')
        b_dict[chrom][pos] = 1
    return var[0]


def build_jobs(samtools, bed, sample):
    prefix = 'ALIGN/' + sample + '/BAM/'
    cmd = samtools + ' depth -Q 20 -b ' + bed + prefix + sample + 'merged.final.bam > ' + sample + '_covered.txt;'
    return cmd


def compile_results(slist):
    cov_dict = {}
    for i in xrange(1, len(slist), 1):
        for line in open(slist[i] + '_covered.txt'):
            info = line.rstrip('\n').split('\t')
            if info[0] not in cov_dict:
                cov_dict[info[0]] = {}
            if info[1] not in cov_dict[info[0]]:
                cov_dict[info[0]][info[1]] = {}
            cov_dict[info[0]][info[1]][slist[i]] = info[2]
    return cov_dict


def calc_pos_cov(table, samtools, out):
    fh = open(table)
    head = next(fh)
    head = head.rstrip('\n').split('\t')
    # create bed file to get coverage
    bed_fn = out + '.bed'
    bed = open(bed_fn, 'w')
    vlist = []
    # in the event an indel happens in one sample and snv in another at same position, don't process twice

    for line in fh:
        vlist.append(create_bed(line, bed))
    bed.close()
    fh.close()
    job_list = []
    src_cmd = '. ~/.novarc;'
    # get bams, then build jobs
    for i in xrange(1, len(head), 1):
        sys.stderr.write(date_time() + 'Getting bam for ' + head[i] + '\n')
        bam = 'ALIGN/' + head[i] + '/BAM/' + head[i] + '.merged.final.bam'
        dl_cmd = src_cmd + 'swift download PDX --prefix ALIGN/' + head[i] + '/BAM/' + head[i] + '.merged.final.ba;'
        subprocess.call(dl_cmd, shell=True)
        # try pdx container, if not, try pancan
        if os.path.isfile(bam):
            job_list.append(build_jobs(samtools, bed_fn, head[i]))
        else:
            sys.stderr.write(date_time() + 'Bam for sample ' + head[i] + ' not in PDX contanier, trying PANCAN\n')
            dl_cmd = src_cmd + 'swift download PANCAN --prefix ALIGN/' + head[i] + '/BAM/' + head[i] \
                     + '.merged.final.ba;'
            subprocess.call(dl_cmd, shell=True)
            if os.path.isfile(bam):
                job_list.append(build_jobs(samtools, bed_fn, head[i]))
            else:
                sys.stderr.write('Could not find bam for ' + head[i] + '\n')
                exit(1)
    sys.stderr.write('Running depth jobs\n')
    job_manager(job_list, '8')
    sys.stderr.write(date_time() + 'Compiling results\n')
    cov_dict = compile_results(head)
    sys.stderr.write(date_time() + 'Writing to output table\n')
    out_fh = open(out + '_variant_coverage_table.txt')
    out_fh.write('\t'.join(head) + '\n')
    for var in vlist:
        out_fh.write(var)
        for i in xrange(1, len(head), 1):
            m = re.search('\S+-(chr\w+)_(\d+)_\w+->\w+', var)
            (chrom, pos) = (m.group(1), m.group(2))
            if head[i] in cov_dict[chrom][pos]:
                out_fh.write('\t' + cov_dict[chrom][pos][head[i]])
            else:
                out_fh.write('\t0')
        out_fh.write('\n')
    out_fh.close()
    sys.stderr.write(date_time() + 'Fin\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Simple total read coverage of a position using a variant table.')
    parser.add_argument('-t', '--table', action='store', dest='table', help='Variant table')
    parser.add_argument('-s', '--samtools', action='store', dest='samtools', help='Location of samtools binary')
    parser.add_argument('-o', '--out', action='store', dest='out', help='Output prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (table, samtools, out) = (inputs.table, inputs.samtools, inputs.out)
    calc_pos_cov(table, samtools, out)
