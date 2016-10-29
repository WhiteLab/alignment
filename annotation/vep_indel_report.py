#!/usr/bin/env python

import argparse
import os
import sys
import re
from pysam import VariantFile
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.log import log


def create_target(cfile):
        cindex = {}
        fh = open(cfile, 'r')
        for line in fh:
            m = re.search('(chr\w+):(\d+)-(\d+)', line)
            try:
                (chrom, start, end) = (m.group(1), m.group(2), m.group(3))
            except:
                sys.stderr.write(line + ' doesn\'t fit format (chr\w+):(\d+)-(\d+), skipping\n')
                continue
            if chrom not in cindex:
                cindex[chrom] = {}
            cindex[chrom][(int(start))] = int(end)
        return cindex


def mark_target(chrom, pos, on_dict):
    f = 0
    if chrom in on_dict:
        for start in sorted(on_dict[chrom]):
            if start <= int(pos) <= on_dict[chrom][start]:
                f = 1
                break
            elif start > int(pos):
                break
    status = "OFF"
    if f == 1:
        status = "ON"
    return status


def calc_pct(a, b):
    # return both formatted and unformatted
    ratio = float(b) / (float(a) + float(b)) * 100
    fmt = "{0:.2f}%".format(ratio)
    return ratio, fmt


def output_highest_impact(chrom, pos, ref, alt, ann_list, loc_dict, tflag, out):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''

    for ann in ann_list:
        impact = ann[loc_dict['IMPACT']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            for ann in rank_dict[impact]:
                # need to add coverage info for indels
                (gene, variant_class, effect, aa, codon, snp_id, ExAC_MAFs, biotype) = (ann[loc_dict['SYMBOL']],
                ann[loc_dict['VARIANT_CLASS']], ann[loc_dict['Consequence']], ann[loc_dict['Amino_acids']],
                ann[loc_dict['Codons']], ann[loc_dict['Existing_variation']], ann[loc_dict['ExAC_MAF']],
                ann[loc_dict['BIOTYPE']])
                # need to parse exac maf to get desired allele freq, not all possible
                ExAC_MAF = ''
                if len(ExAC_MAFs) > 1:
                    maf_list = ExAC_MAFs.split('&')
                    for maf in maf_list:
                        check = re.match(alt + ':(\S+)', maf)
                        if check:
                            ExAC_MAF = check.group(1)
                if f == 0:
                    top_gene = gene
                    f = 1
                    outstring += '\t'.join((chrom, pos, ref, alt, snp_id, ExAC_MAF, gene, variant_class,
                                            effect, impact, biotype, codon, aa, tflag)) + '\n'
                    out.write(outstring)
                if f == 1 and gene != top_gene and impact != 'MODIFIER':
                    outstring += '\t'.join((chrom, pos, ref, alt, snp_id, ExAC_MAF, gene,
                                            effect, impact, biotype, codon, aa, tflag)) + '\n'
                    out.write(outstring)


def gen_report(vcf, c):
    # open out file and index counts, context, etc
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    loc = 'LOGS/' + parts[0] + '.indels.vep_priority.report.log'
    log(loc, date_time() + 'Creating prioritized impact reports for ' + vcf + '\n')
    on_dict = {}
    if c != 'n':
        on_dict = create_target(c)
        log(loc, date_time() + 'Target file given, creating index for on target info\n')
    vcf_in = VariantFile(vcf)

    out = open(parts[0] + '.indels.vep.prioritized_impact.report.xls', 'w')
    desired = {'Consequence': '', 'IMPACT': '', 'SYMBOL': '', 'Amino_acids': '', 'Codons': '', 'Existing_variation': '',
               'ExAC_MAF': '', 'BIOTYPE': '', 'VARIANT_CLASS': ''}

    desc_string = vcf_in.header.info['ANN'].record['Description']
    desc_string = desc_string.lstrip('"')
    desc_string = desc_string.rstrip('"')
    desc_string = desc_string.replace('Consequence annotations from Ensembl VEP. Format: ', '')
    f_pos_list = []
    desc_list = desc_string.split('|')
    ann_size = len(desc_list)
    for i in xrange(0, ann_size, 1):
        if desc_list[i] in desired:
            f_pos_list.append(i)
            desired[desc_list[i]] = i
    out.write('chr\tpos\tcontext\tref\talt\tsnp_ID\tExAC_MAF\tgene\tvariant_class_effect\timpact'
            '\tbiotype\tcodon_change\tamino_acid_change\ton/off-target\n')
    for record in vcf_in.fetch():
        (chrom, pos, ref, alt) = record.contig, str(record.pos), record.ref, record.alts[0]
        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        tflag = 'NA'
        if c != 'n':
            tflag = mark_target(chrom, pos, on_dict)
        output_highest_impact(chrom, pos, ref, alt, ann_list, desired, tflag, out)

    out.close()
    log(loc, date_time() + 'Creating prioritized report for ' + vcf + ' complete!\n')
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report, prioritizing highest impact calls.')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='VEP annotated variant file')
    parser.add_argument('-c', '--custom', action='store', dest='c',
                        help='bed file to mark whether hit was on or off-target. if not desired, enter \'n\' ')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.c)


if __name__ == '__main__':
    main()
