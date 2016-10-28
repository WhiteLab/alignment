#!/usr/bin/env python

import argparse
import os
import sys
import re
from pysam import VariantFile
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.log import log


def create_ind(out):
    res_dict = {}
    mut_out = open(out)
    next(mut_out)
    head = next(mut_out)
    head = head.rstrip('\n').split('\t')
    h_list = (2, 25, 26, 37, 38)
    for line in mut_out:
        info = line.rstrip('\n').split('\t')
        var_tup = '\t'.join((info[0] + info[1] + info[3] + info[4]))
        res_dict[var_tup] = {}
        for i in h_list:
            res_dict[var_tup][head[i]] = info[i]
    mut_out.close()
    return res_dict


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


def output_highest_impact(chrom, pos, ref, alt, ann_list, mut_dict, loc_dict, tflag, out):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''
    var_tup = '\t'.join((chrom , pos , ref , alt))
    (context,  norm_ref_ct, norm_alt_ct, tum_ref_ct, tum_alt_ct) = (mut_dict[var_tup]['context'],
    mut_dict[var_tup]['normal_ref_count'], mut_dict[var_tup]['normal_alt_count'], mut_dict[var_tup]['tumor_ref_count'],
    mut_dict[var_tup]['tumor_alt_count'])
    norm_alt_pct = '0'
    norm_alt_rf = 0.0
    tum_alt_pct = '0'
    tum_alt_rf = 0.0
    tn_ratio = tum_alt_ct
    if int(norm_alt_ct) + int(norm_ref_ct) > 0:
        (norm_alt_rf, norm_alt_pct) = calc_pct(norm_alt_ct, norm_ref_ct)
    if int(tum_alt_ct) + int(tum_ref_ct) > 0:
        (tum_alt_rf, tum_alt_pct) = calc_pct(tum_alt_ct, tum_ref_ct)
    if norm_alt_rf > 0:
        tn_ratio = "{0:.2f}".format(tum_alt_rf/norm_alt_rf)

    for ann in ann_list:
        impact = ann[loc_dict['IMPACT']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            for ann in rank_dict[impact]:
                (gene, effect, aa, codon, snp_id, ExAC_MAF, biotype) = (ann[loc_dict['SYMBOL']],
                ann[loc_dict['Consequence']], ann[loc_dict['Amino_acids']], ann[loc_dict['Codons']],
                ann[loc_dict['Existing_variation']], ann[loc_dict['ExAC_MAF']], ann[loc_dict['BIOTYPE']])
                if f == 0:
                    top_gene = gene
                    f = 1
                    outstring += '\t'.join((chrom, pos, context, ref, alt, norm_ref_ct, norm_alt_ct, norm_alt_pct,
                                            tum_ref_ct, tum_alt_ct, tum_alt_pct, tn_ratio, snp_id, ExAC_MAF, gene,
                                            effect, impact, biotype, codon, aa, tflag)) + '\n'
                out.write(outstring)
                if f == 1 and gene != top_gene and rank != 'MODIFIFER':
                    outstring += '\t'.join((chrom, pos, context, ref, alt, norm_ref_ct, norm_alt_ct, norm_alt_pct,
                                            tum_ref_ct, tum_alt_ct, tum_alt_pct, tn_ratio, snp_id, ExAC_MAF, gene,
                                            effect, impact, biotype, codon, aa, tflag)) + '\n'
                    out.write(outstring)


def gen_report(vcf, out, c):
    # open out file and index counts, context, etc
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    loc = 'LOGS/' + parts[0] + '.vep_priority_report.log'
    log(loc, date_time() + 'Creating prioritized impact reports for ' + vcf + '\n')
    mut_dict = create_ind(out)
    log(loc, date_time() + 'Created index for added mutect info\n')
    on_dict = {}
    if c != 'n':
        on_dict = create_target(c)
        log(loc, date_time() + 'Target file given, creating index for on target info\n')
    vcf_in = VariantFile(vcf)

    out = open(parts[0] + '.vep_prioritized_impact_report.xls', 'w')
    desired = {'Consequence': '', 'IMPACT': '', 'SYMBOL': '', 'Amino_acids': '', 'Codons': '', 'Existing_variation': '',
               'ExAC_MAF': '', 'BIOTYPE': ''}

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
    out.write('chr\tpos\tcontext\tref\talt\tnormal_ref_count\tnormal_alt_count\t%_normal_alt\t'
            'tumor_ref_count\ttumor_alt_count\t%_tumor_alt\tT/N_%_alt_ratio\tsnp_ID\tExAC_MAF\tgene\teffect\timpact'
            '\tbiotype\tcodon_change\tamino_acid_change\ton/off-target\n')
    for record in vcf_in.fetch():
        (chrom, pos, ref, alt) = record.contig, str(record.pos), record.ref, record.alts[0]

        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        tflag = 'NA'
        if c != 'n':
            tflag = mark_target(chrom, pos, on_dict)
        output_highest_impact(chrom, pos, ref, alt, ann_list, mut_dict, desired, tflag, out)

    out.close()
    log(loc, date_time() + 'Creating prioritized report for ' + vcf + ' complete!\n')
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report, prioritizing highest impact calls.')
    parser.add_argument('-v', '--vcf', action='store', dest='vcf',
                        help='VEP annotated variant file')
    parser.add_argument('-o', '--out', action='store', dest='out',
                        help='MuTect output table - has counts for tumor/normal')
    parser.add_argument('-c', '--custom', action='store', dest='c',
                        help='bed file to mark whether hit was on or off-target. if not desired, enter \'n\' ')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.out, args.c)


if __name__ == '__main__':
    main()