#!/usr/bin/env python

import argparse
import os
import sys
import re
from pysam import VariantFile
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from utility.log import log
from report_tools import *
import pdb


def calc_pct(a, b):
    # return both formatted and unformatted
    ratio = float(b) / (float(a) + float(b)) * 100
    fmt = "{0:.2f}%".format(ratio)
    return ratio, fmt


def output_highest_impact(chrom, pos, ref, alt, not_shared, ann_list, loc_dict, tflag, out, ref_flag, call_type):
    rank = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
    top_gene = ''
    f = 0
    # secondary hit flag to avoid excessive repeats
    f1 = 0
    # index annotations by impact rank
    rank_dict = {}
    outstring = ''
    cand_top = []
    cand_next = []
    r_flag = 99
    for ann in ann_list:
        impact = ann[loc_dict['IMPACT']]
        if impact not in rank_dict:
            rank_dict[impact] = []
        rank_dict[impact].append(ann)
    for impact in rank:
        if impact in rank_dict:
            cur_rank = rank.index(impact)
            if cur_rank < r_flag:
                r_flag = cur_rank
            for ann in rank_dict[impact]:
                # need to add coverage info for indels
                (gene, tx_id, variant_class, effect, aa_pos, aa, codon, snp_id, ExAC_MAFs, biotype) = \
                    (ann[loc_dict['SYMBOL']], ann[loc_dict['Feature']], ann[loc_dict['VARIANT_CLASS']],
                     ann[loc_dict['Consequence']], ann[loc_dict['Protein_position']], ann[loc_dict['Amino_acids']],
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
                # Format amino acid change to be oldPOSnew
                if len(aa) > 0:
                    # if a snv or syn, just aaPOS
                    test = aa.split('/')
                    if len(test) == 1:
                        aa += str(aa_pos)
                    else:
                        aa = test[0] + str(aa_pos) + test[1]
                if call_type == 'snv':
                    norm_alt_pct = '0'
                    norm_alt_rf = 0.0
                    tum_alt_pct = '0'
                    tum_alt_rf = 0.0
                    (norm_ref_ct, norm_alt_ct, tum_alt_ct, tum_ref_ct) = (not_shared['norm_ref_ct'],
                                        not_shared['norm_alt_ct'], not_shared['tum_alt_ct'], not_shared['tum_ref_ct'])
                    tn_ratio = tum_alt_ct
                    if int(norm_alt_ct) + int(norm_ref_ct) > 0:
                        (norm_alt_rf, norm_alt_pct) = calc_pct(norm_ref_ct, norm_alt_ct)
                    if int(tum_alt_ct) + int(tum_ref_ct) > 0:
                        (tum_alt_rf, tum_alt_pct) = calc_pct(tum_ref_ct, tum_alt_ct)
                    if norm_alt_rf > 0:
                        tn_ratio = "{0:.2f}".format(tum_alt_rf / norm_alt_rf)
                    try:
                        cur_var = '\t'.join((chrom, pos, ref, alt, str(norm_ref_ct), str(norm_alt_ct), norm_alt_pct,
                        str(tum_ref_ct), str(tum_alt_ct), tum_alt_pct, str(tn_ratio), snp_id, ExAC_MAF, gene, tx_id,
                        variant_class, effect, impact, biotype, codon, aa, tflag)) + '\n'
                    except:
                        pdb.set_trace()
                else:
                    try:
                        cur_var = '\t'.join((chrom, pos, ref, alt, snp_id, ExAC_MAF, gene, tx_id, variant_class, effect,
                                         impact, biotype, codon, aa, tflag)) + '\n'
                    except:
                        pdb.set_trace()

                if ref_flag == 'n':
                    if f == 0:
                        top_gene = gene
                        f = 1
                        outstring += cur_var
                    if f == 1 and gene != top_gene and impact != 'MODIFIER' and f1 != 0:
                        outstring += cur_var
                        f1 = 0
                else:
                    if f < 1 and cur_rank == r_flag:
                        if gene in ref_flag and tx_id == ref_flag[gene]:
                            top_gene = gene
                            outstring += cur_var
                            f = 1
                        else:
                            top_gene = gene
                            f = 0.5
                            cand_top.append(cur_var)

                    if f > 0 and gene != top_gene and impact != 'MODIFIER' and f1 < 1 and cur_rank == (r_flag + 1):
                        if gene in ref_flag and tx_id == ref_flag[gene]:
                            outstring += cur_var
                            f1 = 1
                        else:
                            f1 = 0.5
                            cand_next.append(cur_var)
    if ref_flag != 'n':
        if f == 0.5:
            outstring += cand_top[0]
        if f1 == 0.5:
            outstring += cand_next[0]
    out.write(outstring)


def gen_report(vcf, c, ref_flag):
    # open out file and index counts, context, etc
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    loc = 'LOGS/' + parts[0] + '.snv.strelka.vep_priority.report.log'
    log(loc, date_time() + 'Creating prioritized impact reports for ' + vcf + '\n')
    on_dict = {}
    if c != 'n':
        on_dict = create_target(c)
        log(loc, date_time() + 'Target file given, creating index for on target info\n')
    vcf_in = VariantFile(vcf)
    call_type = 'snv'
    if bool(re.search('indel', fn)):
        out = open(parts[0] + '.indel.strelka.vep.prioritized_impact.report.xls', 'w')
        call_type = 'indel'
    else:
        out = open(parts[0] + '.snv.strelka.vep.prioritized_impact.report.xls', 'w')
    desired = {'Consequence': 0, 'IMPACT': 0, 'SYMBOL': 0, 'Feature': 0, 'Protein_position': 0, 'Amino_acids': 0,
               'Codons': 0, 'Existing_variation': 0, 'ExAC_MAF': 0, 'BIOTYPE': 0, 'VARIANT_CLASS': 0}

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
    if call_type == 'snv':
        out.write('chr\tpos\tref\talt\tnormal_ref_count\tnormal_alt_count\t%_normal_alt\ttumor_ref_count\t'
              'tumor_alt_count\t%_tumor_alt\tT/N_%_alt_ratio\tsnp_ID\tExAC_MAF\tgene\ttranscript_id\t'
              'variant_class_effect\teffect\timpact\tbiotype\tcodon_change\tamino_acid_change\ton/off-target\n')
    else:
        out.write('chr\tpos\tref\talt\tsnp_ID\tExAC_MAF\tgene\ttranscript_id\tvariant_class_effect\teffect\timpact\t'
                  'biotype\tcodon_change\tamino_acid_change\ton/off-target\n')
    if ref_flag != 'n':
        ref_flag = create_index(ref_flag)

    for record in vcf_in.fetch():
        # dict contains what's different between strelka indel and snv reports
        (chrom, pos, ref, alt) = (record.contig, str(record.pos),
        record.ref, record.alts[0])
        if call_type == 'snv':
            not_shared = {'norm_ref_ct': record.samples['NORMAL'][(record.ref + 'U')][0],
                          'norm_alt_ct': record.samples['TUMOR'][(record.ref + 'U')][0],
                           'tum_ref_ct': record.samples['NORMAL'][(record.alts[0] + 'U')][0],
                           'tum_alt_ct': record.samples['TUMOR'][(record.alts[0] + 'U')][0]}
        else:
            not_shared = {}
        ann_list = [_.split('|') for _ in record.info['ANN'].split(',')]
        tflag = 'NA'
        if c != 'n':
            tflag = mark_target(chrom, pos, on_dict)
            # only outputting ON TARGET hits
            if tflag == 'OFF':
                continue
        output_highest_impact(chrom, pos, ref, alt, not_shared, ann_list, desired, tflag, out, ref_flag, call_type)

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
    parser.add_argument('-r', '--reference', action='store', dest='ref', help='Tab-separated reference table with gene '
                    'symbols and refseq + ensembl ids to standardize what transcript is used. Use flag \'n\' to skip')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.c, args.ref)


if __name__ == '__main__':
    main()
