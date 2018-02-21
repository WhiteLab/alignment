#!/usr/bin/env python3

import argparse
import os
import sys
from pysam import VariantFile
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.log import log
from annotation.report_tools import *


def create_mutect_ind(out):
    res_dict = {}
    mut_out = open(out)
    next(mut_out)
    head = next(mut_out)
    head = head.rstrip('\n').split('\t')
    h_list = (2, 25, 26, 37, 38)
    for line in mut_out:
        info = line.rstrip('\n').split('\t')
        var_tup = '\t'.join((info[0], info[1], info[3], info[4]))
        res_dict[var_tup] = {}
        for i in h_list:
            res_dict[var_tup][head[i]] = info[i]
    mut_out.close()
    return res_dict


def output_highest_impact(chrom, pos, ref, alt, ann_list, mut_dict, loc_dict, tflag, out, ref_flag):
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
    var_tup = '\t'.join((chrom, pos, ref, alt))
    (context,  norm_ref_ct, norm_alt_ct, tum_ref_ct, tum_alt_ct) = (mut_dict[var_tup]['context'],
    mut_dict[var_tup]['n_ref_count'], mut_dict[var_tup]['n_alt_count'], mut_dict[var_tup]['t_ref_count'],
    mut_dict[var_tup]['t_alt_count'])
    norm_alt_pct = '0'
    norm_alt_rf = 0.0
    tum_alt_pct = '0'
    tum_alt_rf = 0.0
    tn_ratio = tum_alt_ct
    if int(norm_alt_ct) + int(norm_ref_ct) > 0:
        (norm_alt_rf, norm_alt_pct) = calc_pct(norm_ref_ct, norm_alt_ct)
    if int(tum_alt_ct) + int(tum_ref_ct) > 0:
        (tum_alt_rf, tum_alt_pct) = calc_pct(tum_ref_ct, tum_alt_ct)
    if norm_alt_rf > 0:
        tn_ratio = "{0:.2f}".format(tum_alt_rf/norm_alt_rf)

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
                (gene, tx_id, hgvsg, effect, aa_pos, aa, codon, snp_id, gnomAD_AF, biotype) = (ann[loc_dict['SYMBOL']],
                ann[loc_dict['Feature']], ann[loc_dict['HGVSg']], ann[loc_dict['Consequence']],
                ann[loc_dict['Protein_position']], ann[loc_dict['Amino_acids']], ann[loc_dict['Codons']],
                ann[loc_dict['Existing_variation']], ann[loc_dict['gnomAD_AF']], ann[loc_dict['BIOTYPE']])

                # Format amino acid change to be oldPOSnew
                if len(aa) > 0:
                    # if a snv or syn, just aaPOS
                    test = aa.split('/')
                    if len(test) == 1:
                        aa += str(aa_pos)
                    else:
                        aa = test[0] + str(aa_pos) + test[1]
                cur_var = '\t'.join((chrom, pos, context, ref, alt, norm_ref_ct, norm_alt_ct, norm_alt_pct, tum_ref_ct,
                                     tum_alt_ct, tum_alt_pct, tn_ratio, snp_id, gnomAD_AF, gene, hgvsg, tx_id, effect,
                                     impact, biotype, codon, aa, tflag)) + '\n'
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


def gen_report(vcf, out, c, ref_flag, cache):
    # open out file and index counts, context, etc
    fn = os.path.basename(vcf)
    parts = fn.split('.')
    sample = parts[0]
    loc = 'LOGS/' + sample + '.subsitutions.vep' + cache + '.priority_report.log'
    suffix = '.subsitutions.vep' + cache + '.prioritized_impact.report.xls'
    log(loc, date_time() + 'Creating prioritized impact reports for ' + vcf + '\n')
    mut_dict = create_mutect_ind(out)
    log(loc, date_time() + 'Created index for added mutect info\n')
    on_dict = {}
    if c != 'n':
        on_dict = create_target(c)
        log(loc, date_time() + 'Target file given, creating index for on target info\n')
    vcf_in = VariantFile(vcf)
    out_fn = sample + suffix
    out = open(out_fn, 'w')
    desired = {'Consequence': 0, 'IMPACT': 0, 'SYMBOL': 0, 'Feature': 0, 'HGVSg': 0, 'Protein_position': 0,
               'Amino_acids': 0, 'Codons': 0, 'Existing_variation': 0, 'ExAC_MAF': 0, 'BIOTYPE': 0}

    desc_string = vcf_in.header.info['ANN'].record['Description']
    desc_string = desc_string.lstrip('"')
    desc_string = desc_string.rstrip('"')
    desc_string = desc_string.replace('Consequence annotations from Ensembl VEP. Format: ', '')
    f_pos_list = []
    desc_list = desc_string.split('|')
    ann_size = len(desc_list)
    for i in range(0, ann_size, 1):
        if desc_list[i] in desired:
            f_pos_list.append(i)
            desired[desc_list[i]] = i
    out.write('chr\tpos\tcontext\tref\talt\tnormal_ref_count\tnormal_alt_count\t%_normal_alt\ttumor_ref_count\t'
              'tumor_alt_count\t%_tumor_alt\tT/N_%_alt_ratio\tsnp_ID\tExAC_MAF\tgene\tHGVSg\ttx_id\teffect\timpact\t'
              'biotype\tcodon_change\tamino_acid_change\ton/off-target\n')
    if ref_flag != 'n':
        ref_flag = create_index(ref_flag)
    for record in vcf_in.fetch():
        (chrom, pos, ref, alt) = record.contig, str(record.pos), record.ref, record.alts[0]
        ann_list = [_.split('|') for _ in record.info['ANN']]
        tflag = 'NA'
        if c != 'n':
            tflag = mark_target(chrom, pos, on_dict)
            # only outputting ON TARGET hits
            if tflag == 'OFF':
                continue
        output_highest_impact(chrom, pos, ref, alt, ann_list, mut_dict, desired, tflag, out, ref_flag)

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
    parser.add_argument('-r', '--reference', action='store', dest='ref', help='Tab-separated reference table with gene '
                    'symbols and refseq + ensembl ids to standardize what transcript is used.  Use flag \'n\' to skip')
    parser.add_argument('-n', '--cache', action='store', dest='cache', help='vep version, i.e. 91')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    gen_report(args.vcf, args.out, args.c, args.ref, args.cache)


if __name__ == '__main__':
    main()
