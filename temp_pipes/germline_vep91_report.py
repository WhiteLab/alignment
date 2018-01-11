#!/usr/bin/env python3

import re
from pysam import VariantFile
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from annotation.report_tools import create_index


def output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, ann_list, loc_dict, out, ref_flag):
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
                (gene, tx_id, hgvsg, variant_class, effect, aa_pos, aa, codon, snp_id, ExAC_MAF, biotype, sift,
                 clin_sig) = (ann[loc_dict['SYMBOL']], ann[loc_dict['Feature']], ann[loc_dict['HGVSg']],
                 ann[loc_dict['VARIANT_CLASS']], ann[loc_dict['Consequence']], ann[loc_dict['Protein_position']],
                 ann[loc_dict['Amino_acids']], ann[loc_dict['Codons']], ann[loc_dict['Existing_variation']],
                 ann[loc_dict['gnomAD_AF']], ann[loc_dict['BIOTYPE']],
                 ann[loc_dict['SIFT']], ann[loc_dict['CLIN_SIG']])
                phred = ann[loc_dict['CADD_PHRED']][0]
                if variant_class != 'SNV':
                    phred = ann[loc_dict['CADD_PHRED'][1]]
                # Format amino acid change to be oldPOSnew
                if len(aa) > 0:
                    # if a snv or syn, just aaPOS
                    test = aa.split('/')
                    if len(test) == 1:
                        aa += str(aa_pos)
                    else:
                        aa = test[0] + str(aa_pos) + test[1]
                # need to parse exac maf to get desired allele freq, not all possible
                alt_ct = alt_ct.replace('(', '')
                alt_ct = alt_ct.replace(')', '')
                alt_ct = alt_ct.split(',')[0]
                cur_var = '\t'.join((chrom, pos, ref, alt, alt_ct, tot_ct, gene, hgvsg, tx_id, effect, impact,
                                     biotype, codon, aa, snp_id, variant_class, sift, ExAC_MAF, clin_sig, phred)) + '\n'
                if ref_flag == 'n':
                    if f == 0:
                        top_gene = gene
                        f = 1
                        outstring += cur_var
                    if f == 1 and gene != top_gene and impact != 'MODIFIER' and f1 != 0:
                        outstring += cur_var
                        f1 = 1
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


def gen_report(vcf, sample, ref_flag):
    vcf_in = VariantFile(vcf)
    # run cadd twice over snv and indel file
    out = open(sample + '.germline.vep91.xls', 'w')
    desired = {'Consequence': 0, 'IMPACT': 0, 'SYMBOL': 0, 'Feature': 0, 'HGVSg': 0, 'Protein_position': 0,
               'Amino_acids': 0, 'Codons': 0, 'BIOTYPE': 0, 'SIFT': 0, 'Existing_variation': 0, 'VARIANT_CLASS': 0,
               'gnomAD_AF': 0, 'CLIN_SIG': 0, 'CADD_PHRED': []}

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
            if desc_list[i] == 'CADD_PHRED':
                desired[desc_list[i]].append(i)
            else:
                desired[desc_list[i]] = i
    out.write('CHROM\tPOS\tREF\tAllele\tTotal Allele Count\tTotal Position Coverage\tGene\tHGVSg\tTranscript_id'
              '\tEffect\tIMPACT\tBIOTYPE\tCodons\tAmino_acids\tExisting_variation\tVARIANT_CLASS\tSIFT\tgnomAD_AF'
              '\tCLIN_SIG\tCADD_PHRED\n')
    if ref_flag != 'n':
        ref_flag = create_index(ref_flag)
    for record in vcf_in.fetch():
        (chrom, pos, ref, alt, alt_ct, tot_ct) = (record.contig, str(record.pos), record.ref, record.alts[0],
                                                  str(record.info['TR']), str(record.info['TC']))
        ann_list = [_.split('|') for _ in record.info['ANN']]
        output_highest_impact(chrom, pos, ref, alt, alt_ct, tot_ct, ann_list, desired, out, ref_flag)
    out.close()
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='parse VEP annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='VEP annotated variant file')
    parser.add_argument('-s', '--sample', action='store', dest='sample', help='Sample name')
    parser.add_argument('-r', '--reference', action='store', dest='ref', help='Tab-separated reference table with gene '
                    'symbols and refseq + ensembl ids to standardize what transcript is used.  Use flag \'n\' to skip')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (vcf, sample, ref) = (inputs.infile, inputs.sample, inputs.ref)
    gen_report(vcf, sample, ref)
