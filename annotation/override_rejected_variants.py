#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from analysis.variant_annot_pipe import *
import subprocess
import os


def populate_table(table):
    # hash variants to recover
    var_dict = {}
    for line in open(table):
        # format pair\tchrom\tpos\tref\talt\ttype
        (pair, chrom, pos, ref, alt) = line.rstrip('\n').split('\t')
        if pair not in var_dict:
            var_dict[pair] = {}
        ty = 'snv'
        if len(ref) > 1 or len(alt) > 1:
            ty = 'indel'
        if ty not in var_dict[pair]:
            var_dict[pair][ty] = {}
        var_dict[pair][ty]['\t'.join((chrom, pos, ref, alt))] = 0
    return var_dict


def get_snv_files(cont, analysis, annotation, pair):
    src_cmd = '. /home/ubuntu/.novarc;'
    out = analysis + '/' + pair + '/OUTPUT/' + pair + '.out'
    vcf = analysis + '/' + pair + '/OUTPUT/' + pair + '.vcf'
    ann_table = annotation + '/' + pair + '/OUTPUT/' + pair + '.subsitutions.vep.prioritized_impact.report.xls'
    get_out = src_cmd + 'swift download ' + cont + ' ' + out
    subprocess.call(get_out, shell=True)
    get_vcf = src_cmd + 'swift download ' + cont + ' ' + vcf
    subprocess.call(get_vcf, shell=True)
    get_ann_table = src_cmd + 'swift download ' + cont + ' ' + ann_table
    subprocess.call(get_ann_table, shell=True)
    return out, vcf, ann_table


def recreate_analysis(out, vcf, vdict, pair):
    # open both out and vcf files
    out_fh = open(out)
    vcf_fh = open(vcf)
    new_out = pair + '.curated.out.keep'
    new_vcf = pair + '.curated.vcf.keep'
    temp_vcf = pair + '.temp.vcf'
    temp_vcf_fh = open(temp_vcf, 'w')
    temp_out = pair + '.temp.out'
    temp_out_fh = open(temp_out, 'w')
    new_out_fh = open(new_out, 'w')
    new_vcf_fh = open(new_vcf, 'w')
    head = next(out_fh)
    new_out_fh.write(head)
    temp_out_fh.write(head)
    head = next(out_fh)
    new_out_fh.write(head)
    temp_out_fh.write(head)
    for line in vcf_fh:
        new_vcf_fh.write(line)
        temp_vcf_fh.write(line)
        if line[0:6] == '#CHROM':
            break
    for line in vcf_fh:
        out_cur = next(out_fh)
        vcf_info = line.split('\t')
        if vcf_info[6] == 'PASS':
            new_out_fh.write(out_cur)
            new_vcf_fh.write(line)
        else:
            var = '\t'.join((vcf_info[0], vcf_info[1], vcf_info[3], vcf_info[4]))
            if var in vdict:
                vdict[var] = 1
                out_info = out_cur.rstrip('\n').split('\t')
                out_info[-1] = 'OVERRIDE'
                vcf_info[6] = 'OVERRIDE'
                new_out_fh.write('\t'.join(out_info) + '\n')
                temp_out_fh.write('\t'.join(out_info) + '\n')
                new_vcf_fh.write('\t'.join(vcf_info))
                temp_vcf_fh.write('\t'.join(vcf_info))
    out_fh.close()
    vcf_fh.close()
    new_vcf_fh.close()
    new_out_fh.close()
    temp_vcf_fh.close()
    temp_out_fh.close()
    return vdict, temp_vcf


def override_rejected_variants(config_file, table, ref_mnt):
    (novosort, obj, cont, analysis, annotation, germ_flag, indel_flag, annot_used, wg) = parse_config(config_file)
    sys.stderr.write(date_time() + 'Populating update table\n')
    var_dict = populate_table(table)
    pair_list = 'pairs.txt'
    pair_fh = open(pair_list, 'w')
    ann_table_list = []
    for pair in var_dict:
        sys.stderr.write(date_time() + 'Processing ' + pair + '\n')
        pair_fh.write(pair + '\n')
        if 'snv' in var_dict[pair]:
            (out, vcf, ann_table) = get_snv_files(cont, analysis, annotation, pair)
            ann_table_list.append(ann_table)
            (var_dict[pair]['snv'], temp_vcf) = recreate_analysis(out, vcf, var_dict[pair]['snv'], pair)
    sys.stderr.write(date_time() + 'Annotating recovered variants\n')
    check = vep(config_file, pair_list, ref_mnt, '.temp.vcf', '.snv.curated.vcf', '.temp.out', 'mutect')
    if check == 0:
        sys.stderr.write(date_time() + 'Running VEP failed.  Check parameters\n')
        exit(1)
    sys.stderr.write(date_time() + 'Recombining tables\n')
    # combine new entries into old reports
    for ann_table in ann_table_list:
        f = 0
        fn = os.path.basename(ann_table)
        parts = fn.split('.')
        old = open(ann_table)
        to_integrate = parts[0] + '.subsitutions.vep.prioritized_impact.report.xls'
        update = open(to_integrate)
        new_table = parts[0] + '.subsitutions.vep.curated_reports.xls'
        nt = open(new_table, 'w')
        head = next(old)
        nt.write(head)
        next(update)
        new = next(update)
        new_info = new.split('\t')
        for line in old:
            info = line.split('\t')
            if f == 0 and new_info[0] == info[0] and int(info[1]) > int(new_info[1]):
                nt.write(new)
                try:
                    new = next(update)
                    new_info = new.split('\t')
                except:
                    f = 1
            nt.write(line)
        if f == 0:
            nt.write(new)
        nt.close()
    sys.stderr.write('Fin!\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Recover rejected variants deemed to be real by manual curation and '
                                                 'reincoprorate into reports.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-t', '--table', action='store', dest='table', help='Table with pairs and coords to recover')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, table, ref_mnt) = (inputs.config_file, inputs.table, inputs.ref_mnt)
    override_rejected_variants(config_file, table, ref_mnt)
