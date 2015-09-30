#!/usr/bin/python

import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['snpEff'], config_data['tools']['snpsift'],
            config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['refs']['intervals'])


def custom(config_file, out, ref_mnt, pos, cflag):
    (java, snpeff, snpsift, report, dbsnp, intervals) = parse_config(config_file)
    # create position dict
    fp = open(pos, 'r')
    p_dict = {}
    # track number of entries
    n = 0
    for line in fp:
        line = line.rstrip('\n')
        c = line.split('\t')
        if c[0] not in p_dict:
            p_dict[c[0]] = {}
        p_dict[c[0]][c[1]] = 0

        n += 1
    fp.close()
    # create subfile to annotate
    fc = open('custom.out', 'w')
    fo = open(out, 'r')
    head = next(fo)
    fc.write(head)
    head = next(fo)
    fc.write(head)
    # track found
    x = 0
    for line in fo:
        if x == n:
            break
        line = line.rstrip('\n')
        info = line.split('\t')
        if info[0] in p_dict:
            if info[1] in p_dict[info[0]]:
                fc.write(line + '\n')
                p_dict[info[0]][info[1]] = 1
                x += 1
    fo.close()
    fc.close()
    for chrom in p_dict:
        for p in p_dict[chrom]:
            if p_dict[chrom][p] == 0:
                sys.stderr.write(chrom + '\t' + p + ' not found in out file\n')
    dbsnp = ref_mnt + '/' + dbsnp
    intervals = ref_mnt + '/' + intervals
    run_snpsift = java + ' -jar ' + snpsift + ' annotate ' + dbsnp + ' custom.out > custom.out.sift.vcf'
    subprocess.call(run_snpsift, shell=True)
    run_snpeff = java + ' -jar ' + snpeff + ' eff -t hg19 custom.out.sift.vcf -v > custom.out.eff.vcf'
    subprocess.call(run_snpeff, shell=True)
    run_report = report + ' -i custom.out.eff.vcf -c '
    if cflag == 'n':
        run_report += intervals
    else:
        run_report += 'n'
    run_report += ' > custom.vcf.keep.eff.xls'
    subprocess.call(run_report, shell=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-o', '--out', action='store', dest='out', help='MuTect out file to pasre through')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to mark on/off target')
    parser.add_argument('-p', '--position', action='store', dest='pos', help='File with genomic positions to look for')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, out, ref_mnt, pos, cflag) = (inputs.config_file, inputs.out, inputs.ref_mnt, inputs.pos, inputs.cflag)
    custom(config_file, out, ref_mnt, pos, cflag)
