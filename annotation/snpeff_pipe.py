#!/usr/bin/python
import sys
from utility.date_time import date_time
from utility.job_manager import job_manager
import subprocess
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return (config_data['tools']['java'], config_data['tools']['snpEff'], config_data['tools']['snpsift'],
            config_data['tools']['report'], config_data['refs']['dbsnp'], config_data['refs']['intervals'])


def snpeff_pipe(config_file, sample_pairs, ref_mnt, cflag):
    # edit to grab from config max thread count
    max_t = 8
    (java, snpeff, snpsift, report, dbsnp, intervals) = parse_config(config_file)
    dbsnp = ref_mnt + '/' + dbsnp
    intervals = ref_mnt + '/' + intervals
    fh = open(sample_pairs)
    mk_log_dir = 'mkdir LOGS'
    subprocess.call(mk_log_dir, shell=True)
    cmd_list = []
    run_snpsift = java + ' -jar ' + snpsift + ' annotate ' + dbsnp
    run_snpeff = java + ' -jar ' + snpeff + ' eff -t hg19 '
    for line in fh:
        line = line.rstrip('\n')
        (sample, tumor_id, normal_id) = line.split('\t')
        # run snpsift first, then snpeff
        run_report = report + ' -i ' + sample + '.out.keep.eff.vcf -c '
        if cflag == 'n':
            run_report += intervals
        else:
            run_report += 'n'
        run_report += ' > ' + sample + '.vcf.keep.eff.xls'
        run_snp = run_snpsift + ' ' + sample + '.out.keep > ' + sample + '.out.keep.sift.vcf 2> LOGS/' + sample \
                  + '.snpeff.log;' + run_snpeff + ' ' + sample + '.out.keep.sift.vcf -v > ' + sample \
                  + '.out.keep.eff.vcf  2>> LOGS/' + sample + '.snpeff.log;' + run_report
        cmd_list.append(run_snp)
    job_manager(cmd_list, max_t)
    sys.stderr.write(date_time() + 'SNP annotation  completed!\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='SNP annotation.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-f', '--flag', action='store', dest='cflag',
                        help='\'y\' if whole genome, \'n\' if custom capture to mark on/off target')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt, cflag) = (
        inputs.config_file, inputs.sample_pairs, inputs.ref_mnt, inputs.cflag)
    snpeff_pipe(config_file, sample_pairs, ref_mnt, cflag)
