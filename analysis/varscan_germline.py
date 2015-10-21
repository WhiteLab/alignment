#!/usr/bin/env python
__author__ = 'Miguel'
import sys
import json

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['samtools'], config_data['tools']['varscan'], config_data['refs']['fai'], \
           config_data['refs']['fa_ordered'], config_data['params']['threads']


def varscan_germline(config_file, sample, ref_mnt):
    (samtools, varscan, region, fasta, th) = parse_config(config_file)
    region = ref_mnt + '/' + region
    fasta = ref_mnt + '/' + fasta
    rf = open(region, 'r')
    cmd_list = []
    for line in rf:
        chrom = line.split('\t')
        cmd = samtools + ' mpileup -r ' + chrom[0] + ' -B -f ' + fasta + ' ' + sample +\
              '.merged.final.bam | java -Xmx4000m -jar ' + varscan + ' mpileup2cns --output-vcf 1 --min-var-freq 0.35 --variants 1 > ' + chrom[0] + '.vcf'
        cmd_list.append(cmd)
    rf.close()
    proc = int(th) - 2
    job_manager(cmd_list, str(proc))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Germline variant caller.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample prefix')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Reference drive mount location.  Example would be /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, ref_mnt) = (inputs.config_file, inputs.sample, inputs.ref_mnt)
    varscan_germline(config_file, sample, inputs.ref_mnt)
