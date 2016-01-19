#!/usr/bin/env python
import argparse
import re
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time


class Reporter:
    def __init__(self, infile, c, filt=False):
        self.infile = infile
        self.filter_missense_nonsense_only = filt
        self.c = c
        self.__build_regex()
        self.__identify_columns()
        if c != 'n':
            self.__index_intervals()
        self.outstring = '\t'.join(self.columns) + '\n'
        self.parse_infile()

    def __build_regex(self):
        self.regex = re.compile(r'(?P<eff>[\w\_]+)\(%s\)' % '\|'.join([
            '(?P<effect_impact>[^\|]*)',
            '(?P<functional_class>[^\|]*)',
            '(?P<codon_change>[^\|]*)',
            '(?P<amino_acid_change>[^\|]*)',
            '(?P<amino_acid_length>[^\|]*)',
            '(?P<gene_name>[^\|]*)',
            '(?P<transcript_biotype>[^\|]*)',
            '(?P<gene_coding>[^\|]*)',
            '(?P<transcript_id>[^\|]*)',
            '(?P<exon_intro_rank>[^\|]*)',
            '(?P<genotype_number>[^\|]*)',
            '?(?P<warnings_errors>[^\|]*)?'
        ]))

    def __identify_columns(self):
        self.columns = ['chr', 'pos', 'context', 'ref', 'alt', 'normal_ref_count',
                        'normal_alt_count', '%_normal_alt', 'tumor_ref_count', 'tumor_alt_count',
                        '%_tumor_alt', 'T/N_%_alt_ratio', 'dbSnp_id', 'gene', 'effect', 'coding', 'codon_change',
                        'amino_acid_change', 'amino_acid_length', 'mutect_call']
        if self.c != 'n':
            self.columns.append('on/off-target')
            sys.stderr.write(date_time() + 'Interval list given - marking on/off target hits\n')

    def __index_intervals(self):
        self.index = {}
        fh = open(self.c, 'r')
        for line in fh:
            m = re.search('(chr\w+):(\d+)-(\d+)', line)
            try:
                (chrom, start, end) = (m.group(1), m.group(2), m.group(3))
            except:
                sys.stderr.write(line + ' doesn\'t fit format (chr\w+):(\d+)-(\d+), skipping\n')
                continue
            if chrom not in self.index:
                self.index[chrom] = {}
            self.index[chrom][(int(start))] = int(end)

    def search_index(self, chrom, pos):
        f = 0
        if chrom in self.index:
            for start in sorted(self.index[chrom]):
                if start <= int(pos) <= self.index[chrom][start]:
                    f = 1
                    break
                elif start > int(pos):
                    break
        status = "OFF"
        if f == 1:
            status = "ON"
        return status

    def calc_pct(self, a, b):
        # return both formatted and unformatted
        ratio = float(b) / (float(a) + float(b)) * 100
        fmt = "{0:.2f}%".format(ratio)
        return ratio, fmt

    def parse_infile(self):

        column_refs = ''  # lookup for irregularly placed columns
        sys.stderr.write(date_time() + 'Processing ' + self.infile + '\n')
        try:
            for line in open(self.infile):
                # Skip headers.
                if line[0] == '#':
                    continue
                line = line.strip().split('\t')
                if line[0] == 'contig':
                    column_refs = line
                    continue

                report = list()
                report.append(line[column_refs.index('contig')])  # chromosome
                report.append(line[column_refs.index('0')])  # position
                # if context contains rs id, remove since it will be in another column
                id_check = re.search('(\w+);rs(\d+)', line[column_refs.index('context')])
                if id_check:
                    report.append(id_check.group(1))
                else:
                    report.append(line[column_refs.index('context')])  # context
                report.append(line[column_refs.index('REF_ALLELE')])  # reference
                report.append(line[column_refs.index('ALT_ALLELE')])  # alternative
                report.append(line[column_refs.index('n_ref_count')])  # normal reference count
                report.append(line[column_refs.index('n_alt_count')])  # normal alternative count
                # calculate % normal and tumor coverage, keep both formatted and unforamtted version in order to calculate that ratio of the ratios
                (n_rat, n_fmt) = self.calc_pct(line[column_refs.index('n_ref_count')],
                                               line[column_refs.index('n_alt_count')])
                report.append(n_fmt)
                report.append(line[column_refs.index('t_ref_count')])  # tumor reference count
                report.append(line[column_refs.index('t_alt_count')])  # tumor alternative count
                (t_rat, t_fmt) = self.calc_pct(line[column_refs.index('t_ref_count')],
                                               line[column_refs.index('t_alt_count')])
                report.append(t_fmt)
                t_n_rat = float(line[column_refs.index('t_alt_count')])
                if n_rat > 0:
                    t_n_rat = t_rat / n_rat
                report.append("{0:.2f}".format(t_n_rat))
                # parse context for dbSnp id
                if id_check:
                    report.append(id_check.group(2))
                else:
                    report.append('NA')
                # Interpret the effects of the variant.
                for match in self.regex.finditer(line[7]):
                    # NOTE this is a hack for Kevin's request
                    if self.filter_missense_nonsense_only:
                        if match.group('functional_class') not in ['MISSENSE', 'NONSENSE']:
                            continue

                    self.outstring += '\t'.join(report)
                    self.outstring += '\t%s' % match.group('gene_name')
                    self.outstring += '\t%s' % match.group('eff')
                    self.outstring += '\t%s' % match.group('gene_coding')
                    self.outstring += '\t%s' % match.group('codon_change')
                    self.outstring += '\t%s' % match.group('amino_acid_change')
                    self.outstring += '\t%s' % match.group('amino_acid_length')
                    self.outstring += '\t%s' % line[-1]  # keep?
                    if self.c != 'n':
                        chrom = line[column_refs.index('contig')]
                        pos = line[column_refs.index('0')]
                        flag = self.search_index(chrom, pos)
                        self.outstring += '\t%s' % flag
                    self.outstring += '\n'
                    break

        except IOError:
            print >> sys.stderr, "ERROR: Unable to open -i/--infile: " + self.infile
            exit(1)

    def print_output(self):
        sys.stdout.write(self.outstring)
        sys.stderr.write(date_time() + 'Changed last Report complete\n')
        return 0


def main():
    parser = argparse.ArgumentParser(
        description='parse snpEff annotated output into a digestable report.')
    parser.add_argument('-i', '--infile', action='store', dest='infile',
                        help='snpEff annotated variant file')
    parser.add_argument('-c', '--custom', action='store', dest='c',
                        help='bed file to mark whether hit was on or off-target. if not desired, enter \'n\' ')

    parser.add_argument('-f', '--filter_missense_nonsense_only',
                        action='store_true', dest='f', default=False,
                        help='Apply a filter that only reports NONSENSE and ' \
                             + 'MISSENSE vars')
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    args = parser.parse_args()

    r = Reporter(args.infile, args.c, args.f)
    r.print_output()


if __name__ == '__main__':
    main()
