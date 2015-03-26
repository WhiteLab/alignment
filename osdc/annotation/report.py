#!/usr/bin/env python
import argparse
import re
import sys

class Reporter:

  def __init__(self, infile, filt=False):
    self.infile = infile
    self.filter_missense_nonsense_only = filt

    self.__build_regex()
    self.__identify_columns()
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
      #'?(?P<warnings_errors>[^\|]*)?'
    ]))


  def __identify_columns(self):
    self.columns = ['chr', 'pos', 'context', 'ref', 'alt', 'normal_ref_count',
                    'normal_alt_count', 'tumor_ref_count', 'tumor_alt_count',
                    'dbSnp_id','gene', 'effect', 'coding', 'codon_change',
                    'amino_acid_change', 'amino_acid_length', 'mutect_call']


  def parse_infile(self):
    self.outstring = '\t'.join(self.columns) + '\n'

    column_refs = '' # lookup for irregularly placed columns

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
        report.append(line[column_refs.index('contig')]) # chromosome
        report.append(line[column_refs.index('0')]) # position
        report.append(line[column_refs.index('context')]) # context
        report.append(line[column_refs.index('REF_ALLELE')]) # reference
        report.append(line[column_refs.index('ALT_ALLELE')]) # alternative
        report.append(line[column_refs.index('n_ref_count')]) # normal reference count
        report.append(line[column_refs.index('n_alt_count')]) # normal alternative count
        report.append(line[column_refs.index('t_ref_count')]) # tumor reference count
        report.append(line[column_refs.index('t_alt_count')]) # tumor alternative count
        # parse context for dbSnp id
        id_check=re.search('rs(\d+)',line[column_refs.index('context')])
        if id_check:
          report.append(id_check.group(1))
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
          self.outstring += '\t%s' % line[-1] # keep?
          self.outstring += '\n'

          break

    except IOError:
      print >>sys.stderr, "ERROR: Unable to open -i/--infile: " + self.infile
      exit(1)


  def print_output(self):
    print >>sys.stdout, self.outstring
    return 0


def main():
  parser = argparse.ArgumentParser(
    description='parse snpEff annotated output into a digestable report.')
  parser.add_argument('-i', '--infile', action='store', dest='infile',
                      help='snpEff annotated variant file')
  parser.add_argument('-f', '--filter_missense_nonsense_only',
                      action='store_true', dest='f', default=False,
                      help='Apply a filter that only reports NONSENSE and ' \
                          + 'MISSENSE vars')
  if len(sys.argv) == 1:
    parser.print_help()
    exit(1)
  args = parser.parse_args()

  r = Reporter(args.infile, args.f)
  r.print_output()

if __name__ == '__main__':
  main()

