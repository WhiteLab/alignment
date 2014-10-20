#!/usr/bin/perl -w

# Alignment Pipeline Submission Controller
# Manages the creation, submission, and tracking of PBS scripts for the processing
# of fastq files within one or more directories.

$|=1;

use strict;
use warnings;
use File::Basename;
BEGIN { push @INC, dirname(__FILE__).'/modules' };

# Rate controller
use SubmitterBot;
my $BOT = new SubmitterBot;

# List generators in order of execution.
use flagstats;
use genome_coverage;

my @GENERATORS = (
  flagstats::new(),
  genome_coverage::new(),
);

# Parse command line arguments
if(scalar @ARGV < 3) { print "Usage: <job_name_prefix> <max_nodes> [<file>]+\n" and exit; }
our $JOB_NAME = shift;
our $MAX_NODES = shift;
our %SAMPLES = ();
foreach my $sample (@ARGV) { # e.g. 2011-1502_111228_SN673_0122_AC028YACXX_1_1_sequence.txt.gz
  if($sample =~ /^(\S+).sp.rmdup.srt.bam$/) {
    print "$sample\n";
    $SAMPLES{$1}{f1} = $1.".sp.rmdup.srt.bam";
  }
}

# Global namespace definitions
our $LOGS_DIR    = 'logs';
our $SCRIPTS_DIR = 'scripts';
our $SAMTOOLS    = "/glusterfs/users/mark/src/samtools/samtools";
our $BEDTOOLS    = "/glusterfs/users/mark/src/bedtools/bin/bedtools";
our $REFSEQ_BED  = "/glusterfs/users/mark/data/hg19/targets/nimblegen.merged.bed";

# Ensure directories exist or complain
if ( ! -d $LOGS_DIR ) { mkdir $LOGS_DIR or die $!; }
if ( ! -d $SCRIPTS_DIR ) { mkdir $SCRIPTS_DIR or die $!; }

# Call generators in order specified.
$BOT->job_name($JOB_NAME);
foreach(@GENERATORS) {
  our $gen = $_;
  our @cmds = ();
  foreach our $sample (keys %SAMPLES) {
    my $script = $gen->{'print'}();
    push @cmds, "qsub $script";
  }
  $BOT->batch($gen->{'task'});
  $BOT->max_jobs($MAX_NODES*$gen->{'mult'});
  $BOT->jobs(\@cmds);
  $BOT->submit_jobs();
}
