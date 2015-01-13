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
use fastx;
use bwa_mem_pe;
use picard_sort_pe;
use picard_rmdup_pe;
use flagstats;
use target_coverage;

my @GENERATORS = (
  fastx::new(),
  bwa_mem_pe::new(),
  picard_sort_pe::new(),
  picard_rmdup_pe::new(),
  flagstats::new(),
  target_coverage::new(),
);

# Parse command line arguments
if(scalar @ARGV < 3) { print "Usage: <job_name_prefix> <max_nodes> [<file>]+\n" and exit; }
our $JOB_NAME = shift;
our $MAX_NODES = shift;
our %SAMPLES = ();
foreach my $sample (@ARGV) { # e.g. 2011-1502_111228_SN673_0122_AC028YACXX_1_1_sequence.txt.gz
  if($sample =~ /^(\S+)_1_sequence\.txt\.gz$/) {
    print "$sample\n";
    $SAMPLES{$1}{f1} = $1."_1_sequence.txt.gz";
    $SAMPLES{$1}{f2} = $1."_2_sequence.txt.gz";
  }
}

# Global namespace definitions
our $LOGS_DIR    = 'logs';
our $SCRIPTS_DIR = 'scripts';
our $JAVA        = "/glusterfs/users/mark/src/jdk1.7.0_45/bin/java";
our $FASTX       = "/usr/local/tools/fastx_toolkit-0.0.13/src/fastx_quality_stats/fastx_quality_stats";
our $BWA         = "/glusterfs/users/mark/src/bwa-0.7.8/bwa";
our $SAMTOOLS    = "/glusterfs/users/mark/src/samtools/samtools";
our $BEDTOOLS    = "/glusterfs/users/mark/src/bedtools/bin/bedtools";
our $PICARD      = "/raid/users/jgrundst/TOOLS/picard/dist/picard.jar";
our $SEQPREP     = "/glusterfs/users/mark/src/seqprep/SeqPrep";
our $BWA_REF_FA  = "/glusterfs/users/mark/data/hg19/bwa/0.7.8/hg19.fa";
our $SAM_REF_FA  = "/glusterfs/users/mark/data/hg19/samtools/hg19.fa";
our $CAPTURE_REFSEQ_BED  = "/glusterfs/users/mark/data/hg19/targets/nimblegen.merged.bed";
our $PIGZ        = "/glusterfs/users/jgrundst/bin/pigz";
our $PICARD_TMP  = "PICARD_TMP";
our $MAXMEM      = 2000000000;

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
