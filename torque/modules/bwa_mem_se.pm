use strict;
use warnings;

package bwa_mem_se;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 1;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my @HGACID = split(/_/,$main::sample);
  my $RGRP = "\@RG\\tID:$main::sample\\tLB:$HGACID[0]\\tSM:$HGACID[0]\\tPL:illumina";

  my $task = __PACKAGE__;
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=8\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "($main::BWA mem -t \$PBS_NUM_PPN -R \"$RGRP\" -v 2 $main::BWA_REF_FA $main::sample.me.fq.gz | $main::SAMTOOLS view -bT $main::SAM_REF_FA - > $main::sample.me.bam) > $main::LOGS_DIR/$main::sample.bwa.se.log 2>&1\n";

  close SCRIPT;
  return $script;
}

1;
