use strict;
use warnings;

package picard_rmdup_se;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 2;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my $task = __PACKAGE__;
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=4\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "$main::JAVA -Xmx8g -jar $main::PICARD MarkDuplicates CREATE_INDEX=true TMP_DIR=$main::PICARD_TMP REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 INPUT=$main::sample.me.srt.bam OUTPUT=$main::sample.me.rmdup.srt.bam METRICS_FILE=$main::sample.me.rmdup.srt.metrics VALIDATION_STRINGENCY=LENIENT > $main::LOGS_DIR/$main::sample.picard.rmdup.se.log 2>&1\n";

  close SCRIPT;
  return $script;
}

1;
