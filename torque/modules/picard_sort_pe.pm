use strict;
use warnings;

package picard_sort_pe;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 1;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my $task = __PACKAGE__;
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=8\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "$main::JAVA -Xmx8g -jar $main::PICARD SortSam CREATE_INDEX=true TMP_DIR=$main::PICARD_TMP INPUT=$main::sample.bam OUTPUT=$main::sample.srt.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT > $main::LOGS_DIR/$main::sample.picard.sort.pe.log 2>&1\n";

  close SCRIPT;
  return $script;
}

1;
