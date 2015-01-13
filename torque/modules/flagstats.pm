use strict;
use warnings;

package flagstats;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 8;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my $task = __PACKAGE__;
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=1\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "$main::SAMTOOLS flagstat $main::sample.srt.bam > $main::sample.srt.flagstats\n";
  print SCRIPT "$main::SAMTOOLS flagstat $main::sample.rmdup.srt.bam > $main::sample.rmdup.srt.flagstats\n";

  close SCRIPT;
  return $script;
}

1;
