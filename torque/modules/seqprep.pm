use strict;
use warnings;

package seqprep;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 8;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my $f1 = $main::SAMPLES{$main::sample}{f1};
  my $f2 = $main::SAMPLES{$main::sample}{f2};

  my $task = __PACKAGE__;
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=1\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "$main::SEQPREP -6 -f $f1 -r $f2 -1 $f1.clip.fq.gz -2 $f2.clip.fq.gz -s $main::sample.me.fq.gz > $main::LOGS_DIR/$main::sample.seqprep.log 2>&1\n";

  close SCRIPT;
  return $script;
}

1;
