use strict;
use warnings;

package fastx;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 2;
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
  print SCRIPT "#PBS -l nodes=1:ppn=4\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "gzip -dc $f1 | $main::FASTX -N -o $f1.qs &\n";
  print SCRIPT "gzip -dc $f2 | $main::FASTX -N -o $f2.qs &\n";
  print SCRIPT "wait\n";

  close SCRIPT;
  return $script;
}

1;
