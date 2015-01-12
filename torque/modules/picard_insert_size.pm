use strict;
use warnings;

package picard_insert_size;

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
  print SCRIPT "$main::JAVA -Xmx2g -jar $main::PICARD ",
  	"CollectInsertSizeMetrics I=$main::sample.rmdup.srt.bam ",
	"H=$main::sample.insert_metrics.hist ",
	"O=$main::sample.insert_metrics.txt > ",
	"$main::LOGS_DIR/$main::sample.picard.insert_size.log 2>&1\n"; # NOTE Was using SILENT stringency

  close SCRIPT;
  return $script;
}

1;
