use strict;
use warnings;

package picard_merge;

sub new {
  my $struct = {};
  $struct->{'task'} = __PACKAGE__;
  $struct->{'mult'} = 8;
  $struct->{'print'} = \&print;

  return bless $struct;
}

sub print {
  my $script = sprintf "$main::SCRIPTS_DIR/$task.%s.sh", scalar @main::cmds;
  open SCRIPT, '>', $script or die $!;

  print SCRIPT "#!/usr/bin/env bash\n";
  print SCRIPT "#PBS -N $main::JOB_NAME.$task\n";
  print SCRIPT "#PBS -l nodes=1:ppn=1\n";
  print SCRIPT "#PBS -j oe\n";
  print SCRIPT "#PBS -o $main::LOGS_DIR\n";
  print SCRIPT "cd \$PBS_O_WORKDIR\n";
  print SCRIPT "$main::JAVA -Xmx2g -jar $main::PICARD/MergeSamFiles.jar CREATE_INDEX=true MAX_RECORDS_IN_RAM=6000000 TMP_DIR=$main::PICARD_TMP  O=$main::sample.sp.rmdup.srt.bam  I=$main::sample.me.rmdup.srt.bam I=$main::sample.clip.rmdup.srt.bam SO=coordinate AS=true USE_THREADING=true VALIDATION_STRINGENCY=LENIENT > $main::LOGS_DIR/$main::sample.picard.merge.log 2>&1\n"; # NOTE Was using SILENT stringency

  close SCRIPT;
  return $script;
}

1;
