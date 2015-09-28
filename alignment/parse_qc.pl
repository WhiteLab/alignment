#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 2){
    print STDERR "Usage: $0 {output flag 0 - table, 1 - curl/jason command/output 2 - both}{sample prefix}\n";
    exit(1);
}
my $flag=$ARGV[0];
my $prefix=$ARGV[1];

my %files = archiveLookup();
my $tbl.="BionimbusID\tDate\tMachine\tRun\tBarCode\tLane\tread_length\ttotal_reads\tpost_align_reads\tfraction_aligned\tpost_rmdup_reads\tfraction_rmduped\ttarget_size\taligned_bp_ot\tfraction_aligned_bp_ot\tfraction_sequenced_bp_ot\taverage_ot_cov\tcoverage1\tcoverage2\tcoverage8\tfraction_coverage1\tfraction_coverage2\tfraction_coverage8\taligned_target_enrichment\tsequenced_target_enrichment\tmedian_insert_size\tmedian_absolute_deviation\tmean_insert_size\tinsert_standard_devation\n";

my %destinations;
my %runs;
my $GENOME_SIZE = 3137161264;
my $hostname='localhost';
my $port=9200;

my @samples = sort keys %files;
foreach my $sample (@samples) {
    print STDERR "$sample\n";
    my @set = sort keys %{$files{$sample}};
    
}
my %data = directory_search(".");
output_stats(\%data,\$flag);
if($flag == 0 || $flag == 2){
    open(TBL, ">","$prefix.qc_stats.txt");
    print TBL $tbl;
    close TBL;
}
sub output_stats
{
    #date variables
    my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    (my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,my $isdst)=gmtime();
    my $cur_date="$days[$wday], $mday $months[$mon] ".($year+1900)." $hour:$min:$sec GMT";
    
    my $data_ref = $_[0];
    my $f_ref=$_[1];
    my $f=$$f_ref;
    my %data = %$data_ref;
    my @samples = sort (keys %data);
    foreach my $sample (@samples) {
	my @annot = split (/\_/, $sample);
	my ($db_id,$date,$machine,$run,$barcode,$lane);
	$db_id = $annot[0];
	if($data{$sample}{SRC} eq 'HGAC'){
	    $date = $annot[1];
	    $machine = $annot[2];
	    $run = $annot[3];
	    $barcode = $annot[4];
	    $lane = $annot[5];
	}
	else{
	    ($date,$machine,$run,$barcode)=("NA","NA","NA","NA");
	    $lane = $annot[-1];
	}
	unless ($data{$sample}{LEN}) {
	    $data{$sample}{LEN} = 100;
	}
	unless ($data{$sample}{READS}) {
	    $data{$sample}{READS}="NA";
	}
	if ($data{$sample}{RAW}{MAP}) {
	    # add DB ID, Run date, machine SN, run #, flow cell barcode, lane
	    $tbl.="$db_id\t$date\t$machine\t$run\t$barcode\t$lane\t";
	    
	    # aggregate arguments and structure for curl/JSON command output
	    my $json.="{\n\t\"Bionimbus_id\": \"".$db_id."\",\n\t\"Run\":";
	    $json.=" \"".$run."\",\n\t\"Date\": \"".$date."\",\n\t\"Machine\": \"".$machine."\",\n\t\"BarCode\": \"".$barcode."\",\n\t\"Lane\": \"".$lane."\",\n\t\"read_length\": \"".$data{$sample}{LEN}."\",\n\t\"total_reads\": \"".$data{$sample}{RAW}{TOTAL}."\",\n";
	    # get Total reads, Aligned reads, rmdup aligned reads, Fraction duplicate reads
	    my $ali_frac = sprintf ("%.4f", $data{$sample}{RAW}{MAP}/$data{$sample}{RAW}{TOTAL});
	    my $dup_frac = sprintf ("%.4f", 1 - ($data{$sample}{RMDUP}{MAP} / $data{$sample}{RAW}{TOTAL}));
	    
	    # get # on target bases, Read length, Fraction on target bases, On target Enrichment
	    my $total_bp = $data{$sample}{RAW}{TOTAL} * $data{$sample}{LEN};
	    my $aligned_rmdup_bp = $data{$sample}{RMDUP}{MAP} * $data{$sample}{LEN};
	    my $on_target_reads = $data{$sample}{COV}{TOTAL} / $data{$sample}{LEN};
	    my $on_target_bp = $data{$sample}{COV}{TOTAL};
	    my $frac_aligned_bp_ot = sprintf ("%.4f", $on_target_bp / $aligned_rmdup_bp);
	    my $frac_sequenced_bp_ot = sprintf ("%.4f", $on_target_bp / $total_bp);
	    # add alignment stats to JSON statement
	    # median_insert_size\tmedian_absolute_deviation\tmean_insert_size\tinsert_standard_devation
	    #add date to json object
# Wed, 04 Mar 2015 01:41:44 GMT                                                                                                                                               
	    my $cur_date=$days[$wday].", $mday ".$months[$mon]." ".($year+1900)." ".$hour.":".$min.":".$sec." GMT";
	    $json.="\t\"post_align_reads\": \"".$data{$sample}{RAW}{MAP}."\",\n\t\"fraction_aligned\": \"".$ali_frac."\",\n\t\"post_rmdup_reads\": \"".$data{$sample}{RMDUP}{MAP}."\",\n\t\"fraction_rmduped\": \"".$dup_frac."\",\n\t\"target_size\": \"".$data{$sample}{COV}{TARGET}."\",\n\t\"aligned_bp_ot\": \"".$data{$sample}{COV}{TOTAL}."\",\n\t\"frac_aligned_bp_ot\": \"".$frac_aligned_bp_ot."\",\n\t\"fraction_sequenced_bp_ot\": \"".$frac_sequenced_bp_ot."\",\n\t\"median_insert_size\": \"".$data{$sample}{MI}."\",\n\t\"median_absolute_deviation\": \"".$data{$sample}{MA}."\",\n\t\"mean_insert_size\": \"".sprintf("%.2f",$data{$sample}{XI})."\",\n\t\"insert_standard_deviation\": \"".sprintf("%.2f",$data{$sample}{SI})."\",\n\t\"date_aligned\": \"".$cur_date."\",\n";
	    
	    my $average_ot_coverage = sprintf ("%.2f", $on_target_bp / $data{$sample}{COV}{TARGET});
	    my $on_target_frac = sprintf ("%.4f", $on_target_reads / $data{$sample}{RMDUP}{MAP} );
	    
	    # calculate enrichment
	    
	    my $frac_sizes = sprintf ("%.4f", $GENOME_SIZE / $data{$sample}{COV}{TARGET});
	    my $frac_aligned_reads = sprintf ("%.4f", $aligned_rmdup_bp / $on_target_bp);
	    my $aligned_target_enrichment = sprintf ("%.2f", $frac_sizes / $frac_aligned_reads);
	    
	    my $frac_sequenced_reads = sprintf ("%.4f", $total_bp / $on_target_bp);
	    my $sequenced_target_enrichment = sprintf ("%.2f", $frac_sizes / $frac_sequenced_reads);
	    
	    #HEADER: read_length total_reads post_align_reads fraction_aligned post_rmdup_reads fraction_rmduped target_size aligned_bp_ot fraction_aligned_bp_ot fraction_sequenced_bp_ot average_ot_cov coverage1 coverage2 coverage8 fraction_coverage1 fraction_coverage2 fraction_coverage8 aligned_target_enrichment sequenced_target_enrichment\n";
	    $tbl.="$data{$sample}{LEN}\t";			# read_length
	    $tbl.="$data{$sample}{RAW}{TOTAL}\t";	# total_reads
	    $tbl.="$data{$sample}{RAW}{MAP}\t";		# post_align_reads
	    $tbl.="$ali_frac\t";					# fraction_aligned
	    $tbl.="$data{$sample}{RMDUP}{MAP}\t";	# post_rmdup_reads
	    $tbl.="$dup_frac\t";					# fraction_rmduped
	    $tbl.="$data{$sample}{COV}{TARGET}\t";	# target_size
	    $tbl.="$data{$sample}{COV}{TOTAL}\t";	# aligned_bp_ot
	    $tbl.="$frac_aligned_bp_ot\t";			# fraction_aligned_bp_ot
	    $tbl.="$frac_sequenced_bp_ot\t";		# fraction_sequenced_bp_ot
	    $tbl.="$average_ot_coverage\t";			# average_ot_cov
				# add coverage stats to json
	    $json.="\t\"average_ot_coverage\": \"".$average_ot_coverage."\",\n\t";
	    for('1','2','8') {            			# cov 1x 2x 8x
		$tbl.= $data{$sample}{COV}{"RS$_"} . "\t";
		my $key="RS$_";
		$json.="\"coverage$_\": \"".$data{$sample}{COV}{$key}."\",\n\t";
	    }
	    for('1','2','8') {						# frac cov 1x 2x 8x
		my $fc = sprintf ("%.4f", $data{$sample}{COV}{"RS$_"} / $data{$sample}{COV}{TARGET});
		$tbl.= $fc."\t";
		$json.="\"fraction_coverage$_\": \"".$fc."\",\n\t";
	    }
	    $tbl.= "$aligned_target_enrichment\t";
	    $tbl.= "$sequenced_target_enrichment\t";
	    $tbl.="$data{$sample}{MI}\t$data{$sample}{MA}\t".sprintf("%.2f",$data{$sample}{XI})."\t".sprintf("%.2f",$data{$sample}{SI})."\n";
	    $json.="\"aligned_target_enrichment\": \"".$aligned_target_enrichment."\",\n\t\"sequenced_target_enrichment\": \"".$sequenced_target_enrichment."\"\n}";
	    #	    system($json);
	    if($f == 1 || $f == 2){
		open(CJ, ">>", "$prefix.qc_stats.json");
		print CJ $json;
		close CJ;
	    }
	}
    }
}
sub directory_search
{
    my %data;
    my @dir = `ls $_[0]`;
    foreach my $file (@dir) {
	chomp $file;
	if ($file =~ /(.+[1-8])\.rmdup.srt.bam.flagstats$/ || $file =~ /(.+[1-8]).rmdup.accepted_hits.bam.flagstats$/) {
	    my $sample = $1;
	    %{$data{$sample}{RMDUP}} = parseFS("$_[0]/$file");
	    my @name = split (/_/, $sample);
	    $data{$name[0]}{RUNS}{$_[0]}{$name[-1]}++;
	    
	}
	elsif ($file =~ /(.+[1-8])\.srt\.bam\.flagstats$/ || $file =~ /(.+[1-8]).accepted_hits.bam.flagstats$/) {
	    my $sample = $1;
	    %{$data{$sample}{RAW}} = parseFS("$_[0]/$file");
	}
	elsif ($file =~ /(.+[1-8]).(genome|exome|capture)\.hist$/ || $file =~ /(.+[1-8]).rnaSeq.covWhist$/) {
	    my $sample = $1;
	    ($data{$sample}{COV}{TOTAL}, $data{$sample}{COV}{RS1}, $data{$sample}{COV}{RS2}, $data{$sample}{COV}{RS8}, $data{$sample}{COV}{TARGET}) = parseCoverage("$_[0]/$file");
	}
	elsif ($file =~ /(.+[1-8])_[12]_sequence.txt.gz$/ || $file =~ /(.+[1-8])_sequence.txt.gz$/ || $file=~ /(^\S+)_\D*\d\.f\w*q\.gz$/) {
	    my $sample = $1;
	    my $temp = `gzip -dc $_[0]/$file | head -n 2 | tail -n 1`;
	    chomp $temp;
	    my $length = length($temp);
	    $data{$sample}{LEN} = $length;
	    if($file=~/_sequence.txt.gz$/){
		$data{$sample}{SRC}='HGAC';
	    }
	    else{
		$data{$sample}{SRC}='other';
	    }
	}
	elsif ($file=~/(.+[1-8])\.insert_metrics.hist$/){
	    my $sample=$1;
	    ($data{$sample}{MI},$data{$sample}{MA},$data{$sample}{XI},$data{$sample}{SI})=parseIns($file);
	}
	
    }
    return(%data);
}
sub parseFS
{
    my %data;
    my ($fs)  = @_;
    open IN, $fs or die "Can't open $fs!\n";
    my $line = <IN>;
    chomp $line;
    if($line =~ /^(\d+)\s/) {
	$data{TOTAL} = $1;
    }
    $line = <IN>;
    $line = <IN>;
    chomp $line;
    if($line =~ /^(\d+).*\((\d+\.\d+)/) {
	$data{MAP} = $1;
	$data{PCT} = $2;
    }
    return(%data);
}
sub parseCoverage
{
    my $total_bases = 0;
    my $rs1_bases = my $rs2_bases = my $rs8_bases = 0;
    my $target_size = 0;
    my ($fs)  = @_;
    open IN, $fs or die "Can't open $fs!\n";
    my $target_size_flag = 0;
    while (<IN>) {
	chomp $_;
	my @temp = split (/\s+/, $_);
	if(!$target_size_flag){
	    $target_size = $temp[3];
	    $target_size_flag = 1;
	}
	$total_bases = $total_bases + $temp[1]*$temp[2];
	$target_size = $temp[3];
	if ($temp[1] >= 1) {
	    $rs1_bases = $rs1_bases + $temp[2];
	}
	if ($temp[1] >= 2) {
	    $rs2_bases = $rs2_bases + $temp[2];
	}
	if ($temp[1] >= 8) {
	    $rs8_bases = $rs8_bases + $temp[2];
	}
    }
    return($total_bases, $rs1_bases, $rs2_bases, $rs8_bases, $target_size);
}

sub parseIns
{
    my ($f)=@_;
    open(INS, $f) or die "Can't open insert histogram file $f\n";
    while(my $line=<INS>){
	chomp $line;
	if($line=~/^MEDIAN/){
	    my $stats=<INS>;
	    close INS;
	    my @stats= split /\t/,$stats;
	    return($stats[0],$stats[1],$stats[4],$stats[5]);
	}
    }
}
sub annotIn
{
    my %annotation;
    open IN, $_[0] or die "Can't open $_[0]!\n";
    while (<IN>) {
	chomp $_;
	my @temp = split (/\t/, $_);
	my $db_id = $temp[0];
	$annotation{$db_id}{EXP_TYPE} = $temp[3];
	$annotation{$db_id}{SAMPLE_NAME} = $temp[1];
	$annotation{$db_id}{USER} = $temp[6];
	$annotation{$db_id}{PROJECT} = $temp[2];
    }
    return(%annotation);
}
sub archiveLookup
{
    my %files;
    my @dir = `ls .`;
    foreach my $file (@dir) {
	chomp $file;
	# add flexibility for fastq naming conventions
	if ($file =~ /(.+)_sequence.txt(.gz)*/) {
	    my @temp = split (/_/, $file);
	    my $db_id = $temp[0];
	    my $date = $temp[1];
	    my $mid = $temp[2];
	    my ($year, $month, $day) = (0,0,0);
	    if ($date =~ /(\d\d)(\d\d)(\d\d)/) {
		$year = $1;
		$month = $2;
		$day = $3;
	    }
	    $files{$db_id}{$file}++;
	}
	if($file =~/(^\S+)_\D*\d\.f\w*q\.gz$/){
	    my @temp = split (/_/, $file);
	    my $db_id = $temp[0];
	    my $date = $temp[1];
	    my $mid = $temp[2];
	    my ($year, $month, $day) = (0,0,0);
	    if ($date =~ /(\d\d)(\d\d)(\d\d)/) {
		$year = $1;
		$month = $2;
		$day = $3;
	    }
	    $files{$db_id}{$file}++;
	    
	}
    }
    return(%files);
}
