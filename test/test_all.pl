#!/usr/bin/env perl
use strict;
use File::Path qw (make_path);
use FindBin;
use File::Temp qw(tempdir);
use lib "$FindBin::Bin/../lib";
use System;
use constant ECHO => 1;

my $i = 0;
my $result = 0;
my $debug_flag = "";
my $log_flag = "";
my $executing_path = "$FindBin::Bin";
my $temp_dir = tempdir(CLEANUP => 1);
if (@ARGV[0] =~ /debug/) {
	$temp_dir = $ARGV[0];
	make_path ("debug");
	set_debug(1);
	open FH, ">", "debug.log";
	close FH;
	set_log("debug.log");
	$debug_flag = "-debug";
	$log_flag = "-log " . System::get_log_file();
}

$temp_dir = File::Spec->rel2abs($temp_dir);

##########################################################################################
## Testing format_sra.pl
##########################################################################################
printlog (++$i .". Checking that format_sra works correctly...", ECHO);
$result = system_call ("perl $executing_path/../format_sra.pl -in $executing_path/test_sra.fasta -out $temp_dir/test_db -num 7 $debug_flag $log_flag", 1);
if ($result != 0) {
	print "\nFormat_sra failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = `tail -n +2 $temp_dir/test_db.atram > $temp_dir/test_db.test`;
$result = `diff $executing_path/test_atram.txt $temp_dir/test_db.test > $executing_path/test.results.$i.diff`;
if ($result == 1) {
	print "\nFormat_sra returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

printlog (++$i .". Checking a defective file...", ECHO);

$result = system_call ("perl $executing_path/../format_sra.pl -in $executing_path/test_bad.fasta -out $temp_dir/test_db $debug_flag $log_flag", 1);
if ($result == 0) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

##########################################################################################
## Testing aTRAM.pl
##########################################################################################

printlog (++$i .". Checking that aTRAM works correctly...", ECHO);
$result = system_call ("perl $executing_path/../aTRAM.pl -db $temp_dir/test_db -target $executing_path/testref.fasta -out $temp_dir/test_atram $debug_flag $log_flag", 1);
if ($result != 0) {
	print "\aTRAM failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = `diff $executing_path/test_results_atram.txt $temp_dir/test_atram.results.txt > $executing_path/test.results.$i.diff`;
if ($result == 1) {
	print "\aTRAM returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

##########################################################################################
## Testing AlignmentPipeline.pl
##########################################################################################

printlog (++$i .". Checking that AlignmentPipeline works correctly...", ECHO);

# make some test target files and sample files:
open FH, ">", "$temp_dir/test.samples";
print FH "test\t$temp_dir/test_db.atram";
close FH;

open FH, ">", "$temp_dir/test.targets";
print FH "protein\t$executing_path/protref.fasta\n";
print FH "region\t$executing_path/testref.fasta\n";
print FH "bad\t$executing_path/badref.fasta\n";
print FH "complete\t$executing_path/completeref.fasta\n";
close FH;

$result = system_call ("perl $executing_path/../Pipelines/AlignmentPipeline.pl -samples $temp_dir/test.samples -targets $temp_dir/test.targets -out $temp_dir/test_ap -iter 5 $debug_flag $log_flag");
if ($result != 0) {
	print "\nAlignmentPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

if (-e "$temp_dir/test_ap/results.txt") {
	$result = `diff $executing_path/test_results_ap.txt $temp_dir/test_ap/results.txt > $executing_path/test.results.$i.diff`;
	if ($result == 1) {
		print "\nAlignmentPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
		exit;
	}
} else {
	print "\nAlignmentPipeline did not execute. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

##########################################################################################
## Testing BasicPipeline.pl
##########################################################################################

printlog (++$i .". Checking that BasicPipeline works correctly...", ECHO);
$result = system_call ("perl $executing_path/../Pipelines/BasicPipeline.pl -samples $temp_dir/test.samples -targets $temp_dir/test.targets -out $temp_dir/test_bp -iter 5 $debug_flag $log_flag");
if ($result != 0) {
	print "\nBasicPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = `grep -h '>' $temp_dir/test_bp/test/*.best.fasta > $temp_dir/test_bp/results.txt`;
if ((-s "$temp_dir/test_bp/results.txt") > 0) {
	$result = `diff $executing_path/test_results_bp.txt $temp_dir/test_bp/results.txt > $executing_path/test.results.$i.diff`;
	if ($result == 1) {
		print "\nBasicPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
		exit;
	}
} else {
	print "\nBasicPipeline did not execute. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";


print "\nAll tests successfully passed.\n\n";
system_call("rm -r $executing_path/test.results.*");
