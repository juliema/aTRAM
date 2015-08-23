#!/usr/bin/env perl
use strict;
use File::Path qw (make_path);
use FindBin;
use File::Temp qw(tempdir);
use lib "$FindBin::Bin/../lib";
use System;
use constant ECHO => 1;
use Configuration;

Configuration::initialize();


my $i = 0;
my $result = 0;
my $debug_flag = "";
my $log_flag = "";
my $executing_path = "$FindBin::Bin";
my $temp_dir = tempdir(CLEANUP => 1);
if (@ARGV[0] =~ /debug/) {
	$temp_dir = "debug";
	make_path ("debug");
set_debug(1);
$debug_flag = "-debug";
}

$temp_dir = File::Spec->rel2abs($temp_dir);

open FH, ">", "$temp_dir/debug.log";
close FH;
set_log("$temp_dir/debug.log");
$log_flag = "-log " . get_log_file();

##########################################################################################
## Testing Configuration
##########################################################################################
printlog (++$i .". Checking that configuration is correct...", ECHO);
unless (get_bin("blastn") ne "") {
	fail_with_msg ("Configuration is not valid...did you run configure.pl?");
}
printlog ("OK", ECHO);

##########################################################################################
## Testing format_sra.pl
##########################################################################################
printlog (++$i .". Checking that format_sra works correctly on an interleaved file...", ECHO);
$result = run_command ("$executing_path/../format_sra.pl", "-in $executing_path/test_sra.fasta -out $temp_dir/test_db -num 7 $debug_flag $log_flag", {"no_exit"=>1});
if ($result != 0) {
	fail_with_msg ("Format_sra failed.");
}

`tail -n +2 $temp_dir/test_db.atram > $temp_dir/test_db.test`;
`diff $executing_path/test_format.txt $temp_dir/test_db.test > $executing_path/test.results.$i.diff`;
if ($? != 0) {
	fail_with_msg ("Format_sra returned incorrect results.");
}
`rm $executing_path/test.results.$i.diff`;

printlog ("OK", ECHO);

printlog (++$i .". Checking that format_sra works correctly on two separate files...", ECHO);
$result = run_command ("$executing_path/../format_sra.pl", "-1input $executing_path/test_sra1.fasta -2input $executing_path/test_sra2.fasta -out $temp_dir/test_db -num 7 $debug_flag $log_flag", {"no_exit"=>1});
if ($result != 0) {
	fail_with_msg ("Format_sra failed.");
}

`tail -n +2 $temp_dir/test_db.atram > $temp_dir/test_db.test`;
`diff $executing_path/test_format.txt $temp_dir/test_db.test > $executing_path/test.results.$i.diff`;
if ($? != 0) {
	fail_with_msg ("Format_sra returned incorrect results.");
}
`rm $executing_path/test.results.$i.diff`;

printlog ("OK", ECHO);

printlog (++$i .". Checking a defective file...", ECHO);

$result = run_command ("$executing_path/../format_sra.pl", "-in $executing_path/test_bad.fasta -out $temp_dir/test_db_bad $debug_flag $log_flag", {"no_exit"=>1});

if ($result == 0) {
	fail_with_msg ("Test failed.");
}

printlog ("OK", ECHO);

##########################################################################################
## Testing aTRAM.pl
##########################################################################################

printlog (++$i .". Checking that aTRAM works correctly...", ECHO);
$result = run_command ("$executing_path/../aTRAM.pl" ,"-db $temp_dir/test_db -target $executing_path/testref.fasta -out $temp_dir/test_atram $debug_flag $log_flag", {"no_exit"=>1});
if ($result != 0) {
	fail_with_msg ("aTRAM failed.");
}

$result = `grep -c "." $temp_dir/test_atram.results.txt`;
chomp $result;
if ($result != 10) {
	fail_with_msg ("aTRAM returned incorrect results.");
}

printlog ("OK", ECHO);

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

##########################################################################################
## Testing BasicPipeline.pl
##########################################################################################

printlog (++$i .". Checking that BasicPipeline works correctly...", ECHO);
$result = run_command ("$executing_path/../Pipelines/BasicPipeline.pl", "-samples $temp_dir/test.samples -targets $temp_dir/test.targets -out $temp_dir/test_bp -iter 5 $debug_flag $log_flag", {"no_exit"=>1});
if ($result != 0) {
	fail_with_msg ("BasicPipeline died in execution.");
}

$result = `grep -h '>' $temp_dir/test_bp/test/*.best.fasta > $temp_dir/test_bp/results.txt`;
if ((-s "$temp_dir/test_bp/results.txt") > 0) {
	$result = `grep -c "." $temp_dir/test_bp/results.txt`;
	chomp $result;
	if ($result != 12) {
		fail_with_msg ("BasicPipeline returned incorrect results.");
	}
} else {
	fail_with_msg ("BasicPipeline did not execute.");
}

printlog ("OK", ECHO);

##########################################################################################
## Testing AlignmentPipeline.pl
##########################################################################################

printlog (++$i .". Checking that AlignmentPipeline works correctly...", ECHO);

$result = run_command ("$executing_path/../Pipelines/AlignmentPipeline.pl", "-samples $temp_dir/test.samples -targets $temp_dir/test.targets -out $temp_dir/test_ap -iter 5 $debug_flag $log_flag", {"no_exit"=>1});
if ($result != 0) {
	fail_with_msg ("AlignmentPipeline died in execution.");
}

if (-e "$temp_dir/test_ap/results.txt") {
	$result = `grep -c "." $temp_dir/test_ap/results.txt`;
	chomp $result;
	if ($result != 3) {
		fail_with_msg ("AlignmentPipeline returned incorrect results.");
	}
} else {
	fail_with_msg ("AlignmentPipeline did not execute.");
}

printlog ("OK", ECHO);

printlog ("\nAll tests successfully passed.\n", ECHO);
exit 0;

sub fail_with_msg {
	my $msg = shift;

	close_log();
	print "$msg\nLast logged output:\n";

	my $result = `tail $temp_dir/debug.log`;
	$result =~ s/.*: $i\..*?$//ms;
	print $result;

	exit 1;
}
