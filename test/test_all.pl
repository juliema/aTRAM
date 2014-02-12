#!/usr/bin/env perl
use strict;
# use File::Spec;
use FindBin;

my $i = 1;
my $result = 0;
my $executing_path = "$FindBin::Bin";

select STDOUT;
print $i++ .". Checking that format_sra works correctly...";
unless (system_call ("cp $executing_path/test_sra.fasta $executing_path/test_inst.fasta") == 0) {
	die "Couldn't find test_sra.fasta";
}


$result = system_call ("perl $executing_path/../format_sra.pl -in $executing_path/test_inst.fasta -out $executing_path/test_inst -num 7");
if ($result == 1) {
	print "\nFormat_sra failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}
print "OK\n";

print $i++ . ". Checking that AlignmentPipeline works correctly...";
$result = system_call ("perl $executing_path/../Pipelines/AlignmentPipeline.pl -samples $executing_path/testsamples.txt -targets $executing_path/testtargets.txt -out $executing_path/test_inst");
if ($result == 1) {
	print "\nAlignmentPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = system_call ("diff $executing_path/test_results1.txt $executing_path/test_inst/results.txt > test_inst.results1.diff");
if ($result == 1) {
	print "\nAlignmentPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

print $i++ . ". Checking that BasicPipeline works correctly...";
$result = system_call ("perl $executing_path/../Pipelines/AlignmentPipeline.pl -samples $executing_path/testsamples.txt -targets $executing_path/testtargets.txt -out $executing_path/test_inst");
if ($result == 1) {
	print "\nBasicPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = system_call ("cat $executing_path/test_inst/test/*.results.txt > $executing_path/test_inst/results.txt");
$result = system_call ("diff $executing_path/test_results2.txt $executing_path/test_inst/results.txt > test_inst.results2.diff");
if ($result == 1) {
	print "\nBasicPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}
print "OK\n";

print "\nAll tests successfully passed.\n\n";
system_call("rm -r $executing_path/test_inst*");

sub system_call {
	my $cmd = shift;
	open my $saveout, ">&STDOUT";
	open my $saveerr, ">&STDERR";
	open STDOUT, '>', File::Spec->devnull();
	open STDERR, '>', File::Spec->devnull();

	my $exit_val = eval {
		system ($cmd);
	} >> 8;

	open STDOUT, ">&", $saveout;
	open STDERR, ">&", $saveerr;

	return $exit_val;
}
