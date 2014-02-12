#!/usr/bin/env perl
use strict;
# use File::Spec;
use FindBin;
use File::Temp qw(tempdir);

my $i = 0;
my $result = 0;
my $executing_path = "$FindBin::Bin";
my $temp_dir = tempdir(CLEANUP => 1);


##########################################################################################
## Testing format_sra.pl
##########################################################################################

print ++$i .". Checking that format_sra works correctly...";
unless (system_call ("cp $executing_path/test_sra.fasta $temp_dir/test_inst.fasta") == 0) {
	die "Couldn't find test_sra.fasta";
}

$result = system_call ("perl $executing_path/../format_sra.pl -in $temp_dir/test_inst.fasta -out $temp_dir/test_inst -num 7");
if ($result == 1) {
	print "\nFormat_sra failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = system_call ("tail -n +2 $temp_dir/test_inst.atram > $temp_dir/test_inst.test");
$result = system_call ("diff $executing_path/test_atram.txt $temp_dir/test_inst.test > $executing_path/test_inst.results.$i.diff");
if ($result == 1) {
	print "\nFormat_sra returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

##########################################################################################
## Testing AlignmentPipeline.pl
##########################################################################################

print ++$i . ". Checking that AlignmentPipeline works correctly...";

# make some test target files and sample files:
open FH, ">", "$temp_dir/test_inst.samples";
print FH "test\t$temp_dir/test_inst.atram";
close FH;

open FH, ">", "$temp_dir/test_inst.targets";
print FH "protein\t$executing_path/protref.fasta\n";
print FH "region\t$executing_path/testref.fasta\n";
print FH "bad\t$executing_path/badref.fasta\n";
print FH "complete\t$executing_path/completeref.fasta\n";
close FH;

$result = system_call ("perl $executing_path/../Pipelines/AlignmentPipeline.pl -samples $temp_dir/test_inst.samples -targets $temp_dir/test_inst.targets -out $temp_dir/test_inst");
if ($result == 1) {
	print "\nAlignmentPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = system_call ("diff $executing_path/test_results1.txt $temp_dir/test_inst/results.txt > $executing_path/test_inst.results.$i.diff");
if ($result == 1) {
	print "\nAlignmentPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "OK\n";

##########################################################################################
## Testing BasicPipeline.pl
##########################################################################################

print ++$i . ". Checking that BasicPipeline works correctly...";
$result = system_call ("perl $executing_path/../Pipelines/BasicPipeline.pl -samples $temp_dir/test_inst.samples -targets $temp_dir/test_inst.targets -out $temp_dir/test_inst");
if ($result == 1) {
	print "\nBasicPipeline died in execution. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

$result = system_call ("cat $temp_dir/test_inst/test/*.results.txt > $temp_dir/test_inst/results.txt");
$result = system_call ("diff $executing_path/test_results2.txt $temp_dir/test_inst/results.txt > $executing_path/test_inst.results.$i.diff");
if ($result == 1) {
	print "\nBasicPipeline returned incorrect results. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}
print "OK\n";


print "\nAll tests successfully passed.\n\n";
system_call("rm -r $executing_path/test_inst.results.*");

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
