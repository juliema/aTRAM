#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/lib";
use Configuration;

# find or make config.txt:
my $config_file = "$FindBin::Bin/lib/config.txt";
unless (-e $config_file) {
	# if the config_file isn't existing, make an empty one for filling.
	open CONFIG_FH, ">", $config_file;
	close CONFIG_FH;
}

Configuration::initialize();
open CONFIG_FH, ">", $config_file;
print CONFIG_FH "# Enter the full path for the software binary below:\n";

my $sw_ready = 1;
my $i = 1;
my $result = 0;
print "==== aTRAM checklist ====\n\n";

print $i++ .". Checking for required software...\n";
my @req_software = @Configuration::req_software;

foreach my $sw (@req_software) {
	my $fullpath = Configuration::find_bin($sw); # see if $sw has already been located.

	if ($fullpath eq "") {
		print "   ...$sw couldn't be found on this system.\n";
		$sw_ready = 0;
	} else {
		$result = system_call("$fullpath 2>&1 1>/dev/null");
		if ($result == 127) {
			print "   ...$sw was not found at $fullpath.\n";
			$sw_ready = 0;
		} else {
			print "   ...$sw is present.\n";
		}
	}

	print CONFIG_FH "$sw=$fullpath\n";
}

print $i++ .". Checking for assembly software...\n";
my $assemblers = Configuration::get_assemblers();
my $assembler_present = 0;
foreach my $assembler (keys %$assemblers) {
	my $assembler_ready = 0;
	print "   For assembler $assembler:\n";
	my $assembly_software = $assemblers->{$assembler};
	foreach my $sw (@$assembly_software) {
		my $fullpath = Configuration::find_bin($sw); # see if $sw has already been located.

		if ($fullpath eq "") {
			print "      ...$sw couldn't be found on this system.\n";
		} else {
			$result = system_call("$fullpath --version");
			if ($result == 127) {
				print "      ...$sw was not found at $fullpath.\n";
			} else {
				print "      ...$sw is present.\n";
				$assembler_ready++;
			}
		}

		print CONFIG_FH "$sw=$fullpath\n";
	}
	if ($assembler_ready == @$assembly_software) {
		print "   Assembler $assembler is ready to use.\n";
		$assembler_present++;
	} else {
		print "   Assembler $assembler cannot find all its binaries.\n";
	}
}
if ($assembler_present == 0) {
	$sw_ready == 0;
}

close CONFIG_FH;
if ($sw_ready == 0) {
	print "Software required by some parts of aTRAM were not found on this system.\n";
	print "If software is installed but not included in \$PATH, edit the appropriate line in config.txt.\n";
	print "\nContinue checks? [Y/n]\n";
	my $userpath = <STDIN>;
	if ($userpath =~ /[nN]\n/) {
		exit;
	}
}
my $executing_path = $FindBin::Bin;

print $i++ .". Checking that format_sra works correctly...\n";
unless (system_call ("cp $executing_path/test/test_good.fastq $executing_path/test_inst.fastq") == 0) {
	die "Couldn't find test_good.fastq";
}

$result = system_call ("perl $executing_path/format_sra.pl $executing_path/test_inst.fastq");
if ($result == 1) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print $i++ .". Checking a defective file...\n";
unless (system_call ("cp $executing_path/test/test_bad.fasta $executing_path/test_inst.fasta") == 0) {
	die "Couldn't find test_bad.fasta";
}

$result = system_call ("perl $executing_path/format_sra.pl $executing_path/test_inst.fasta");
if ($result == 0) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/aTRAM/issues.\n";
	exit;
}

print "Looks good! You are ready to aTRAM it up!\n";
system_call("rm $executing_path/test_inst.*");

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

