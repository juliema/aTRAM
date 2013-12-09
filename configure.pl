#!/usr/bin/perl
use strict;
use File::Basename;
use File::Spec;


my $i = 1;
print "==== sTRAM checklist ====\n\n";

print $i++ .". Checking for required software...\n";
my @req_software = qw(blastn makeblastdb gawk);
my @assembly_software = qw (velveth velvetg Trinity.pl);

open CONFIG_FH, ">", "config.txt";
print CONFIG_FH "# Enter the full path for the software binary below:\n";
my $sw_ready = 1;
foreach my $sw (@req_software) {
	my $result = system_call("$sw 2>&1 1>/dev/null");
	my $fullpath = which ($sw);

	if ($result == 127) {
		print "   ...$sw is not installed.\n";
		$sw_ready = 0;
	} else {
		print "   ...$sw is present.\n";
	}
	print CONFIG_FH "$sw=$fullpath\n";
}

if ($sw_ready == 0) {
	print "You need to install some software before you can run sTRAM. \n";
	print "If software is installed but not included in PATH, edit the appropriate line in config.txt.\n";
	print "\nContinue checks? [Y/n]\n";
	my $userpath = <STDIN>;
	if ($userpath =~ /[nN]\n/) {
		exit;
	}
}

print $i++ .". Checking for assembly software...\n";
foreach my $sw (@assembly_software) {
	my $result = system_call("$sw 2>&1 1>/dev/null");
	my $fullpath = which ($sw);

	if ($result == 127) {
		print "   ...$sw is not installed.\n";
		$sw_ready = 0;
	} else {
		print "   ...$sw is present.\n";
	}
	print CONFIG_FH "$sw=$fullpath\n";
}

close CONFIG_FH;
if ($sw_ready == 0) {
	print "I couldn't find some binaries for assembly software. At least one assembler should be installed.\n";
	print "If it is installed but not included in PATH, edit the appropriate line in config.txt.\n";
	print "\nContinue checks? [Y/n]\n";
	my $userpath = <STDIN>;
	if ($userpath =~ /[nN]\n/) {
		exit;
	}
}

print $i++ .". Checking that makelibrary works correctly...\n";

my $executing_path = dirname(__FILE__);

unless (system_call ("cp $executing_path/test/test_good.fastq $executing_path/test_inst.fastq") == 0) {
	die "Couldn't find test_good.fastq";
}

my $result = system_call ("perl $executing_path/makelibrary.pl $executing_path/test_inst.fastq");
if ($result == 1) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/TRAM/issues.\n";
	exit;
}

print $i++ .". Checking a defective file...\n";
unless (system_call ("cp $executing_path/test/test_bad.fasta $executing_path/test_inst.fasta") == 0) {
	die "Couldn't find test_bad.fasta";
}

my $result = system_call ("perl $executing_path/makelibrary.pl $executing_path/test_inst.fasta");
if ($result == 0) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/TRAM/issues.\n";
	exit;
}

print "Looks good! You are ready to sTRAM it up!\n";
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

sub which {
	my $cmd = shift;

	my @pathlist = File::Spec->path();
	foreach my $path (@pathlist) {
		my $cmdpath = File::Spec->catpath("", $path, $cmd);
		if (-x $cmdpath) {
			return "$cmdpath";
		}
	}
	return "";
}
