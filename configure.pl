#!/usr/bin/perl
use strict;
# use ExtUtils::Installed;
use File::Spec;
use File::Basename;


# my @req_modules = qw(IPC::System::Simple);
#
# @req_modules = sort @req_modules;
my $i = 1;
print "==== sTRAM checklist ====\n\n";

# print $i++ . ". Checking for required modules...\n";
# my $installed_modules = ExtUtils::Installed->new;

# foreach my $mod (@req_modules) {
# 	my $res = $installed_modules->validate($mod);
# }

print $i++ .". Checking for required software...\n";
my @req_software = qw(velveth velvetg blastn makeblastdb gawk);

my $sw_ready = 1;
foreach my $sw (@req_software) {
	my $result = system_call("$sw 2>&1 1>/dev/null");

	if ($result == 127) {
		print "   ...$sw is not installed.\n";
		$sw_ready = 0;
	} else {
		print "   ...$sw is present.\n";
	}
}

if ($sw_ready == 0) {
	print "You need to install some software before you can run sTRAM. Continue checks? [Y/n]\n";
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
