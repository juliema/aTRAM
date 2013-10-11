#!/usr/bin/perl
use strict;
use ExtUtils::Installed;
use File::Spec;
# use IPC::System::Simple qw(run system capture EXIT_ANY);


my @req_modules = qw(IPC::System::Simple);

@req_modules = sort @req_modules;
my $i = 1;

print $i++ . ". Checking for required modules...\n";
my $installed_modules = ExtUtils::Installed->new;
# my @modules = $installed_modules->modules();

foreach my $mod (@req_modules) {
	my $res = $installed_modules->validate($mod);
}

print $i++ .". Checking for required software...\n";
my @req_software = qw(velveth velvetg blastn makeblastdb);

my $sw_ready = 1;
foreach my $sw (@req_software) {
	open my $saveout, ">&STDOUT";
	open STDOUT, '>', File::Spec->devnull();

	my $result = system("$sw 2>&1 1>/dev/null");
	#Restore STDOUT
	open STDOUT, ">&", $saveout;

	$result = ($result >> 8);
	if ($result == 127) {
		print "   ...package $sw is not installed.\n";
		$sw_ready = 0;
	# 	my $userpath = <STDIN>;
	# 	if ($userpath =~ /[qQeE]\n/) {
	# 		print "quit";
	# 	} else {
	# 	print "  you said $userpath";
	# 	}
	} else {
		print "   ...package $sw is present.\n";
	}
}

if ($sw_ready == 0) {
	die "You need to install some software before you can run sTRAM.\n";
}

print $i++ .". Checking that prepare_files works correctly...\n";

open my $saveout, ">&STDOUT";
open STDOUT, '>', File::Spec->devnull();
system ("cp test_good.fastq test.fastq");
eval {
	system ("perl ../0-prepare_files.pl test.fastq");
};
open STDOUT, ">&", $saveout;


if ($@) {
	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/TRAM/issues.\n";
	exit;
}
# print $i++ .". Checking a defective file...\n";
#
# open my $saveout, ">&STDOUT";
# open my $saveerr, ">&STDERR";
# open STDERR, '>', File::Spec->devnull();
# system ("cp test_bad.fasta test.fasta");
# eval {
# 	system ("perl ../0-prepare_files.pl test_bad.fasta");
# };
# open STDOUT, ">&", $saveout;
# open STDERR, ">&", $saveerr;
#
# if ($@) {
# 	print "...0-prepare_files passes tests.\n";
# }
#
#
print "Looks good! You are ready to sTRAM it up!\n";
