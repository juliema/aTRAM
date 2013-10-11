#!/usr/bin/perl
use strict;
use ExtUtils::Installed;
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
my $temp = `which velveth`;
if ($temp == 0) {
	print "   Velveth is not installed in the usual place. If it is installed and accessible with a specific path, please enter it: ";
	my $userpath = <STDIN>;
	if ($userpath =~ /[qQeE]\n/) {
		print "quit";
	} else {
	print "  you said $userpath";
	}
} else {
	print "velveth is present\n";
}


# print $i++ .". Checking that prepare_files works correctly...\n";
# system ("cp test_good.fastq test.fastq");
#
# eval {
# 	system ("perl ../0-prepare_files.pl test.fastq");
# };
# if ($@) {
# 	print "Test failed. Please contact the developers with details of this failure at https://github.com/juliema/TRAM/issues.\n";
# 	exit;
# }
#
# system ("cp test_bad.fasta test.fasta");
# eval {
# 	my $result = system ("perl ../0-prepare_files.pl test.fasta");
# };
#
# if ($@) {
# 	print "...0-prepare_files passes tests.\n";
# }
