#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw(tempfile cleanup);
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
require Subfunctions;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $short_read_archive = "";
my $output_file = "";
my $help = 0;
my $debug = 0;
my $numlibraries = 0;
my $max_processes = 4;

GetOptions ('input=s' => \$short_read_archive,
            'output=s' => \$output_file,
            'numlibraries=i' => \$numlibraries,
            'debug' => \$debug,
			'max_processes|processes=i' => \$max_processes,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive) {
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
}

unless ($output_file) {
    $output_file = $short_read_archive;
}

my $tempdir = dirname ("$output_file");
my $log_file = "$output_file.log";
open my $log_fh, ">", $log_file or die "couldn't open $log_file\n";
set_log ($log_fh);

my $libsize = (-s $short_read_archive);
if ($numlibraries == 0) {
	$numlibraries = int($libsize / 5e8);
	if (($numlibraries % 2) == 0) {
		$numlibraries++;
	}
	printlog ("$short_read_archive is $libsize bytes; we will make $numlibraries libraries.");
} else {
	printlog ("making $numlibraries libraries.");
}

my @primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151);
my $multiplier = $primes[$libsize % (@primes)];

my @tempfiles = ();

my @out1_fhs = ();
my @out2_fhs = ();
my @out1_bucketfiles = ();
my @out2_bucketfiles = ();
my @out1_sortedfiles = ();
my @out2_sortedfiles = ();

# setting up tempfiles for sorting:
for (my $i=0; $i<$numlibraries; $i++) {
	my ($fh, $filename) = tempfile("$output_file.bucket.$i.1.XXXX", UNLINK => 1);
	push @out1_bucketfiles, $filename;
	push @out1_fhs, $fh;

	($fh, $filename) = tempfile("$output_file.bucket.$i.2.XXXX", UNLINK => 1);
	push @out2_bucketfiles, $filename;
	push @out2_fhs, $fh;

	($fh, $filename) = tempfile("$output_file.sorted.$i.1.XXXX", UNLINK => 1);
	push @out1_sortedfiles, $filename;

	($fh, $filename) = tempfile("$output_file.sorted.$i.2.XXXX", UNLINK => 1);
	push @out2_sortedfiles, $filename;
}

# Divide fasta/fastq short reads into buckets for sorting.
printlog ("Dividing fasta/fastq file into buckets for sorting.");
my $name = "";
my $seq = "";
my $seqlen = 0;
my $mapkey = 0;

open SEARCH_FH, "<", $short_read_archive;
while (my $line = readline SEARCH_FH) {
	chomp $line;
	if ($line =~ /^[@>](.*?)([\s\/])([12])/) {
		if ($name ne "") {
			if ($name =~ /\/1/) {
				print {$out1_fhs[$mapkey]} ">$name,$seq\n";
			} elsif ($name =~ /\/2/) {
				print {$out2_fhs[$mapkey]} ">$name,$seq\n";
			}
		}
		$name = "$1\/$3";
		$mapkey = get_mapkey($name);
		$seq = "";
	} elsif ($line =~ /^\+/){
		# is this a fastq quality line? eat chars to the length of the full sequence.
		while ($seqlen > 0) {
			$line = readline SEARCH_FH;
			chomp $line;
			$seqlen = $seqlen - length($line);
		}
	} else {
		$seq .= $line;
		$seqlen = length ($seq);
	}
}

if ($name =~ /\/1/) {
	print {$out1_fhs[$mapkey]} ">$name,$seq\n";
} elsif ($name =~ /\/2/) {
	print {$out2_fhs[$mapkey]} ">$name,$seq\n";
}
close SEARCH_FH;
my @pids = ();
printlog ("starting sort.");

for (my $i=0; $i<$numlibraries; $i++) {
	close $out1_fhs[$i];
	close $out2_fhs[$i];
	push @pids, fork_cmd ("sort -t',' -k 1 -T $tempdir @out1_bucketfiles[$i] > @out1_sortedfiles[$i]");
	if (@pids >= ($max_processes - 1)) {
		# don't spawn off too many threads at once.
		wait_for_forks(\@pids);
	}
}
wait_for_forks(\@pids);

for (my $i=0; $i<$numlibraries; $i++) {
    push @pids, fork_cmd ("sort -t',' -k 1 -T $tempdir @out2_bucketfiles[$i] > @out2_sortedfiles[$i]");
	if (@pids >= ($max_processes - 1)) {
		# don't spawn off too many threads at once.
		wait_for_forks(\@pids);
	}
}
wait_for_forks(\@pids);
printlog ("sorted.");

for (my $i=0; $i<$numlibraries; $i++) {
	printlog ("Making $output_file.$i.1.fasta.");
	open my $in_fh, "<", "@out1_sortedfiles[$i]" or exit_with_msg ("couldn't read @out1_sortedfiles[$i]");
	open my $out_fh, ">", "$output_file.$i.1.fasta" or exit_with_msg ("couldn't create $output_file.$i.1.fasta");
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;

	printlog ("Making $output_file.$i.2.fasta.");
	open my $in_fh, "<", "@out2_sortedfiles[$i]" or exit_with_msg ("couldn't read @out2_sortedfiles[$i]");
	open my $out_fh, ">", "$output_file.$i.2.fasta" or exit_with_msg ("couldn't create $output_file.$i.2.fasta");
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;
}

printlog ("Making blastdbs.");
for (my $i=0; $i<$numlibraries; $i++) {
	# make the blast db from the first of the paired end files
	push @pids, fork_cmd ("makeblastdb -in $output_file.$i.1.fasta -dbtype nucl -out $output_file.$i.db");
}

wait_for_forks(\@pids);

cleanup();
foreach my $tempfile (@tempfiles) {
	printlog ("removing $tempfile");
	system ("rm $tempfile");
}
printlog ("Finished");


sub get_mapkey {
	my $key = shift;
	$key =~ s/\/\d//;
	$key =~ s/#.+$//;
	$key =~ tr/0-9//dc;
	$key =~ /.*(\d{8})$/;
	$key = $1 * $multiplier;
	my $hash = $key % $numlibraries;
	return $hash;
}

__END__

=head1 NAME

makelibrary.pl

=head1 SYNOPSIS

makelibrary.pl -input short_read_archive [-output library_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and prepares it for aTRAM.

=head1 OPTIONS

 -input:   short read archive.
 -output:  optional: prefix of output library (default is the same as -input).

=cut

