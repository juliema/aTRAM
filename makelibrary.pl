#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;


if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $short_read_archive = "";
my $output_file = "";
my $numlibraries = 1;
my $help = 0;

GetOptions ('input=s' => \$short_read_archive,
            'output=s' => \$output_file,
            'number=i' => \$numlibraries,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with makelibrary.pl) and a target gene in fasta form.");
}

unless ($output_file) {
    $output_file = $short_read_archive;
}

my $executing_path = dirname(__FILE__);

my $fastq_input = 0;

# if the sra is a fastq file, make it fasta.
if ($short_read_archive =~ /\.f.*q$/) { # if it's a fastq file:
	print "fastq file inputted...converting to fasta.\n";
	$fastq_input = 1;
	# un-interleave fastq file into fasta:
	open FH, "<", $short_read_archive or exit_with_msg("couldn't open fasta file");

	open OUT_FH, ">", "$short_read_archive.fasta" or exit_with_msg ("couldn't create result file");

	my $fs = readline FH;
	while ($fs) {
		$fs =~ s/@/>/;
		print OUT_FH $fs;
		$fs = readline FH;
		print OUT_FH $fs;
		$fs = readline FH;
		# toss quality lines
		if ($fs !~ /^\+/) {
			print "$fs";
			exit_with_msg ("$short_read_archive is not a properly-formed fastq file: line ". $. ." is problematic.");
		}
		$fs = readline FH;
		$fs = readline FH;
	}

	close FH;
	close OUT1_FH;

}

# now we're working with a fasta file for sure.
my $working_sra = $short_read_archive;
if ($fastq_input == 1) {
	$working_sra = "$short_read_archive.fasta";
}

# sort fasta short-read file
print "sorting fasta file.\n";
system ("bash $executing_path/lib/sort_fasta.sh $working_sra") == 0 or exit_with_msg ("Couldn't find sort_fasta.sh. This script needs to be in the same directory as the rest of TRAM");

if (-z "$working_sra.sorted.fasta") {
	exit_with_msg ("Sort failed. Are you sure $working_sra exists?");
}
# un-interleave fasta file into two paired files:
print "un-interleaving $working_sra.sorted.fasta into paired files.\n";
open FH, "<", "$working_sra.sorted.fasta" or exit_with_msg ("couldn't open fasta file");

open OUT1_FH, ">", "$working_sra.1.fasta" or exit_with_msg ("couldn't create result file");
open OUT2_FH, ">", "$working_sra.2.fasta" or exit_with_msg ("couldn't create result file");

while (1) {
	my $name1 = readline FH;
	my $seq1 = readline FH;
	chomp $name1;
	chomp $seq1;
	if (($name1 eq "") || ($seq1 eq "")) { last; }
	if ($name1 =~ />(.*?)\/1/) {
		my $seqname = $1;
		my $name2 = readline FH;
		my $seq2 = readline FH;
		chomp $name2;
		chomp $seq2;
		if (($name2 eq "") || ($seq2 eq "")) { last; }
		if ($name2 =~ /$seqname\/2/) {
# 			print "$name1 and $name2 are a pair\n";
			print OUT1_FH "$name1\n$seq1\n";
			print OUT2_FH "$name2\n$seq2\n";
		} else {
# 			print "$name2 wasn't the pair of $name1\n";
			next;
		}
	} else {
		# print "$name1 doesn't seem to be a paired read.";
		next;
	}
}


close FH;
close OUT1_FH;
close OUT2_FH;

# make the blast db from the first of the paired end files
print "Making blastdb from first of paired fasta files.\n";
system ("makeblastdb -in $working_sra.1.fasta -dbtype nucl -out $working_sra.db");

if (((-s "$working_sra.1.fasta") + (-s "$working_sra.2.fasta")) == (-s "$working_sra")) {
	print "Files prepared successfully, cleaning up.\n";
	system ("rm $working_sra.sorted.fasta; rm $working_sra;");
} else {
	exit_with_msg ("Something went wrong; leaving intermediate files alone.\n");
}

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

__END__

=head1 NAME

sTRAM.pl

=head1 SYNOPSIS

makelibrary.pl -input short_read_archive [-output library_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and prepares it for sTRAM.

=head1 OPTIONS

 -input:   short read archive.
 -output:  optional: prefix of output library (default is the same as -input).
 -number:  optional: the number of partitioned libraries to make.

=cut

