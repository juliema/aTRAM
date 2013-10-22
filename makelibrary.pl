#!/usr/bin/perl
use strict;
use File::Basename;


if (@ARGV == 0) {
	exit_with_msg ("Usage: 0-prepare_files.pl short_read_archive.fasta/fastq");
}
my $short_read_archive = shift;

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

my $name1 = readline FH;
while ($name1) {
	my $seq1 = readline FH;
	my $name2 = readline FH;
	my $seq2 = readline FH;
	$name1 =~ />(.*?)\/1/;
	my $name = $1;
	if ($name2 !~ /$name\/2/) {
		exit_with_msg ("$working_sra.sorted.fasta does not contain paired reads: read $name\/1 is not followed by $name\/2 at line " . ($. - 4));
	}
	print OUT1_FH $name1 . $seq1;
	print OUT2_FH $name2 . $seq2;
	$name1 = readline FH;
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
