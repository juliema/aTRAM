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
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
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

my @out1_fhs = ();
my @out2_fhs = ();

for (my $i=1; $i<=$numlibraries; $i++) {
	open my $out1_fh, ">", "$working_sra.$i.1.fasta" or exit_with_msg ("couldn't create result file");
	open my $out2_fh, ">", "$working_sra.$i.2.fasta" or exit_with_msg ("couldn't create result file");
	push @out1_fhs, $out1_fh;
	push @out2_fhs, $out2_fh;
}

my $count = 0;
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
			# these two sequences are a pair: write them out.
			my $i = $count % $numlibraries;
			$count++;
			my $out1_fh = $out1_fhs[$i];
			my $out2_fh = $out2_fhs[$i];
			print $out1_fh "$name1\n$seq1\n";
			print $out2_fh "$name2\n$seq2\n";
		} else {
			# the next sequence isn't the pair of the first sequence. Skip these.
			next;
		}
	} else {
		# this sequence doesn't have the name format of a paired read. Skip it.
		next;
	}
}

close FH;

for (my $i=0; $i<$numlibraries; $i++) {
	my $out1_fh = $out1_fhs[$i];
	my $out2_fh = $out2_fhs[$i];
	close $out1_fh;
	close $out2_fh;
}

for (my $i=1; $i<=$numlibraries; $i++) {
	# make the blast db from the first of the paired end files
	print "Making blastdb from first of paired fasta files.\n";
	system ("makeblastdb -in $working_sra.$i.1.fasta -dbtype nucl -out $working_sra.$i.db");
}

system ("rm $working_sra.sorted.fasta");
if ($fastq_input == 1) {
	system ("rm $working_sra");
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

