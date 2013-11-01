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
	print "" . timestamp() . ": fastq file inputted...converting to fasta.\n";
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

my $tempfile1 = "$working_sra.temp.1";
my $tempfile2 = "$working_sra.temp.2";
my $tempdir = dirname ("$working_sra");
my $outfile = "$working_sra.sorted.fasta";

# sort fasta short-read file
print "" . timestamp() . ": sorting fasta file.\n";

open SEARCH_FH, "<", $working_sra;
open TEMP1_FH, ">", $tempfile1;
open TEMP2_FH, ">", $tempfile2;
my $name = "";
my $seq = "";
my $line;
while ($line = readline SEARCH_FH) {
	chomp $line;
	if ($line =~ />(.*)/) {
		if ($name ne "") {
			if ($line =~ /\/1/) {
				print TEMP1_FH ">$name,$seq\n";
			} elsif ($name =~ /\/2/) {
				print TEMP2_FH ">$name,$seq\n";
			}
		}
		$name = $1;
		$seq = "";
	} else {
		$seq .= $line;
	}
}

if ($line =~ /\/1/) {
	print TEMP1_FH ">$name,$seq\n";
} elsif ($name =~ /\/2/) {
	print TEMP2_FH ">$name,$seq\n";
}
close SEARCH_FH;
close TEMP1_FH;
close TEMP2_FH;

my $tempsort1 = "$working_sra.sort.1";
my $tempsort2 = "$working_sra.sort.2";

print "" . timestamp() . ": running sort on $tempsort1.\n";
system ("sort -t',' -k 1 --parallel=8 -T $tempdir $tempfile1 > $tempsort1");
print "" . timestamp() . ": running sort on $tempsort2.\n";
system ("sort -t',' -k 1 --parallel=8 -T $tempdir $tempfile2 > $tempsort2");
print "" . timestamp() . ": sorted.\n";
system ("rm $tempfile1");
system ("rm $tempfile2");

if (-z "$tempsort1") {
	exit_with_msg ("Sort failed. Are you sure $working_sra exists?");
}
# un-interleave fasta file into two paired files:
print "" . timestamp() . ": validating paired files.\n";
# open FH, "<", "$outfile" or exit_with_msg ("couldn't open fasta file");
open FH1, "<", "$tempsort1";
open FH2, "<", "$tempsort2";

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
	my $line1 = readline FH1;
	chomp $line1;
	my ($name1, $seq1) = split(/,/,$line1);
	if (($name1 eq "") || ($seq1 eq "")) { last; }
	if ($name1 =~ />(.*?)\/1/) {
		my $seqname = $1;
		my $line2 = readline FH2;
		chomp $line2;
		my ($name2, $seq2) = split(/,/,$line2);
		if (($name2 eq "") || ($seq2 eq "")) { last; }
		my $i = $count % $numlibraries;
		$count++;
		my $out1_fh = $out1_fhs[$i];
		my $out2_fh = $out2_fhs[$i];
		print $out1_fh "$name1\n$seq1\n";
		print $out2_fh "$name2\n$seq2\n";
	} else {
		# this sequence doesn't have the name format of a paired read. Skip it.
		next;
	}
}

close FH1;
close FH2;

for (my $i=0; $i<$numlibraries; $i++) {
	my $out1_fh = $out1_fhs[$i];
	my $out2_fh = $out2_fhs[$i];
	close $out1_fh;
	close $out2_fh;
}

for (my $i=1; $i<=$numlibraries; $i++) {
	# make the blast db from the first of the paired end files
	print "" . timestamp() . ": Making blastdb from first of paired fasta files.\n";
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

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return "$date $time";
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

