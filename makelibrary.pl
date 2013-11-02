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

my $tempfile1 = "$short_read_archive.temp.1";
my $tempfile2 = "$short_read_archive.temp.2";
my $tempdir = dirname ("$short_read_archive");
my $outfile = "$short_read_archive.sorted.fasta";

# sort fasta short-read file
print "" . timestamp() . ": sorting fasta file.\n";

open SEARCH_FH, "<", $short_read_archive;
open TEMP1_FH, ">", $tempfile1;
open TEMP2_FH, ">", $tempfile2;
my $name = "";
my $seq = "";
my $line;
my $seqlen = 0;
while ($line = readline SEARCH_FH) {
	chomp $line;
	if ($line =~ /[@>](.*?)([\s\/])([12])/) {
		if ($name ne "") {
			if ($name =~ /\/1/) {
				print TEMP1_FH ">$name,$seq\n";
			} elsif ($name =~ /\/2/) {
				print TEMP2_FH ">$name,$seq\n";
			}
		}
		hash_to_bucket($name);
		$name = "$1\/$3";
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
	print TEMP1_FH ">$name,$seq\n";
} elsif ($name =~ /\/2/) {
	print TEMP2_FH ">$name,$seq\n";
}
close SEARCH_FH;
close TEMP1_FH;
close TEMP2_FH;
my $tempsort1 = "$short_read_archive.sort.1";
my $tempsort2 = "$short_read_archive.sort.2";

print "" . timestamp() . ": running sort on $tempsort1.\n";
system ("sort -t',' -k 1 -T $tempdir $tempfile1 > $tempsort1");
print "" . timestamp() . ": running sort on $tempsort2.\n";
system ("sort -t',' -k 1 -T $tempdir $tempfile2 > $tempsort2");
print "" . timestamp() . ": sorted.\n";
system ("rm $tempfile1");
system ("rm $tempfile2");

if (-z "$tempsort1") {
	exit_with_msg ("Sort failed. Are you sure $short_read_archive exists?");
}
# un-interleave fasta file into two paired files:
print "" . timestamp() . ": validating paired files.\n";
# open FH, "<", "$outfile" or exit_with_msg ("couldn't open fasta file");
open FH1, "<", "$tempsort1";
open FH2, "<", "$tempsort2";

my @out1_fhs = ();
my @out2_fhs = ();

for (my $i=1; $i<=$numlibraries; $i++) {
	open my $out1_fh, ">", "$short_read_archive.$i.1.fasta" or exit_with_msg ("couldn't create result file");
	open my $out2_fh, ">", "$short_read_archive.$i.2.fasta" or exit_with_msg ("couldn't create result file");
	push @out1_fhs, $out1_fh;
	push @out2_fhs, $out2_fh;
}

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
		my $i = hash_to_bucket($name);
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
	print "" . timestamp() . ": Making blastdb from $short_read_archive.$i.1.fasta.\n";
	system ("makeblastdb -in $short_read_archive.$i.1.fasta -dbtype nucl -out $short_read_archive.$i.db");
}

system ("rm $tempsort1");
system ("rm $tempsort2");

print "" . timestamp() . ": Finished\n";

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

sub hash_to_bucket {
	my $key = shift;
	my $bucket;

	if ($key =~ /(.+)\/\d/) {
		$key = $1;
	}

	$key =~ s/\D//g;
	$bucket = $key % $numlibraries;
	return $bucket;
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

