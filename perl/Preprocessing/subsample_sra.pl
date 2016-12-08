#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $short_read_archive = "";
my $output_file = "";
my $help = 0;
my $frac = 0;
my $numlibraries = 1000;

GetOptions ('input=s' => \$short_read_archive,
            'output=s' => \$output_file,
            'fraction=f' => \$frac,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless ($short_read_archive) {
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
}

unless ($output_file) {
    pod2usage(-msg => "Must specify an output name for the subsampled fasta file.");
}

my $libsize = (-s $short_read_archive);

my $maxlib = $numlibraries * $frac;

my $libsizeMB = $libsize / 1e6;
$libsizeMB =~ s/(\d*)\.(\d{2}).*/\1.\2/;

my $subsetMB = int($libsizeMB * $frac);

if ($short_read_archive =~ /\.f.*q/) {
	$subsetMB = $subsetMB/2;
}

print "$short_read_archive is $libsizeMB MB. $output_file.fasta should be ~$subsetMB MB.\n";

my @primes = (29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151);
my $multiplier = $primes[$libsize % (@primes)];

my $name = "";
my $seq = "";
my $seqlen = 0;
my $key = 0;

open OUT_FH, ">", "$output_file.fasta";

open SEARCH_FH, "<:crlf", $short_read_archive;
while (my $line = readline SEARCH_FH) {
	chomp $line;
	if ($line =~ /^[@>](.*?)([\s\/])([12])/) {
		if ($name ne "") {
			if ($key < $maxlib) {
				if ($name =~ /\/1/) {
					print OUT_FH ">$name\n$seq\n";
				} elsif ($name =~ /\/2/) {
					print OUT_FH ">$name\n$seq\n";
				}
			}
		}
		$name = "$1\/$3";
		$key = subsample($name);
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

if ($key < $maxlib) {
	if ($name =~ /\/1/) {
		print OUT_FH ">$name\n$seq\n";
	} elsif ($name =~ /\/2/) {
		print OUT_FH ">$name\n$seq\n";
	}
}
close SEARCH_FH;
close OUT_FH;

print "Finished: $output_file.fasta is " . ((-s "$output_file.fasta") / 1e6) . " MB.\n";

sub subsample {
	my $key = shift;
	$key =~ s/\/\d//;
	$key =~ s/#.+$//;
	$key =~ tr/0-9//dc;
	$key =~ /.*(\d{8})$/;
	$key = $1 * $multiplier;
	return $key % $numlibraries;
}

__END__

=head1 NAME

subsample_sra.pl

=head1 SYNOPSIS

subsample_sra.pl -input short_read_archive -output library_name -frac fraction

Takes a random fraction of a fasta or fastq file of paired-end short reads.

=head1 OPTIONS

 -input:   short read archive.
 -output:  prefix of output library.
 -frac:    fraction to subsample.

=cut

