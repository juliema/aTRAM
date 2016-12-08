#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;
use Mapreduce;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $short_read_archive = "";
my $help = 0;
my $numshards = 0;
my $srasize = 0;
my $multiplier = 0;
my $debug = 0;

GetOptions ('input=s' => \$short_read_archive,
            'number=i' => \$numshards,
            'size=i' => \$srasize,
            'multiplier=i' => \$multiplier,
            'debug|verbose' => \$debug,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive) {
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
}

if ($srasize == 0) {
	$srasize = (-s $short_read_archive);
}

my $srasizeMB = $srasize / 1e6;
$srasizeMB =~ s/(\d*)\.(\d{2}).*/$1.$2/;

if ($multiplier == 0) {
	$multiplier = set_multiplier ($srasize);
	print ("Using multiplier $multiplier based on size $srasize.\n");
} else {
	$multiplier = set_multiplier ($multiplier);
	print ("Using multiplier $multiplier based on user input.\n");
}

# if the user didn't specify how many shards to make, we should make as many as we need so that they average 500MB each.
if ($numshards == 0) {
	# if it's a fastq file, the file is twice the size that it would be if it were a fasta.
	if ($short_read_archive =~ /\.f.*q/) {
		$srasize = $srasize/2;
	}
	$numshards = int($srasize / 5e8);
	print ("$short_read_archive is $srasizeMB MB; we will make $numshards shards.\n");
}

# declare how many shards we'll be making.
my $total_shards = set_total_shards ($numshards);
if ($total_shards == 0) {
	print "No shards to be made.\n";
	exit;
} else {
	print ("Making $total_shards shards.\n");
}

open FH, "<:crlf", $short_read_archive;
my @keys = ();
for (my $i=0;$i<$numshards; $i++) {
	$keys[$i] = 0;
}

my $totalkeys = 0;
foreach my $line (<FH>) {
	if ($line =~ /^[@>](.*)$/) {
		my $name = $1;
		my ($key, $calc) = map_to_shard ($name, $debug);
		if ($debug) {
			print "$name maps to $key: $calc\n";
		}
		$keys[$key]++;
		$totalkeys++;
	}
}

for (my $i=0;$i<$numshards; $i++) {
	print "shard $i has\t$keys[$i] keys\n";
}

print "average shard should be " . int ($totalkeys / $numshards) . " keys in size.\n\n";


close FH;


__END__

=head1 NAME

test_format_sra.pl

=head1 SYNOPSIS

test_format_sra.pl -input short_read_archive [-output aTRAM_db_name] [-number int]

Tests format_sra for a particular set of Illumina names. Takes various parameters, but at a minimum, an input file with fasta/q names.

=head1 OPTIONS

 -input:      short read archive.
 -number:     optional: number of shards to create (default is however many are required for each to be ~500MB).
 -size:       optional: the actual size in bytes of the sra, if the input file is a subset of the sra.
 -multiplier: optional: a number to be used as the mapping multiplier.

=cut

