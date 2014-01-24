#!usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
require Subfunctions;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $help = 0;
my $samplefile = "";
my $targetfile = "";
my $outfile = "";
my $kmer = 31;
my $iter = 10;
my $frac = 0.01;
my $ins_length = 400;

GetOptions ('samples=s' => \$samplefile,
            'targets=s' => \$targetfile,
            'kmer=i' => \$kmer,
            'iter=i' => \$iter,
			'frac=f' => \$frac,
			'ins_length=i' =>  \$ins_length,
			'output|outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($outfile eq "") {
	$outfile = "result";
}


if (($targetfile eq "") || ($samplefile eq "")) {
    pod2usage(-msg => "Must specify a list of aTRAM databases and a list of target sequences.", -verbose => 1);
	exit;
}

my $samples = {};
my @samplenames = ();
open FH, "<", "$samplefile" or die "couldn't open $samplefile";
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$samples->{$1} = $2;
		push @samplenames, $1;
	}
}
close FH;

my $targets = {};
my @targetnames = ();
open FH, "<", $targetfile;
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$targets->{$1} = $2;
		push @targetnames, $1;
	}
}
close FH;

open LOG_FH, ">", "$outfile.log";
print LOG_FH "target\tsample\tcontig\tbitscore\tpercentcoverage\n";
foreach my $target (@targetnames) {
	open FH, ">", "$outfile.$target.exons.fasta";
	truncate FH, 0;
	close FH;
	open FH, ">", "$outfile.$target.full.fasta";
	truncate FH, 0;
	close FH;

	# for each sample:
	foreach my $sample (@samplenames) {
		my $outname = "$outfile.$target.$sample";
		print "$target $sample\n";
		system_call ("perl ~/TRAM/aTRAM.pl -reads $samples->{$sample} -target $targets->{$target} -iter $iter -ins_length $ins_length -frac $frac -assemble Velvet -out $outname -kmer $kmer -complete");
		system_call ("rm -r $outname.Velvet");
		# run percentcoverage to get the contigs nicely aligned
		system_call ("perl ~/TRAM/test/PercentCoverage.pl $targets->{$target} $outname.best.fasta $outname");

		# find the one best contig (one with fewest gaps)
		system_call ("blastn -task blastn -query $outname.exons.fasta -subject $targets->{$target} -outfmt '6 qseqid bitscore' -out $outname.blast");
		open FH, "<", "$outname.blast";
		my $contig = "";
		my $score = 0;
		foreach my $line (<FH>) {
			if ($line =~ /(\S+)\s+(\S+)$/) {
				if ($1 =~ /reference/) {
					next;
				}
				if ($2 > $score) {
					$contig = $1;
					$score = $2;
				}
			}
		}
		close FH;
		system_call ("rm $outname.blast");

		open FH, "<", "$outname.Table.txt";
		$contig = "";
		my $percent = 0;
		foreach my $line (<FH>) {
			if ($line =~ /(\S+)\t(\S+)\t(\S+)$/) {
				if ($1 =~ /$target/) {
					next;
				}
				if ($3 > $percent) {
					$contig = $1;
					$percent = $3;
				}
			}
		}
		close FH;

		$percent =~ s/^(\d+\.\d{2}).*$/\1/;
		print LOG_FH "$target\t$sample\t$contig\t$score\t$percent\n";
		if ($contig ne "") {
			# pick this contig from the fasta file
			my ($taxa, $taxanames) = parse_fasta ("$outname.exons.fasta");
			# write this contig out to the target.fasta file, named by sample.
			open FH, ">>", "$outfile.$target.exons.fasta";
			print "adding $contig to $target.exons.fasta\n";
			print FH ">$sample\n$taxa->{$contig}\n";
			close FH;
			($taxa, $taxanames) = parse_fasta ("$outname.best.fasta");
			# write this contig out to the target.fasta file, named by sample.
			open FH, ">>", "$outfile.$target.full.fasta";
			print "adding $contig from $outname.best.fasta to $target.full.fasta\n";
			print FH ">$sample\n$taxa->{$contig}\n";
			close FH;
		}
	}
}

close LOG_FH;



__END__

=head1 NAME

AlignmentPipeline.pl

=head1 SYNOPSIS

AlignmentPipeline.pl -input short_read_archive [-output aTRAM_db_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and creates an aTRAM database for use by aTRAM.pl.

=head1 OPTIONS

 -samples:    tab-delimited list of aTRAM databases: for each row, "aTRAM_db_name   aTRAM_db_location".
 -targets:    tab-delimited list of target FASTA sequences: for each row, "gene_name   fasta_location".
 -kmer:       optional: kmer value for Velvet (default is 31).
 -iter:       optional: number of aTRAM iterations (default is 5).
 -frac:       optional: fraction of aTRAM database to use (default is 1.0).
 -ins_length: optional: insert length of Illumina short read library (default is 400).
 -output:   optional: prefix of aTRAM database (default is the same as -input).

=cut

