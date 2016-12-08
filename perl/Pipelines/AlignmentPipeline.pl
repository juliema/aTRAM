#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;
use Parsing;
use Postprocessing;
use Configuration;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "$0 " . join (" ", @ARGV) . "\n";

my $atrampath = "$FindBin::Bin/..";
my $help = 0;
my $samplefile = "";
my $targetfile = "";
my $outdir = ".";
my $kmer = 31;
my $iterations = 10;
my $fraction = 1;
my $ins_length = 400;
my $debug = 0;
my $protein = 0;
my $complete = 0;
my $processes = 0;
my $assembler = "";
my $aligner = "";
my $log_file = "";
my $max_memory = 0;

GetOptions ('samples=s' => \$samplefile,
            'targets=s' => \$targetfile,
            'kmer=i' => \$kmer,
            'iterations=i' => \$iterations,
			'fraction=f' => \$fraction,
			'ins_length=i' =>  \$ins_length,
			'output|outdir=s' => \$outdir,
			'processes=i' => \$processes,
			'debug|verbose' => \$debug,
			'complete' => \$complete,
			'assembler=s' => \$assembler,
			'aligner=s' => \$aligner,
			'max_memory|memory=i' => \$max_memory,
            'log_file=s' => \$log_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($outdir eq "") {
	$outdir = "result";
}

set_debug ($debug);

if (($targetfile eq "") || ($samplefile eq "")) {
    pod2usage(-msg => "Must specify a list of aTRAM databases and a list of target sequences.", -verbose => 1);
	exit 2;
}

$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);

if ($log_file eq "") {
	$log_file = File::Spec->catfile($outdir, "pipeline.log");
}

set_log($log_file);

printlog ("Running $runline");

Configuration::initialize();

my $samples = {};
my @samplenames = ();
open FH, "<:crlf", "$samplefile" or die "Couldn't open sample file $samplefile";
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$samples->{$1} = $2;
		push @samplenames, $1;
	}
}
close FH;

if (@samplenames == 0) {
	printlog ("Sample file $samplefile doesn't contain a list.");
	die;
}

my $targets = {};
my @targetnames = ();
open FH, "<:crlf", $targetfile or die "Couldn't open target file $targetfile";
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$targets->{$1} = $2;
		push @targetnames, $1;
	}
}
close FH;

if (@targetnames == 0) {
	printlog ("Target file $targetfile doesn't contain a list.");
	die;
}

my $atram_dir = File::Spec->catfile($outdir, "aTRAM");
make_path($atram_dir);
my $align_dir = File::Spec->catfile($outdir, "alignments");
make_path($align_dir);

open TABLE_FH, ">", File::Spec->catfile($outdir, "results.txt");
print TABLE_FH "target\tsample\tcontig\tbitscore\tpercentcoverage\n";
foreach my $target (@targetnames) {
	my ($targetseqs, undef) = parsefasta ("$targets->{$target}");
	my $refseq = join ("",(values %$targetseqs));
	if (is_protein($refseq)) {
		# if any of the seqs in the target fasta are protein sequences, break.
		printlog ("Target fasta file $targets->{$target} is not a DNA sequence, cannot perform alignment.");
		next;
	}

	my $trimmed_file = File::Spec->catfile($align_dir, "$target.trimmed.fasta");
	open FH, ">", $trimmed_file;
	truncate FH, 0;
	print FH ">reference\n$refseq\n";
	close FH;

	my $full_contigs_file = File::Spec->catfile($align_dir, "$target.full.fasta");
	open FH, ">", $full_contigs_file;
	truncate FH, 0;
	close FH;

	# for each sample:
	foreach my $sample (@samplenames) {
		my $outname = "$target.$sample";
		printlog ("Processing $outname...");

		my $memory_flag = "";
		if ($max_memory > 0) { $memory_flag = "-max_memory $max_memory"; }

		my $complete_flag = "";
		if ($complete == 1) { $complete_flag = "-complete"; }

		my $processes_flag = "";
		if ($processes > 0) { $processes_flag = "-processes $processes"; }

		my $debug_flag = "";
		if ($debug == 1) { $debug_flag = "-debug"; }

		if ($assembler ne "") { $assembler = "-assemble $assembler"; }

		my $atram_outname = File::Spec->catfile($atram_dir, $outname);
		my $atram_result = run_command ("$atrampath/aTRAM.pl", "-reads $samples->{$sample} -target $targets->{$target} -iter $iterations -ins_length $ins_length -frac $fraction $assembler -out $atram_outname -kmer $kmer $complete_flag $processes_flag $memory_flag $debug_flag -log $log_file", {"no_exit"=>1}); # don't exit because we want to capture a nonzero result.

		if ($atram_result) {
			printlog ("aTRAM found no contigs matching $target for $sample.");
			next;
		}
		# run percentcoverage to get the contigs nicely aligned
		my $comparefile = "$atram_outname.best.fasta";
		if (($complete == 1) && (-e "$atram_outname.complete.fasta")) {
			$comparefile = "$atram_outname.complete.fasta";
		}

		my $contigs = percentcoverage ($targets->{$target}, $comparefile, $atram_outname, $aligner);

		if (!(defined $contigs)) {
			printlog ("percentcoverage failed for $targets->{$target}, $comparefile.\n");
			exit -1;
		}

		my $refseq = delete $contigs->{reference};

		open EXON_FH, ">", "$atram_outname.trimmed.fasta";
		my $results = {};
		my @result_names = ();

		print EXON_FH ">reference\n$refseq\n";
		my $total_length = length $refseq;
		foreach my $contig (keys %$contigs) {
			print EXON_FH ">$contig\n$contigs->{$contig}\n";
			my $gaps = ($contigs->{$contig} =~ tr/N-//);
			my $total = $total_length - $gaps;
			my $percent = $total / $total_length;
			push @result_names, $contig;
			$results->{$contig}->{total} = $total;
			$results->{$contig}->{percent} = $percent;
		}

		close EXON_FH;

		# find the one best contig (one with fewest gaps)
		if ($protein == 0) {
			run_command (get_bin("blastn"), "-task blastn -query $atram_outname.trimmed.fasta -subject $targets->{$target} -outfmt '6 qseqid bitscore' -out $atram_outname.blast");
		}
		open FH, "<:crlf", "$atram_outname.blast";
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
		`rm $atram_outname.blast`;

		$contig = "";
		my $percent = 0;
		foreach my $c (@result_names) {
			if ($results->{$c}->{percent} > $percent) {
				$contig = $c;
				$percent = $results->{$c}->{percent};
			}
		}

		$percent =~ s/^(\d+\.\d{2}).*$/\1/;
		print TABLE_FH "$target\t$sample\t$contig\t$score\t$percent\n";
		if ($contig ne "") {
			# pick this contig from the fasta file
			my ($taxa, $taxanames) = parsefasta ("$atram_outname.trimmed.fasta");
			# write this contig out to the target.fasta file, named by sample.
			open FH, ">>", $trimmed_file;
			printlog ("adding $contig to $target.trimmed.fasta");
			print FH ">$sample\n$taxa->{$contig}\n";
			close FH;
			($taxa, $taxanames) = parsefasta ("$atram_outname.best.fasta");
			# write this contig out to the target.fasta file, named by sample.
			open FH, ">>", $full_contigs_file;
			printlog ("adding $contig from $atram_outname.best.fasta to $target.full.fasta");
			print FH ">$sample\n$taxa->{$contig}\n";
			close FH;
		}
	}
}

close TABLE_FH;

printlog ("Finished executing $runline");


__END__

=head1 NAME

AlignmentPipeline.pl

=head1 SYNOPSIS

AlignmentPipeline.pl -samples samplefile -targets targetfile -output outputprefix

Runs aTRAM on a list of aTRAM databases and a list of target sequences, returning a
fasta file for each target sequence with the best contig for each aTRAM database.

For each target sequence, the best match is determined by aligning the best contigs with
the target sequence to determine how much of the target sequence is covered, and then
re-blasting the aligned contigs against the target sequence to find the best-scoring match.

=head1 OPTIONS

 -samples:    tab-delimited list of aTRAM databases: for each row, "aTRAM_db_name   aTRAM_db_location".
 -targets:    tab-delimited list of target FASTA sequences: for each row, "gene_name   fasta_location".
 -kmer:       optional: kmer value for Velvet (default is 31).
 -iterations: optional: number of aTRAM iterations (default is 5).
 -fraction:   optional: fraction of aTRAM database to use (default is 1.0).
 -ins_length: optional: insert length of Illumina short read library (default is 400).
 -output:     output file prefix.

=cut

