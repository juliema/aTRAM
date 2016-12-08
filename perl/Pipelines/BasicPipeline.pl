#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "$0 " . join (" ", @ARGV) . "\n";

my $atrampath = "$FindBin::Bin/..";
my $help = 0;
my $samplefile = "";
my $targetfile = "";
my $outdir = ".";
my $log_file = "";
my $kmer = 31;
my $iterations = 20;
my $fraction = 1;
my $ins_length = 400;
my $debug = 0;
my $protein = 0;
my $complete = 0;
my $processes = 0;
my $max_memory = 0;

GetOptions ('samples=s' => \$samplefile,
            'targets=s' => \$targetfile,
            'kmer=i' => \$kmer,
            'iterations=i' => \$iterations,
			'fraction=f' => \$fraction,
			'ins_length=i' =>  \$ins_length,
			'output|outdir=s' => \$outdir,
			'debug|verbose' => \$debug,
			'protein' => \$protein,
			'complete' => \$complete,
			'processes=i' => \$processes,
			'max_memory|memory=i' => \$max_memory,
            'log_file=s' => \$log_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

set_debug ($debug);

$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);

if (($targetfile eq "") || ($samplefile eq "")) {
    pod2usage(-msg => "Must specify a list of aTRAM databases and a list of target sequences.", -verbose => 1);
}

if ($log_file eq "") {
	$log_file = File::Spec->catfile($outdir, "pipeline.log");
}

set_log($log_file);

printlog ("Running $runline");

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


# for each sample:
foreach my $sample (@samplenames) {
	my $sampledir = File::Spec->catfile($outdir, "$sample");
	make_path($sampledir);
	print "made $sampledir\n";
	my $targetcount = 0;
	foreach my $target (@targetnames) {
		$targetcount++;
		my $outname = File::Spec->catfile($sampledir, "$target");
		if (-e "$outname.results.txt") {
			printlog ("$targetcount: $target.results.txt already exists.");
			next;
		} else {
			printlog ("$targetcount: $target $sample");
			my $debug_flag = "";
			if ($debug == 1) { $debug_flag = "-debug"; }
			my $processes_flag = "";
			if ($processes > 0) { $processes_flag = "-processes $processes"; }
			my $protein_flag = "";
			if ($protein == 1) { $protein_flag = "-protein"; }
			my $complete_flag = "";
			if ($complete == 1) { $complete_flag = "-complete"; }
			my $memory_flag = "";
			if ($max_memory > 0) { $memory_flag = "-max_memory $max_memory"; }
			my $atram_result = run_command ("$atrampath/aTRAM.pl", "-reads $samples->{$sample} -target $targets->{$target} -iter $iterations -ins_length $ins_length -frac $fraction -assemble Velvet -out $outname -kmer $kmer $complete_flag $protein_flag $processes_flag $memory_flag $debug_flag -log $log_file", {"no_exit"=>1}); # don't exit because we want to capture a nonzero result.

			if ($atram_result == 255) {
				printlog ("aTRAM of $outname found no contigs.");
				next;
			} elsif ($atram_result != 0) {
				printlog ("unexpected error $atram_result.");
				die;
			}
		}
	}
}

printlog ("Finished executing $runline");




__END__

=head1 NAME

BasicPipeline.pl

=head1 SYNOPSIS

BasicPipeline.pl -samples samplefile -targets targetfile -output outputdir

Runs aTRAM on a list of aTRAM databases and a list of target sequences.

=head1 OPTIONS

 -samples:    tab-delimited list of aTRAM databases: for each row, "aTRAM_db_name   aTRAM_db_location".
 -targets:    tab-delimited list of target FASTA sequences: for each row, "gene_name   fasta_location".
 -kmer:       optional: kmer value for Velvet (default is 31).
 -iterations: optional: number of aTRAM iterations (default is 5).
 -fraction:   optional: fraction of aTRAM database to use (default is 1.0).
 -ins_length: optional: insert length of Illumina short read library (default is 400).
 -output:     output directory.

=cut

