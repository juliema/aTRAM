#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $atrampath = "$FindBin::Bin/..";
my $help = 0;
my $samplefile = "";
my $targetfile = "";
my $outdir = ".";
my $kmer = 31;
my $iter = 20;
my $frac = 1;
my $ins_length = 400;
my $debug = 0;
my $protein = 0;
my $complete = 0;

GetOptions ('samples=s' => \$samplefile,
            'targets=s' => \$targetfile,
            'kmer=i' => \$kmer,
            'iter=i' => \$iter,
			'frac=f' => \$frac,
			'ins_length=i' =>  \$ins_length,
			'output|outdir=s' => \$outdir,
			'debug|verbose' => \$debug,
			'protein' => \$protein,
			'complete' => \$complete,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

set_debug ($debug);

$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);

if (($targetfile eq "") || ($samplefile eq "")) {
    pod2usage(-msg => "Must specify a list of aTRAM databases and a list of target sequences.", -verbose => 1);
	exit;
}

open my $log_fh, ">", File::Spec->catfile($outdir, "pipeline.log");;
set_log($log_fh);

my $samples = {};
my @samplenames = ();
open FH, "<", "$samplefile" or die "Couldn't open sample file $samplefile";
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$samples->{$1} = $2;
		push @samplenames, $1;
	}
}
close FH;

if (@samplenames == 0) {
	die "Sample file $samplefile doesn't contain a list";
}

my $targets = {};
my @targetnames = ();
open FH, "<", $targetfile or die "Couldn't open target file $targetfile";
foreach my $line (<FH>) {
	if ($line =~ /(.+?)\t(.+)/) {
		$targets->{$1} = $2;
		push @targetnames, $1;
	}
}
close FH;

if (@targetnames == 0) {
	die "Target file $targetfile doesn't contain a list";
}

# for each sample:
foreach my $sample (@samplenames) {
	my $sampledir = File::Spec->catfile($outdir, "$sample");
	make_path($sampledir);
	print "made $sampledir\n";
	foreach my $target (@targetnames) {
		my $outname = File::Spec->catfile($sampledir, "$target");
		if (-e "$outname.results.txt") {
			printlog ("$target $sample");
			my $protein_flag = "";
			my $complete_flag = "";
			if ($protein == 1) { $protein_flag = "-protein"; }
			if ($complete == 1) { $complete_flag = "-complete"; }
			my $atram_result = system_call ("perl $atrampath/aTRAM.pl -reads $samples->{$sample} -target $targets->{$target} -iter $iter -ins_length $ins_length -frac $frac -assemble Velvet -out $outname -kmer $kmer $complete_flag $protein_flag", 1);

			if ($atram_result == 255) {
				printlog ("$atram_result: aTRAM of $outname found no contigs.");
				next;
			} elsif ($atram_result != 0) {
				die "unexpected error $atram_result";
			}
		} else {
			printlog ("$target.results.txt already exists.");
		}
	}
}


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
 -iter:       optional: number of aTRAM iterations (default is 5).
 -frac:       optional: fraction of aTRAM database to use (default is 1.0).
 -ins_length: optional: insert length of Illumina short read library (default is 400).
 -output:     output directory.

=cut

