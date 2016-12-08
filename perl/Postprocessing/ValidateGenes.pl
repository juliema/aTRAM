#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec qw (rel2abs);
use File::Temp qw (tempfile);
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;
use Configuration;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "$0 " . join (" ", @ARGV) . "\n";

my $atrampath = "$FindBin::Bin/..";
my $help = 0;
my $inputdir = "";
my $reffile = "";
my $outdir = ".";
my $debug = 0;
my $processes = 1;
my $blastdb = "";

GetOptions ('input=s' => \$inputdir,
            'reference=s' => \$reffile,
			'output|outdir=s' => \$outdir,
			'debug|verbose' => \$debug,
			'processes=i' => \$processes,
			'database|db=s' => \$blastdb,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if (!(-e $inputdir)) {
    pod2usage(-msg => "Must specify a directory of fasta files.", -verbose => 1);
}

if (!(-e $reffile)) {
    pod2usage(-msg => "Must specify a valid reference fasta file.", -verbose => 1);
}

if ($blastdb ne "") {
	my $info = `blastdbcmd -db $blastdb -info`;
	$info =~ /Database: (.+?)\n/;
	if ($reffile !~ /$1/) {
		pod2usage(-msg => "$1 is not a database created from $reffile.", -verbose => 1);
	}
}

set_debug ($debug);

$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);

my($refname, undef, undef) = File::Basename::fileparse($reffile, qw(.fasta .fa));

Configuration::initialize();

my $log_file = File::Spec->catfile($outdir, "pipeline.log");
set_log($log_file);

printlog ("Running $runline");

my $crossfile = "$blastdb.crossref";
if ($blastdb eq "") {
	$blastdb = File::Spec->catfile($outdir, "$refname.blastdb");
	$crossfile = File::Spec->catfile($outdir, "$refname.blastdb.crossref");
	run_command (get_bin('makeblastdb'), "-in $reffile -dbtype nucl -out $blastdb");
	run_command (get_bin('blastn'), "-query $reffile -db $blastdb -num_threads $processes -evalue 1e-50 -out $crossfile -outfmt '6 qseqid sseqid evalue length'");
}

open CROSS_FH, "<:crlf", $crossfile;
my $curr_gene = "";
my $cross_matrix = {};
my @cross_keys = ();
foreach my $line (<CROSS_FH>) {
	if ($line =~ /(.+?)\t(.+?)\t(.+?)\t(.+)$/) {
		my ($qseqid, $sseqid, $evalue, $length) = ($1, $2, $3, $4);
		if ($qseqid ne $curr_gene) {
			$curr_gene = $qseqid;
			push @cross_keys, $curr_gene;
			$cross_matrix->{$curr_gene}->{$sseqid}->{evalue} = $evalue;
			$cross_matrix->{$curr_gene}->{$sseqid}->{"length"} = $length;
		} else {
			if (!(exists $cross_matrix->{$curr_gene}->{$sseqid})) {
				$cross_matrix->{$curr_gene}->{$sseqid}->{evalue} = $evalue;
				$cross_matrix->{$curr_gene}->{$sseqid}->{"length"} = $length;
			}
		}
	}
}

opendir (my $INPUT_DH, "$inputdir");
my @contigfiles = readdir $INPUT_DH;
closedir $INPUT_DH;
print( "opening $inputdir, there were " . @contigfiles . " files\n");

open RESULTS_FH, ">", File::Spec->catfile($outdir, "results.txt");
foreach my $contigfile (@contigfiles) {
	my $filepath = File::Spec->rel2abs (File::Spec->catfile ($inputdir, $contigfile));
	if ($contigfile =~ /(.+)\.(.+)\.fasta/) {
		my $genename = $1;
		my $blastfile = File::Spec->rel2abs (File::Spec->catfile ($outdir, "$contigfile.blast"));
		run_command (get_bin('blastn'), "-query $filepath -db $blastdb -num_threads $processes -evalue 1e-50 -out $blastfile -outfmt '6 sseqid qseqid slen evalue length pident'");
		my $matchmatrix = {};
		my @matrixkeys = ();
		open FH, "<:crlf", $blastfile;
		my $line = readline FH;
		if ($line =~ /(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+)$/) {
			my ($sseqid, $qseqid, $slen, $evalue, $length, $pident) = ($1, $2, $3, $4, $5, $6);
			my $hit_string = "";
			if ($genename eq $sseqid) {
				$hit_string = "self\t$qseqid";
				print "best hit for $genename is itself\n";
			} elsif (exists $cross_matrix->{$genename}->{$sseqid}) {
				delete $cross_matrix->{$genename}->{$sseqid};
				$hit_string = "$sseqid, " . join (", ", keys %{$cross_matrix->{$genename}});
				print "best hit for $genename was $sseqid, this is a match\n";
			} else {
				print "$sseqid is not a match for the bait $genename\n";
			}
			if ($hit_string ne "") {
				print RESULTS_FH "$genename\t$evalue\t$length\t$slen\t$pident\t$hit_string\n";
			}
		}
		close FH;
		`rm $blastfile`;
	}
}

printlog ("Finished $runline");

close RESULTS_FH;


__END__

=head1 NAME

ValidateGenes.pl

=head1 SYNOPSIS

ValidateGenes.pl -input inputdir [-reference reffile | -db blastdb] -output outputdir

Given a fasta file of CDS sequences that were used as bait for aTRAM,
validates that the aTRAM results are in fact good matches for the target gene used.

=head1 OPTIONS

 -input:      directory containing aTRAM results to be validated against the original bait set.
 -reference:  FASTA file of the gene bait sequences used to aTRAM the input files.
 -processes:  optional: number of processes to be used.
 -output:     output directory (will be created if it doesn't exist).
 -database:   optional: if ValidateGenes.pl has been run before on this reffile, this is the name of the blastdb that was created from that run.

=cut

