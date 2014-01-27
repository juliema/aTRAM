#!/usr/bin/env perl
use strict;
use File::Temp qw/ tempfile /;
use Module::Load;
load Assembler;
use Subfunctions;

# Assembler modules need to know:
	# where to find the short reads (pass this in as a file name)
	# what the assembly parameters are. (pass this in as a hash)
# Assembler modules should return a file name for the resulting contigs.

package Trinity;

sub assembler {
	my $self = shift;
	my $short_read_file = shift;
	my $params = shift;

	my $jm = "1G";

	my $path = Assembler->find_bin("Trinity.pl");

	my ($kmer, $tempdir, $longreads, $ins_length, $exp_cov, $min_contig_len) = 0;
	if ((ref $params) =~ /HASH/) {
        if (exists $params->{"jm"}) {
			$jm = $params->{"jm"};
		}
		if (exists $params->{"tempdir"}) {
			$tempdir = $params->{"tempdir"};
		}
		if (exists $params->{"longreads"}) {
			$longreads = $params->{"longreads"};
		}
	}
	# using Trinity.pl
	print "\tassembling with Trinity...\n";
# perl ~/packages/trinityrnaseq_r20131110/Trinity.pl --seqType fa --single Pop_delt_psbA_atpA.1.blast.fasta --run_as_paired --JM 10G
	Subfunctions::system_call ("$path --seqType fa --single $short_read_file --run_as_paired --JM $jm --output $tempdir");

	my ($contigs, $contignames) = Subfunctions::parsefasta ("$tempdir/Trinity.fasta");
	open FH, ">", "$tempdir/results.fasta";
	foreach my $contig (@$contignames) {
		print "$contig\n";
		print FH ">$contig\n$contigs->{$contig}\n";
	}
	close FH;
	return $contigs;
}

sub write_contig_file {
	my $self = shift;
	my $contigs = shift;
	my $renamefile = shift;
	my $prefix = shift;

	if ($prefix) {
		$prefix = "$prefix.";
	} else {
		$prefix = "";
	}

	open OUTFH, ">", $renamefile;
	foreach my $contigname (keys $contigs) {
		my $sequence = $contigs->{$contigname};
		# >comp2_c0_seq1 len=4637 path=[7895:0-4636]
		$contigname =~ s/^comp(\d+)_(c\d+)_seq. len=(\d+).*$/$prefix$1$2_len_$3/;
		print OUTFH ">$contigname\n$sequence\n";
	}
	close OUTFH;
}

return 1;
