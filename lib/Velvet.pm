#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile /;
require Subfunctions;

	# Assembler modules need to know:
	 	# where to find the short reads (pass this in as a file name)
	 	# what the assembly parameters are. (pass this in as a hash)
	# Assembler modules should return a file name for the resulting contigs.


sub assembler {
	my $short_read_file = shift;
	my $params = shift;
	my $log_fh = shift;

	my ($kmer, $tempdir, $longreads, $ins_length, $exp_cov, $min_contig_len) = 0;
	if (ref $params =~ /HASH/) {
		if ($params->{"kmer"}) {
			$kmer = $params->{"kmer"};
		}
		if ($params->{"tempdir"}) {
			$tempdir = $params->{"tempdir"};
		}
		if ($params->{"longreads"}) {
			$longreads = $params->{"longreads"};
		}
		if ($params->{"ins_length"}) {
			$ins_length = $params->{"ins_length"};
		}
		if ($params->{"exp_cov"}) {
			$exp_cov = $params->{"exp_cov"};
		}
		if ($params->{"min_contig_len"}) {
			$min_contig_len = $params->{"min_contig_len"};
		}
	}

	# using velvet
	print "\tassembling with Velvet...\n";
	if ($longreads != 0) {
		system_call ("velveth $tempdir $kmer -fasta -shortPaired $short_read_file -long $longreads", $log_fh);
	} else {
		system_call ("velveth $tempdir $kmer -fasta -shortPaired $short_read_file", $log_fh);
	}
	system_call ("velvetg $tempdir -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth $min_contig_len", $log_fh);

	return "$tempdir/contigs.fa";
}

return 1;
