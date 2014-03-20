#!/usr/bin/env perl
package Velvet;
use strict;
use System;
use Parsing;
use Configuration;

# Assembler modules need to know:
	# where to find the short reads (pass this in as a file name)
	# what the assembly parameters are. (pass this in as a hash)
# Assembler modules should return a hash of the resulting contigs.

# Hash of assembler's required binaries
my $binaries = {velveth => "velveth", velvetg => "velvetg"};

sub get_binaries {
	return $binaries;
}

sub assembler {
	my $self = shift;
	my $short_read_file = shift;
	my $params = shift;

	Configuration::initialize();

	my ($kmer, $tempdir, $ins_length, $exp_cov, $min_contig_len, $output_file) = 0;
	my $longreads = "";

	if ((ref $params) =~ /HASH/) {
        if (exists $params->{"kmer"}) {
			$kmer = $params->{"kmer"};
		}
		if (exists $params->{"tempdir"}) {
			$tempdir = $params->{"tempdir"};
		}
		if (exists $params->{"longreads"}) {
			$longreads = $params->{"longreads"};
		}
		if (exists $params->{"ins_length"}) {
			$ins_length = $params->{"ins_length"};
		}
		if (exists $params->{"exp_cov"}) {
			$exp_cov = $params->{"exp_cov"};
		}
		if (exists $params->{"min_contig_len"}) {
			$min_contig_len = $params->{"min_contig_len"};
		}
		if (exists $params->{"output"}) {
			$output_file = $params->{"output"};
		}
		if (exists $params->{"log_file"}) {
			set_log($params->{"log_file"});
		}
	}
	# using velvet

	# truncate Velvet log file if it already exists
	truncate "$tempdir/Log", 0;

	if ($longreads ne "") {
		run_command (get_bin($binaries->{velveth}), "$tempdir $kmer -fasta -shortPaired $short_read_file -long $longreads", 1);
	} else {
		run_command (get_bin($binaries->{velveth}), "$tempdir $kmer -fasta -shortPaired $short_read_file", 1);
	}
	run_command (get_bin($binaries->{velvetg}), "$tempdir -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth $min_contig_len", 1);
	my ($contigs, undef) = parsefasta ("$tempdir/contigs.fa");

	# copy Velvet log output to logfile.
	open LOGFH, "<:crlf", "$tempdir/Log";
	printlog ("Velvet log:");
	foreach my $line (<LOGFH>) {
		chomp $line;
		printlog ($line);
	}
	printlog ("end Velvet log");
	close LOGFH;

	open OUTFH, ">", $output_file;
	foreach my $contigname (keys %$contigs) {
		my $sequence = $contigs->{$contigname};
		#NODE_41_length_2668_cov_4.901050
		$contigname =~ s/^NODE_(\d+)_length_(\d+)_cov_(\d+\.\d).*$/$1_len_$2_cov_$3/;
		print OUTFH ">$contigname\n$sequence\n";
	}
	close OUTFH;
	return $contigs;
}

return 1;
