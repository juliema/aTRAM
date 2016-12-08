#!/usr/bin/env perl
package Trinity;
use strict;
use File::Temp qw/ tempfile /;
use System;
use Parsing;
use Configuration;

# Assembler modules need to know:
	# where to find the short reads (pass this in as a file name)
	# what the assembly parameters are. (pass this in as a hash)
# Assembler modules should return a hash of the resulting contigs.

# Hash of assembler's required binaries
my $binaries = {trinity => "Trinity.pl"};

sub get_binaries {
	return $binaries;
}

sub assembler {
	my $self = shift;
	my $short_read_file = shift;
	my $params = shift;

	my $jm = "1G";

	Configuration::initialize();

	my ($kmer, $tempdir, $longreads, $ins_length, $exp_cov, $min_contig_len, $output_file) = 0;
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
		if (exists $params->{"output"}) {
			$output_file = $params->{"output"};
		}
	}
	# using Trinity.pl
	run_command (get_bin($binaries->{trinity}), "--seqType fa --single $short_read_file --run_as_paired --JM $jm --output $tempdir");

	my ($contigs, undef) = parsefasta ("$tempdir/Trinity.fasta");
	#### removing the trinity folder so trinity will now do more than one iteration
	`rm -r $tempdir/`;

	open OUTFH, ">", $output_file;
	foreach my $contigname (keys %$contigs) {
		my $sequence = $contigs->{$contigname};
		# >comp2_c0_seq1 len=4637 path=[7895:0-4636]
		$contigname =~ s/^comp(\d+)_(c\d+)_seq. len=(\d+).*$/$1$2_len_$3/;
		print OUTFH ">$contigname\n$sequence\n";
	}
	close OUTFH;
	return $contigs;
}

return 1;
