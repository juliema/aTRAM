#!/usr/bin/env perl
package Mapreduce;
use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(set_multiplier get_multiplier map_to_shard set_total_shards get_total_shards get_max_shard);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

our $total_shards = 0;
our @primes = (29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151);
our $multiplier = 0;

#### MAPREDUCE FUNCTIONS

sub set_multiplier {
	my $factor = shift;
	$multiplier = $primes[$factor % (@primes)];
	return $multiplier;
}

sub get_multiplier {
	if ($multiplier == 0) {
		$multiplier = $primes[0];
	}
	return $multiplier;
}

sub map_to_shard {
	my $name = shift;

	$name =~ s/\/\d//;
	$name =~ s/#.+$//;
	$name =~ tr/0-9//dc;

	$name =~ /.*(\d{8})$/;
	$name = $1 * get_multiplier();

	return $name % $total_shards;
}

sub set_total_shards {
	$total_shards = shift;
	return $total_shards;
}

sub get_total_shards {
	my $dbname = shift;

	if ($total_shards == 0) {
		my $num = 0;
		while (-e "$dbname.$num.1.fasta") {
			$num++;
		}
		$total_shards = $num;
	}
	return $total_shards;
}

sub get_max_shard {
	return get_total_shards() - 1;
}

return 1;
