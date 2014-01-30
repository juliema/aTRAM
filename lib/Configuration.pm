#!/usr/bin/env perl
package Configuration;
use strict;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw( initialize find_bin init_module );
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw( $binaries);
}

our $binaries = {};

sub initialize {
	open FH, "<", "$FindBin::Bin/config.txt" or die "Couldn't find $FindBin::Bin/config.txt. Did you run $FindBin::Bin/configure.pl?";
	foreach my $line (<FH>) {
		$line =~ s/(#.*)$//;
		if ($line =~ /(.*)=(.*)$/) {
			my $name = $1;
			my $path = $2;
			$binaries->{$name} = "$path";
		}
	}
}

sub find_bin {
	my $bin = shift;

	if (exists $binaries->{$bin}) {
		return "$binaries->{$bin}";
	}
	return "";
}

sub init_module {
	my $mod_bins = shift;
	foreach my $b (keys %$mod_bins) {
		$mod_bins->{$b} = find_bin($mod_bins->{$b});
		if ($mod_bins->{$b} eq "") {
			return 0;
		}
	}
	return 1;
}

return 1;
