#!/usr/bin/env perl
package Assembler;
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
	our @EXPORT      = qw( parse_config find_bin );
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

our %assemblers = {};


sub parse_config {
	open FH, "<", "$FindBin::Bin/config.txt";
	foreach my $line (<FH>) {
		$line =~ s/(#.*)$//;
		if ($line =~ /(.*)=(.*)$/) {
			my $name = $1;
			my $path = $2;
			$assemblers{$name} = "$path";
		}
	}
}

sub find_bin {
	my $cmd = shift;

	if (exists $assemblers{$cmd}) {
		return "$assemblers{$cmd}";
	}
	return "$cmd";
}

return 1;
