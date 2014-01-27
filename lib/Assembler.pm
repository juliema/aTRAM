#!/usr/bin/env perl
use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";

our %assemblers = {};

package Assembler;

sub parse_config {
	my $self = shift;
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
	my $self = shift;
	my $cmd = shift;

	if (exists $assemblers{$cmd}) {
		return "$assemblers{$cmd}";
	}
	return "$cmd";
}

return 1;
