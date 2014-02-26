#!/usr/bin/env perl
package Configuration;
use strict;
use File::Spec;
use Cwd qw(realpath);

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

our $binaries = {};

sub initialize {
	# find the lib path in @INC:
	my $config_file = "";
	my $libpath = File::Spec->catfile("aTRAM", "lib");
	foreach my $path (@INC) {
		$path = realpath ($path);
		print "looking at $path\n";
		if ($path =~ /$libpath/) {
			$config_file = File::Spec->catfile($path, "config.txt");
			last;
		}
	}

	my $fh;
	if (!(defined (open $fh, "<", $config_file))) {
		die "Couldn't find $config_file. Did you run configure.pl?";
	}
	foreach my $line (<$fh>) {
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
