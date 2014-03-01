#!/usr/bin/env perl
package Configuration;
use strict;
use File::Spec;
use File::Find;
use File::Basename;
use Module::Load;
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
our @req_software = qw(blastn tblastn blastx tblastx makeblastdb);
our $assemblers = {};
our $assembler_dir = "";
our $config_file = "";

sub initialize {
	# find the lib path in @INC:
	my $libpath = File::Spec->catfile("aTRAM", "lib");
	foreach my $path (@INC) {
		$path = realpath ($path);
		if ($path =~ /$libpath/) {
			$config_file = File::Spec->catfile($path, "config.txt");
			$assembler_dir = File::Spec->catdir($path, "Assembler");
			last;
		}
	}

	if ((-s $config_file) > 0) {
		open my $fh, "<", $config_file or die "Couldn't find $config_file. Did you run configure.pl?";
		foreach my $line (<$fh>) {
			$line =~ s/(#.*)$//;
			if ($line =~ /(.*)=(.*)$/) {
				my $name = $1;
				my $path = $2;
				$binaries->{$name} = "$path";
			}
		}
	}
}

sub get_req_software {
	return \@req_software;
}

sub get_assemblers {
	if (%$assemblers) {
		return $assemblers;
	}

	find ( {wanted => \&assembler_bins, no_chdir => 1} , "$assembler_dir");
	return $assemblers;
}

sub find_bin {
	my $bin = shift;

	if (exists $binaries->{$bin}) {
		return "$binaries->{$bin}";
	} else {
		# if we don't have a path for $sw, ask the system.
		my @pathlist = File::Spec->path();
		foreach my $path (@pathlist) {
			my $cmdpath = File::Spec->catpath("", $path, $bin);
			if (-x $cmdpath) {
				return "$cmdpath";
			}
		}
	}
	return "";
}

sub assembler_bins {
	if ( $File::Find::name ne $assembler_dir) {
		my $modname = basename ($_);
		$modname =~ s/\.pm//;
		load "Assembler::$modname";
		my $bins = $modname->get_binaries();
		my @binnames = values %$bins;
		$assemblers->{$modname} = \@binnames;
	}
	return;
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
