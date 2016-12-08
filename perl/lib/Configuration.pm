#!/usr/bin/env perl
package Configuration;
use strict;
use File::Spec;
use File::Find;
use File::Basename;
use Module::Load;
use Cwd qw(realpath);
use System;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(get_req_software get_optional_software check_module get_assemblers get_bin get_config_file get_atram_path);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

my $binaries = {};
my @req_software = qw(blastn tblastn blastx tblastx makeblastdb);
my @opt_software = qw(muscle mafft);
my $assemblers = {};
my $assembler_dir = "";
my $config_file = "";
my $atrampath = "";

sub initialize {
	if (%$binaries) {
		# we've already initialized.
		return;
	}
	$config_file = get_config_file();
	$atrampath = get_atram_path();
	$assembler_dir = File::Spec->catdir($atrampath, "lib", "Assembler");
	if ((-s $config_file) > 0) {
		open my $fh, "<:crlf", $config_file or die "Couldn't find $config_file. Did you run configure.pl?";
		foreach my $line (<$fh>) {
			$line =~ s/(#.*)$//;
			if ($line =~ /(.+)=(.+)$/) {
				my $name = $1;
				my $path = $2;
				$binaries->{$name} = "$path";
			}
		}
	}
}

sub get_config_file {
	if ($atrampath eq "") {
		$atrampath = get_atram_path();
	}
	$config_file = File::Spec->catfile($atrampath, "config.txt");

	# perhaps it is in the old location in lib: move it to here, then carry on.
	if (!(-e $config_file)) {
		my $oldconfigfile = File::Spec->catfile($atrampath, "lib", "config.txt");
		if (-e $oldconfigfile) {
			run_command ("mv", "$oldconfigfile $config_file");
		}
	}

	return $config_file;
}

sub get_atram_path {
	# find this module's path in %INC:
	my $modpath = realpath($INC{'Configuration.pm'});

	# we know that this module is in the lib directory, which is one dir inside the main aTRAM dir.
	my @pathpieces = File::Spec->splitdir($modpath);
	pop @pathpieces; # this is the file Configuration.pm
	pop @pathpieces; # this is lib
	$atrampath = File::Spec->catdir(@pathpieces);
	return $atrampath;
}

sub get_req_software {
	return \@req_software;
}

sub get_optional_software {
	return \@opt_software;
}

sub get_assemblers {
	if (%$assemblers) {
		return $assemblers;
	}

	find ( {wanted => \&assembler_bins, no_chdir => 1} , "$assembler_dir");
	return $assemblers;
}

sub get_bin {
	my $bin = shift;
	if (exists $binaries->{$bin}) {
		return "$binaries->{$bin}";
	} else {
		printlog ("The binary \"$bin\" was not found in the config.txt file.");
		return "";
	}
}

sub find_bin {
	my $bin = shift;

	if (exists $binaries->{$bin}) {
		return "$binaries->{$bin}";
	} else {
		# if we don't have a path for $sw, ask the system.
		my $min_bin = basename($bin);
		my @pathlist = File::Spec->path();

		foreach my $path (@pathlist) {
			my $cmdpath = File::Spec->catpath("", $path, $min_bin);
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

sub check_module {
	my $mod = shift;
	my $mod_bins = $assemblers->{$mod};
	foreach my $b (@$mod_bins) {
		my $bin = find_bin($b);
		if ($bin eq "") {
			return 0;
		}
	}
	return 1;
}

return 1;
