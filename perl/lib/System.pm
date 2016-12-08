#!/usr/bin/env perl
package System;
use strict;
use File::Temp qw/ tempfile /;
use Parsing;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(timestamp exit_with_msg fork_cmd wait_for_forks printlog run_command debug set_debug set_log get_log_file close_log get_version);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

my $debug = 0;
my $log_fh = 0;
my $log_file = "";
my $version = "v1.01+";

sub get_version {
	return $version;
}

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return "$date $time";
}

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

sub fork_cmd {
	my $cmd = shift;
	my $args = shift;

	# check to make sure the command even exists.
	if ($cmd eq "") {
		printlog ("Can't fork: no command given.");
		die;
	}

    my $child_pid = fork();
    unless ($child_pid) { #child process
		exec ("$cmd $args");
    } else { #parent process
		printlog ("Forking pid $child_pid: $cmd");
        return $child_pid;
    }
}

sub wait_for_forks {
    while (@{$_[0]} > 0) {
     	debug ("waiting for " . join(", ", @{$_[0]}) );
		my $item = pop @{$_[0]};
        waitpid $item, 0;
        if ($? != 0) {
        	printlog ("child process $item died with error $?");
        	@{$_} = ();
        	die 3;
        }
    }
    return 0;
}

sub run_command {
	my $cmd = shift;
	my $args = shift;
	my $params = shift;

	my $no_exit = 0;
	# handle params, if provided.
	if ((ref $params) =~ /HASH/) {
		$no_exit = $params->{"no_exit"};
	}

	# check to make sure the command even exists.
	if ($cmd eq "") {
		printlog ("No command given.");
		die "No command given.";
	}

	open my $saveout, ">&STDOUT";
	open my $saveerr, ">&STDERR";

	if ($debug == 0) {
		# if we're not debugging, dump the system's output to devnull.
		open STDOUT, '>', File::Spec->devnull();
		open STDERR, '>', File::Spec->devnull();
	} elsif (defined get_log_file()) {
		# if a log file has been specified, dump to that.
		printlog ("Running command \"$cmd $args\"");
		open STDOUT, ">>", get_log_file();
		open STDERR, ">>", get_log_file();
	}

	my $result = system("$cmd $args");

	# unwind the redirects if we had pointed them at the log file.
	if (defined get_log_file()) {
		close STDOUT;
		close STDERR;
	}

	# make sure that we return STDOUT and STDERR to their original states.
	open STDOUT, ">&", $saveout;
	open STDERR, ">&", $saveerr;

	# Now we should unpack the exit value.
	my $exit_val = $result >> 8;
	my $signal = $? & 255;

	if (($exit_val != 0) && ($no_exit != 1)) {
		# if the command returned nonzero and the user didn't specify no_exit, print message and die.
		print "\nCommand \"$cmd\" exited with $exit_val\n";
		exit $exit_val;
	}

	if ($signal != 0) {
		# signal caught, so we should die.
		print "\nSignal $signal ($?) caught on command \"$cmd\"\n";
		exit $exit_val;
	}
	printlog ("Returning $exit_val from \"$cmd\"");
	return $exit_val;
}

sub debug {
	my $msg = shift;
	if ($debug) {
		printlog ("DEBUG: $msg");
	}
}

sub set_debug {
	my $debug_new = shift;
	$debug = $debug_new;
}

sub set_log {
	my $log_new = shift;
	$log_file = $log_new;
	if ($log_new ne "") {
		open my $temp_fh, ">>", $log_file or print "Couldn't open $log_file\n";
		$log_fh = $temp_fh;
		debug ("Setting log file to $log_new");
	} else {
		$log_file = "";
		$log_fh = 0;
	}
}

sub get_log_file {
	if ($log_file eq "") {
		return undef;
	}
	return $log_file;
}

sub close_log {
	if ($log_fh) {
		close $log_fh;
	}
}

sub printlog {
	my $msg = shift;
	my $echo = shift;

	if ($echo == 1) {
		my $oldpipe = $|;
		select(STDOUT);
        $| = 1;
		print $msg. "\n";
		$| = $oldpipe;
	}
	$msg = timestamp() . ": " . $msg . "\n";
	if ($log_fh) {
        select($log_fh);
        $| = 1;
		print $log_fh $msg;
		select(STDOUT);
	}
}

return 1;
