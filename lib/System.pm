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
	our @EXPORT      = qw(timestamp exit_with_msg fork_cmd wait_for_forks printlog system_call debug set_debug set_log);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

our $debug = 0;
our $log_fh = 0;
our $log_file = "";

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
	print $log_fh ("\t$cmd\n");
    my $child_pid = fork();
    unless ($child_pid) { #child process
		exec ($cmd);
    } else { #parent process
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
        }
    }
    return;
}

sub system_call {
	my $cmd = shift;
	my $log_me = shift;

	open my $saveout, ">&STDOUT";
	open my $saveerr, ">&STDERR";

	if ($debug == 0) {
		# if we're not debugging, dump the system's output to devnull.
		open STDOUT, '>', File::Spec->devnull();
		open STDERR, '>', File::Spec->devnull();
	}

	if (defined get_log_file()) {
		# if a log file has been specified, dump to that.
		printlog ("Running command \"$cmd\"");
		open STDOUT, ">>", get_log_file();
		open STDERR, ">>", get_log_file();
	}

	my $exit_val = eval {
		system ($cmd);
	};

	# unwind the redirects.
	if (defined get_log_file()) {
		close STDOUT;
		close STDERR;
		printlog ("Returning ".($exit_val >> 8)." from \"$cmd\"");
	}

	open STDOUT, ">&", $saveout;
	open STDERR, ">&", $saveerr;

	if ($? == 2) {
		# user signaled kill, so we should die.
		print "System call \"$cmd\" exited with $exit_val\n";
		exit $exit_val >> 8;
	}

	if (($exit_val != 0) && !(defined $log_me)) {
		print "System call \"$cmd\" exited with $exit_val\n";
		exit $exit_val >> 8;
	}
	return $exit_val >> 8;
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
	open my $temp_fh, ">>", $log_file or print "Couldn't open $log_file\n";
	$log_fh = $temp_fh;
	debug ("Setting log file to $log_new");
}

sub get_log_file {
	if ($log_file eq "") {
		return undef;
	}
	return $log_file;
}

sub printlog {
	my $msg = shift;

	$msg = timestamp() . ": " . $msg . "\n";
	if ($log_fh) {
        select($log_fh);
        $| = 1;
		print $log_fh $msg;
		select(STDOUT);
	} else {
		print $msg;
	}
}

return 1;
