use strict;

our $debug = 0;

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
	my $log_fh = shift;

	print $log_fh ("\t$cmd\n");
    my $child_pid = fork();
    unless ($child_pid) { #child process
		exec ($cmd);
    } else { #parent process
        return $child_pid;
    }
}

sub wait_for_forks {
    while (@{@_[0]} > 0) {
    	my $item = pop @_[0];
        waitpid $item, 0;
    }
    return;
}

sub system_call {
	my $cmd = shift;
	my $log_fh = shift;

	unless ($log_fh) {
		$log_fh = &STDOUT;
	}

	print $log_fh ("\t$cmd\n");
	my ($saveout, $saveerr);
	if ($debug == 0) {
		open $saveout, ">&STDOUT";
		open $saveerr, ">&STDERR";
		open STDOUT, '>', File::Spec->devnull();
		open STDERR, '>', File::Spec->devnull();
	}
	my $exit_val = eval {
		system ($cmd);
	};

	if ($debug == 0) {
		open STDOUT, ">&", $saveout;
		open STDERR, ">&", $saveerr;
	}

	if ($exit_val != 0) {
		print "System call \"$cmd\" exited with $exit_val\n";
		exit;
	}
	return $exit_val;
}

sub debug {
	my $msg = shift;
	if ($debug) {
		print "$msg";
	}
}

my set_debug {
	my $debug_new = shift;
	$debug = $debug_new;
}

sub sortfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my (undef, $tempfile) = tempfile(UNLINK => 1);
	system ("gawk '{if (NF==0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | sort -t',' -k 1 | gawk '{print \">\" \$1 \"$separator\" \$2}' FS=\",\" > $outfile");
}

sub make_hit_matrix {
	my $blast_file = shift;
	my $hit_matrix = shift;

	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		my ($contig, $baitseq, $score, $qstart, $qend, $sstart, $send, $qlen, undef) = split(/\s+/,$line);
		my $strand = 1;
		if ($qend > $qstart) {
			$strand = -1;
		}
		my $currscore = $hit_matrix->{$contig}->{$baitseq};
		if ($score =~ /(\d+\.\d\d)/) {
			$score = $1;
		}
		$hit_matrix->{$contig}->{"length"} = $qlen;
		if ($currscore == undef) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
		} else {
			if (abs($currscore) < $score) {
			$hit_matrix->{$contig}->{$baitseq} = $strand * $score;
			}
		}
	}
	close BLAST_FH;

}

sub count_partial_libraries {
	my $libname = shift;

	my $num = 0;
	while (-e "$libname.$num.1.fasta") {
		$num++;
	}
	return $num;
}

return 1;
