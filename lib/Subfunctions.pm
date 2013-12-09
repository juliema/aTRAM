use strict;
use File::Temp qw/ tempfile /;
use FindBin;
use lib "$FindBin::Bin/lib";

our $debug = 0;
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
		print "found $cmd: at $assemblers{$cmd}\n";
	}
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

sub set_debug {
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

sub flattenfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my (undef, $tempfile) = tempfile(UNLINK => 1);
	system ("gawk '{if (NF==0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | gawk '{print \">\" \$1 \"$separator\" \$2}' FS=\",\" > $outfile");
}

sub make_hit_matrix {
	my $blast_file = shift;
	my $hit_matrix = shift;

	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		my ($contig, $baitseq, $score, $qstart, $qend, $sstart, $send, $qlen, undef) = split(/\s+/,$line);
		my $strand = 1;
		if ((($qend-$qstart) / ($send-$sstart)) < 0) {
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

sub process_hit_matrix {
	my $raw_hit_matrix = shift;
	my $targets = shift;
	my $bitscore = shift;
	my $contiglength = shift;
	my $hit_matrix = shift;
	my $contig_names = shift;

	my $high_score = 0;

	# clean up the hit matrix: only keep hits that meet the bitscore threshold.

	foreach my $contig (keys $raw_hit_matrix) {
		my $contig_high_score = 0;
		my $total = 0;
		$raw_hit_matrix->{$contig}->{"strand"} = 1;
		foreach my $baitseq (@$targets) {
			my $partscore = abs($raw_hit_matrix->{$contig}->{$baitseq});
			if ($partscore) {
				my $partstrand = ($raw_hit_matrix->{$contig}->{$baitseq})/$partscore;
				if ($partscore > 0) {
					# separate out the score and the strand for this part:
					$raw_hit_matrix->{$contig}->{$baitseq} = $partscore;
					$total += $partscore;
					if ($partscore > $contig_high_score) {
						# if this is the best score for the contig, set the contig_high_score and set the strand to this strand.
						$contig_high_score = $partscore;
						$raw_hit_matrix->{$contig}->{"strand"} = $partstrand;
					}
				}
			}
		}
		$raw_hit_matrix->{$contig}->{"total"} = $total;
		if ($total > $high_score) {
			$high_score = $total;
		}
		if (($total >= $bitscore) && ($raw_hit_matrix->{$contig}->{"length"} >= $contiglength)) {
			$hit_matrix->{$contig} = $raw_hit_matrix->{$contig};
			push @$contig_names, $contig;
		}
	}
	return $high_score;
}

sub count_partial_libraries {
	my $libname = shift;

	my $num = 0;
	while (-e "$libname.$num.1.fasta") {
		$num++;
	}
	return $num;
}

sub is_protein {
	my $sequence = shift;
	if ($sequence =~ /[EFILPQ]/) {
		return 1;
	} else {
		return 0;
	}
}

sub percentcoverage {
	my $reffile = shift;
	my $contigfile = shift;
	my $gene = shift;

	###### cat files
	my (undef, $catfile) = tempfile(UNLINK => 1);
	system ("cat $reffile $contigfile > $catfile");
	##### muscle alignment
	system ("muscle -in $catfile -out $gene.muscle.fasta");

	my (undef, $fastafile) = tempfile(UNLINK => 1);
	flattenfasta("$gene.muscle.fasta", $fastafile, ",");

	# parse the output file: save the reference as a separate sequence, put the others into an array.
	my $refseq = "";
	my $contigs = {};
	open FH, "<", $fastafile;
	while (my $line = readline FH) {
		$line =~ />(.*?),(.*)$/;
		my $name = $1;
		my $seq = $2;
		if ($name =~ /$gene/) {
			$refseq = $seq;
		} else {
			$contigs->{$name} = $seq;
		}
	}
	close FH1;

	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($refseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys $contigs) {
			$contigs->{$contig} =~ /(.{$leftlength})(.{$gaplength})(.*)/;
			$contigs->{$contig} = "$1$3";
		}

		$refseq = "$left$remainder";
	}

	####### Print out EXON file
	####### Print OUT Table

	open TABLE_FH, ">", "$gene.Table.txt";
	open EXON_FH, ">", "$gene.exons.fasta";

	print EXON_FH ">$gene\n$refseq\n";
	print TABLE_FH "contig\ttotal\tpercent\n";
	my $total_length = length $refseq;
	foreach my $contig (keys $contigs) {
		print EXON_FH ">$contig\n$contigs->{$contig}\n";
		my $gaps = ($contigs->{$contig} =~ tr/N-//);
		my $total = $total_length - $gaps;
		my $percent = $total / $total_length;
		print TABLE_FH "$contig\t$total\t$percent\n";
	}

	close TABLE_FH;
	close EXON_FH;
}


return 1;
