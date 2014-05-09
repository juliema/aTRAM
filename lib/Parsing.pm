#!/usr/bin/env perl
package Parsing;
use strict;
use File::Temp qw/ tempfile /;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(parsefasta sortfasta flattenfasta make_hit_matrix process_hit_matrix split_seq is_protein);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

sub parsefasta {
	my $fastafile = shift;

	my $taxa = {};
	my @taxanames = ();
	open fileIN, "<:crlf", "$fastafile";
	my $taxonlabel = "";
	my $sequence = "";
	while (my $input = readline fileIN) {
		if ($input =~ /^>(.+)\s*$/) {
			$taxonlabel = $1;
			$taxonlabel =~ s/\s+/_/g;
			push @taxanames, $taxonlabel;
		} else {
			$input =~ /^\s*(.+)\s*$/;
			$taxa->{$taxonlabel} .= $1;
		}
	}

	close fileIN;
	return $taxa, \@taxanames;
}

sub sortfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my ($taxa, $taxanames) = parsefasta ($fastafile);
	my @sortedtaxa = sort @$taxanames;
	open fileOUT, ">", "$outfile";
	foreach my $taxon (@sortedtaxa) {
		print fileOUT ">$taxon$separator$taxa->{$taxon}\n";
	}
	close fileOUT;
}

sub flattenfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my ($taxa, $taxanames) = parsefasta ($fastafile);
	open fileOUT, ">", "$outfile";
	foreach my $taxon (@$taxanames) {
		print fileOUT ">$taxon$separator$taxa->{$taxon}\n";
	}
	close fileOUT;
}

sub make_hit_matrix {
	my $blast_file = shift;
	my $hit_matrix = shift;
	my $iter = shift;

	open BLAST_FH, "<:crlf", $blast_file;
	while (my $line = readline BLAST_FH) {
		my ($contig, $baitseq, $score, $qstart, $qend, $sstart, $send, $qlen, undef) = split(/\s+/,$line);
		if (!(defined $contig)) {
			return undef;
		}
		my $strand = 1;
		if ((($qend-$qstart) / ($send-$sstart)) < 0) {
			$strand = -1;
		}
		my $currscore = $hit_matrix->{$contig}->{$baitseq};
		if ($score =~ /(\d+\.\d\d)/) {
			$score = $1;
		}
		$hit_matrix->{$contig}->{"length"} = $qlen;
		$hit_matrix->{$contig}->{"iteration"} = $iter;
		if (!(defined $currscore)) {
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

	foreach my $contig (keys %$raw_hit_matrix) {
		my $contig_high_score = 0;
		my $total = 0;
		$raw_hit_matrix->{$contig}->{"strand"} = 1;
		foreach my $baitseq (@$targets) {
			my $partscore = $raw_hit_matrix->{$contig}->{$baitseq};
			if (defined $partscore) {
				$partscore = abs($partscore);
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
			# check to see if this contig is already in the hit_matrix:
			if (exists $hit_matrix->{$contig}) {
				if (($hit_matrix->{$contig}->{"length"} == $raw_hit_matrix->{$contig}->{"length"}) && ($hit_matrix->{$contig}->{"total"} == $raw_hit_matrix->{$contig}->{"total"})) {
					# it's the same, don't add it.
					next;
				}
			}
			$hit_matrix->{$contig} = $raw_hit_matrix->{$contig};
			push @$contig_names, $contig;
		}
	}
	return $high_score;
}

sub split_seq {
    my $seq = shift;
    my $start = shift;
    my $end = shift;
    my $max = 30000;
        my $seqlen = length ($seq);
        my $startseq = "";
        my $regionseq = "";
        my $endseq = "";

        my $currstart = $start-1;
        my $currend = $end;
        while ($currstart > $max) {
                $seq =~ /^(.{$max})(.*)$/;
                $startseq .= $1;
                $seq = $2;
                $currstart -= $max;
                $currend -= $max;
        }
        if ($currstart > 0) {
                $seq =~ /^(.{$currstart})(.*)$/;
                $startseq .= $1;
                $seq = $2;
                $currstart = 1;
                $currend = $end - (length ($startseq));
        }

        my $regionsize = $end - $start + 1;
        while ($regionsize > $max) {
                $seq =~ /^(.{$max})(.*)$/;
                $regionseq .= $1;
                $seq = $2;
                $currstart -= $max;
                $currend -= $max;
                $regionsize -= $max;
        }
        if ($regionsize > 0) {
                $seq =~ /^(.{$regionsize})(.*)$/;
                $regionseq .= $1;
                $endseq = $2;
        }
        return ($startseq, $regionseq, $endseq);
}

sub is_protein {
	my $sequence = shift;
	if ($sequence =~ /[EFILPQ]/) {
		return 1;
	} else {
		return 0;
	}
}

return 1;
