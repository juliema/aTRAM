#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

my $debug = 0;
if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

#required: input file names
my $short_read_archive = 0;
my $search_fasta = 0;

#optional parameters:
my $help = 0;
my $log_file = 0;
my $use_ends = 0;
my $protein = 0;
my $blast_file = 0;

#parameters with modifiable default values
my $output_file = 0;
my $ins_length = 300;
my $iterations = 5;
my $start_iter = 1;
my $exp_cov = 30;

GetOptions ('reads=s' => \$short_read_archive,
            'target=s' => \$search_fasta,
            'insert_length|ins_length=i' => \$ins_length,
            'exp_coverage|expected_coverage=i' => \$exp_cov,
            'iterations=i' => \$iterations,
            'start_iteration=i' => \$start_iter,
            'log_file=s' => \$log_file,
            'use_ends' => \$use_ends,
            'output=s' => \$output_file,
            'protein' => \$protein,
            'blast=s' => \$blast_file,
            'debug' => \$debug,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive and $search_fasta) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with 0-prepare_files.pl) and a target gene in fasta form.");
}

# check to make sure that the specified short read archive exists:
unless ((-e "$short_read_archive.1.fasta") && (-e "$short_read_archive.2.fasta")) {
    pod2usage(-msg => "Short read archive does not seem to be in the format made by 0-prepare_files.pl. Did you specify the name correctly?");
}

print $runline;

my $executing_path = dirname(__FILE__);
my $cmd;

my $log_fh;
if ($log_file) {
	open $log_fh, ">", $log_file or die "couldn't open $log_file\n";
} else {
	open $log_fh, ">", "$output_file.log" or die "couldn't open $output_file.log\n";
}

print $log_fh $runline;

unless ($output_file) {
    $output_file = $short_read_archive;
}

unless ($blast_file) {
	(undef, $blast_file) = tempfile(UNLINK => 1);
}

my (undef, $targetdb) = tempfile(UNLINK => 1);
my (undef, $sort_file) = tempfile(UNLINK => 1);
my ($TARGET_FH, $target_fasta) = tempfile();

my $start_seq = "";
my $end_seq = "";
my $hit_matrix = ();

# process the target sequence file to look for the start seq and end seq.
my @target_seqs = ();
open SEARCH_FH, "<", $search_fasta;
my $name = "";
my $seq = "";
while (my $line=readline SEARCH_FH) {
	if ($line =~ />(.*)/) {
		if ($name ne "") {
			push @target_seqs, "$name,$seq";
		}
		$name = $1;
		chomp $name;
		$seq = "";
	} else {
		$seq .= $line;
		chomp $seq;
	}
}
push @target_seqs, "$name,$seq";

close SEARCH_FH;

my $len = 100;
if ($protein ==1) {
	$len = 30;
}

# check first target seq, see if we need to make just the first bit separated for checking completeness.
my $firstseq = shift @target_seqs;
if ($firstseq =~ /(.*),(.{$len})(.*)/) {
	# split into two smaller seqs:
	$start_seq = "sTRAM_target_start";
	unshift @target_seqs, "$1,$3";
	unshift @target_seqs, "$start_seq,$2";
} else {
	# it's short enough to just go as-is
	$firstseq =~ /(.*),/;
	$start_seq = $1;
	unshift @target_seqs, $firstseq;
}

# check last target seq, see if we need to make just the first bit separated for checking completeness.
my $lastseq = pop @target_seqs;
if ($lastseq =~ /(.*),(.*)(.{$len})/) {
	# split into two smaller seqs:
	$end_seq = "sTRAM_target_end";
	push @target_seqs, "$1,$2";
	push @target_seqs, "$end_seq,$3";
} else {
	# it's short enough to just go as-is
	$lastseq =~ /(.*),/;
	$end_seq = $1;
	push @target_seqs, $lastseq;
}

# okay, let's re-assemble the file for the target fasta.
foreach my $line (@target_seqs) {
	$line =~ /(.*),(.*)/;
	print $TARGET_FH ">$1\n$2\n";
}
close $TARGET_FH;

# make a database from the target so that we can compare contigs to the target.
if ($protein == 1) {
	$cmd = "makeblastdb -in $target_fasta -dbtype prot -out $targetdb.db -input_type fasta";
} else {
	$cmd = "makeblastdb -in $target_fasta -dbtype nucl -out $targetdb.db -input_type fasta";
}
system_call ($cmd);

if ($start_iter > 1) {
	my $x = $start_iter-1;
	$search_fasta = "$output_file.$x.contigs.fa";
}

for (my $i=$start_iter; $i<=$iterations; $i++) {
	print ("interation $i starting...\n");
	print $log_fh ("interation $i starting...\n");

	# 1. blast to find any short reads that match the target.
	print "\tblasting short reads...\n";
	if (($protein == 1) && ($i == 1)) {
		system_call ("tblastn -max_target_seqs 100000000 -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i");
	} else {
		system_call ("blastn -task blastn -evalue 10e-10 -max_target_seqs 100000000 -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i");
	}

	# did we not find any reads? Go ahead and quit.
	if ((-s "$output_file.blast.$i") == 0) {
		die "No similar reads were found.\n";
	}

	# 2 and 3. find the paired end of all of the blast hits.
	print "\tgetting paired ends...\n";
	system_call("perl $executing_path/lib/pairedsequenceretrieval.pl $short_read_archive.#.fasta $output_file.blast.$i $output_file.blast.$i.fasta");

	# 4. assemble the reads using velvet
	print "\tassembling with Velvet...\n";
	system_call("velveth $output_file.velvet 31 -fasta -shortPaired $output_file.blast.$i.fasta");
	system_call("velvetg $output_file.velvet -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth 200");

	# we'll use the resulting contigs as the query for the next iteration.
	$search_fasta = "$output_file.$i.contigs.fa";
	system_call("cp $output_file.velvet/contigs.fa $search_fasta");
	print "starting with this: \n";
	system_call("grep '>' $search_fasta");

	# 5. now we filter out the contigs to look for just the best ones.
	print "\tfiltering contigs...\n";
	if ($protein == 1) {
		system_call ("blastx -db $targetdb.db -query $search_fasta -out $blast_file -outfmt '6 qseqid sseqid bitscore'");
	} else {
		system_call ("tblastx -db $targetdb.db -query $search_fasta -out $blast_file -outfmt '6 qseqid sseqid bitscore'");
	}

	# 6. we want to keep the contigs that have a bitscore higher than 100.
# 	system_call("gawk '{if (\$2 ~ \/$start_seq\/) print \$0\"\\tstart\"; else print \$0;}' $blast_file > $sort_file");
# 	system_call("gawk '{if (\$2 ~ \/$end_seq\/) print \$0\"\\tend\"; else print \$0;}' $sort_file > $blast_file");
# 	system_call("gawk '{print \$3\"\\t\"\$1;}' $blast_file | sort -n -r | gawk '{if (\$1 > 100) print \$2;}' | sort | uniq > $sort_file");
	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		$line =~ /(.*?)\s+(.*?)\s+(.*)\s/;
		my $contig = $1;
		my $baitseq = $2;
		my $score = $3;
		my $currscore = $hit_matrix->{$contig}->{$baitseq};
		if ($currscore == undef) {
			$hit_matrix->{$contig}->{$baitseq} = $score;
		} else {
			if ($currscore < $score) {
				$hit_matrix->{$contig}->{$baitseq} = $score;
			}
		}
	}
	close BLAST_FH;

	foreach my $contig (keys $hit_matrix) {
		my $contig_high_score = 0;
		foreach my $baitseq (keys $hit_matrix->{$contig}) {
			my $score = $hit_matrix->{$contig}->{$baitseq};
			if ($score > $contig_high_score) {
				$contig_high_score = $score;
			}
		}
		if ($contig_high_score < 100) {
			delete $hit_matrix->{$contig};
		}
	}

	open SORT_FH, ">", $blast_file;
	foreach my $contig (keys $hit_matrix) {
		print SORT_FH $contig . "\n";
	}
	system_call("sort -n $blast_file > $sort_file");
	system_call("perl $executing_path/lib/findsequences.pl $search_fasta $sort_file $search_fasta");
	print "ending with this: \n";
	system_call("grep '>' $search_fasta");
	# save off these resulting contigs to the ongoing contigs file.
	open CONTIGS_FH, ">>", "$output_file.all.fasta";
	print CONTIGS_FH `cat $search_fasta | gawk '{sub(/>/,">$i\_"); print \$0}'`;
	close CONTIGS_FH;

	# if we flagged to use just the ends of the contigs, clean that up.
	if ($use_ends != 0) {
		system_call("bash $executing_path/lib/sort_contigs.sh $search_fasta");

		open FH, "<", "$search_fasta.sorted.tab";
		my @contigs = <FH>;
		close FH;
		my @new_contigs = ();

		for (my $j=0; $j<@contigs; $j++) {
			my ($len,$name,$seq) = split (/\t/,$contigs[$j]);
			if ($len > ($ins_length * 5)) {
				my $name_start = $name. "_start";
				my $name_end = $name . "_end";
				$seq =~ /^(.{$ins_length}).*(.{$ins_length})$/;
				my $seq_start = $1;
				my $seq_end = $2;
				push @new_contigs, "$ins_length\t$name_start\t$seq_start";
				push @new_contigs, "$ins_length\t$name_end\t$seq_end";
			} else {
				push @new_contigs, $contigs[$j];
			}
		}

		open FH, ">", "$search_fasta";
		foreach my $line (@new_contigs) {
			chomp $line;
			$line =~ /(.*?)\t(.*?)\t(.*)/;
			print FH ">$2\n$3\n";
		}
		close FH;
	}
}

close $log_fh;

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

sub system_call {
	my $cmd = shift;
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

__END__

=head1 NAME

sTRAM.pl

=head1 SYNOPSIS

sTRAM.pl -reads shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

=head1 OPTIONS

  -reads:     		short read archive (already run through makelibrary.pl).
  -target:          fasta file with sequences of interest.
  -output:	        optional: the prefix for the pipeline's output files (default name is the same as -reads).
  -ins_length:	    optional: the size of the fragments used in the short-read library (default 300).
  -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
  -iterations:      optional: the number of pipeline iterations (default 5).
  -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
  -log_file:        optional: a file to store output of the pipeline.
  -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.


=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

