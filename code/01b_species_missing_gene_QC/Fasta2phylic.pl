#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy;

my %opts;
getopts('i:o:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i Genome_folder -o output_folder\n";
my $output_folder = $opts{'o'} or die "use: $0 -i Genome_folder -o output_folder\n";

chdir "./$input_folder" or die "Cannot enter folder $input_folder: $!\n";
my @all_seq_files = glob "*.fas";
chdir "../";

mkdir "./$output_folder" or die "Cannot make folder $output_folder: $!\n";
foreach my $seq (@all_seq_files){
	my %sequ = &read_fasta("$input_folder/$seq");
	my $seq_name = $seq;
	$seq_name=~ s/\.fas/.phy/;
	open OUT, ">$output_folder/$seq_name" or die $!;
	my @seqs = keys %sequ;
	my $seq_one = $sequ{$seqs[0]};
	my $seq_length = length($seq_one);
	my $seq_number = scalar (@seqs);
	print OUT "$seq_number $seq_length CG\n";
	foreach my $i (sort keys %sequ){
		$sequ{$i} =~ s/!/-/g;
		print OUT "$i\n$sequ{$i}\n";
	}
	close OUT;
}


sub read_fasta{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	#open SEQ_INFO, ">>Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
	my %sequences;
	my $id;
	my $count = 0;
	my $seq_length;
	while (<READ_SEQ>){
	chomp;
	s/\s+//;
	if ($_ =~ /^>(?<seq_name>.+$)/){
		$id  = $+{seq_name};
		$id =~ s/carlito_syrichta/tarsius_syrichta/;
		$count += 1;
	}else{
		$sequences{$id} .= $_;
		$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences.\n";
	#print SEQ_INFO "Successfuly read @_, in total $count sequences.\n";
	#close SEQ_INFO;
	#print "\n";
	%sequences;
}
