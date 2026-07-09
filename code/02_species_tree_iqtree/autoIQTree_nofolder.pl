#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
#version 3
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy qw(move);

my %opts;
getopts('o:', \%opts);
my $output_folder = $opts{'o'} or die "use: $0 -o output_fas_folder\n";



mkdir "./${output_folder}_contree" or die "Cannot create folder ${output_folder}_contree: $!\n";
mkdir "./${output_folder}_treefile" or die "Cannot create folder ${output_folder}_treefile: $!\n";
mkdir "./${output_folder}_best_scheme" or die "Cannot create folder ${output_folder}_best_scheme: $!\n";
mkdir "./failed" or die "Cannot create folder finished: $!\n";
#my $output_folder = "Tree_phy";
#chdir "$output_folder" or die "Cannot enter $output_folder: $!\n";
my @all_seq = glob "*.fas";
#chdir "../";

my %seq_length;
my $count = 0;
open LENGTH, "./seq_length.txt" or die "Cannot open seq_length.txt: $!\n";
while (<LENGTH>){
	chomp;
	my @lines = split (/\s+/, $_);
	$seq_length{$lines[0]} = $lines[1];
	$count ++;
}
say $count;



foreach my $i (@all_seq){
	open PART1, "./_partition.txt" or die "Cannot open _partition.txt: $!\n";
	open PART2, ">>./partition.txt" or die "Cannot create partition.txt: $!\n";
	my $seq = $i;
	$seq =~ s/\.fas//;
	while (<PART1>) {
		chomp;
		s/seq_length/$seq_length{$seq}/g;
		s/\r//;
		print PART2 "$_\n";
	}
	#iqtree -nt AUTO -s A1CF.fas -sp partition.txt -m MFP+MERGE -bb 1000
	system "iqtree2 -s $i -redo -sp partition.txt -m MFP+MERGE -B 1000 -T AUTO";
	#system "./raxmlHPC-PTHREADS-SSE3 -T 8 -# 100 -o sarcophilus_harrisii,monodelphis_domestica -n out -q partition.txt -M -f a -m GTRGAMMAI -x $random -p $random -s ./$output_folder/$i";
	$i =~ s/\.fas//;
	move "./partition.txt.contree", "./${output_folder}_contree/$i.contree.nwk";
	move "./partition.txt.treefile", "./${output_folder}_treefile/$i.treefile.nwk";
	move "./partition.txt.best_scheme.nex", "./${output_folder}_best_scheme/$i.best_scheme.nex";
	if (!(-f "./${output_folder}_best_scheme/$i.best_scheme.nex")){
		move "./$i.fas", "./failed/$i.fas";
	}
	
	system "rm partition.txt";
	#$random ++;
}


