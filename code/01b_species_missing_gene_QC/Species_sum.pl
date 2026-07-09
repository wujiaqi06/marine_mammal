#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use FindBin;
my $script_dir = $FindBin::RealBin;

my %opts;
getopts('i:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i Genome_folder\n";
my $output_file = $input_folder.".length.txt";

chdir $input_folder or die $!;
my @files = glob("*.fas");
chdir $script_dir;

open OUT, ">$output_file" or die $!;
my %species;
foreach my $i (@files){
	my %seq = &read_fasta("$input_folder/$i") or die $!;
	foreach my $i (sort keys %seq){
		$species{$i} = "";
	}
	my @seq = values %seq;
	my $seq_length = length($seq[0]);
	$i =~ s/\.fas//;
	print OUT "$i\t$seq_length\n";
}
close OUT;

open OUT1, ">$input_folder.all_species.txt";
foreach my $i (sort keys %species){
	print OUT1 "$i\n";
}


sub read_fasta{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	open SEQ_INFO, ">>Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
	my %sequences;
	my $id;
	my $count = 0;
	my $seq_length;
	while (<READ_SEQ>){
	chomp;
	s/\r//;
	if ($_ =~ /^>(?<seq_name>.+$)/){
		$id  = $+{seq_name};
		$count += 1;
	}else{
		$sequences{$id} .= $_;
		$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences.\n";
	print SEQ_INFO "Successfuly read @_, in total $count sequences.\n";
	close SEQ_INFO;
	#print "\n";
	%sequences;
}
