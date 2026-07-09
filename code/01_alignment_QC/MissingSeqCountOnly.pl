#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i Genome_folder -o output_folder\n";
my $output_folder = $opts{'o'} or die "use: $0 -i Genome_folder -o output_folder\n";

#my $input_folder = "test";
##my $output_folder = "test2";
#mkdir "./$output_folder" or die "Cannot make folder $output_folder: $!\n";
chdir "./$input_folder" or die "Cannot enter folder $input_folder: $!\n";
my @sequences = glob "*.fas";
#chdir "../";

open SPECIES, "../species.txt" or die "Cannot open file species.txt: $!";
my @species_name;
while (<SPECIES>){
	chomp;
	push @species_name, $_;
}
#print "@species_name\n";
open MISS, ">.././$output_folder.MissingInfo.txt" or die "Can not make file output_folder_MissingSum.txt: $!\n";
open MISSCOUNT, ">.././$output_folder.MissingCount.txt" or die "Can not make file output_folder_MissingSum.txt: $!\n";
print MISSCOUNT "Sequence\tMissingNumber\n";
print MISS "Genes\tMissingSpecies\n";
foreach my $seq (@sequences){
	#open OUTF, ">>.././$output_folder/$seq" or die "Cannot make output file $seq: $!\n";
	print MISS "$seq\t";
	my %read_sequence = &read_fasta($seq);
	my @seq_name = sort keys %read_sequence;
	my $count = 0;
	foreach my $species (@species_name){
		if (!(exists ($read_sequence{$species}))){
			#$read_sequence{$species} = "-" x 12;
			#print OUTF ">$species\n$read_sequence{$species}\n";
			print MISS "$species\t";
			$count ++;
		} else {
			#print OUTF ">$species\n$read_sequence{$species}\n";
		}
	}
	print MISS "\n";
	print MISSCOUNT "$seq\t$count\n";
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
