#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
# 2018.4.18
# Usage example: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align -d 0.3 -t phy
# Usage example2: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align_phy -d 0.3 -t fas
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy;

my %opts;
getopts('i:o:r:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -r remove_info\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -r remove_info\n";
my $remove_info = $opts{'r'} or die "use: perl $0 -i Genome_folder -o output_folder -r remove_info\n";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas';
#say @all_seq_file;
chdir "../";

open REM, "./$remove_info" or die "Cannot read $remove_info: $!\n";
my %rm_seq;
while (<REM>){
	chomp;
	my @lines = split ("\t", $_);
	$rm_seq{$lines[0]} .= "$_\n";
}

foreach my $seq (sort @all_seq_file){
	#print "seq is $seq\n";
	my @rm_info;
	if (exists $rm_seq{$seq}){
		#print "Find!\n";
		@rm_info = grep { /\S/ } split ("\n", $rm_seq{$seq});
		#print "rm_info is @rm_info\n";
	} else {
		copy ("./$input_folder/$seq", "./$output_folder/$seq");
	}

	if (scalar (@rm_info) != 0){
		my %read_sequences = &read_fasta("./$input_folder/$seq");
		my @seq_name = sort keys %read_sequences;
		foreach my $i (@rm_info){
			my @lines = split ("\t", $i);
			my $rm_length = $lines[3] - $lines[2] + 1;
			if ($rm_length % 3 != 0){
				$rm_length = int ($rm_length / 3) + 1;
			}
			my $replace = "-" x $rm_length;
			substr ($read_sequences{$lines[1]}, $lines[2]-1, $rm_length, $replace);
			#print "length is $rm_length\n";
		}
		open OUT, ">./$output_folder/$seq" or die "Cannot create $output_folder/$seq: $!\n";
		foreach my $i (@seq_name){
			print OUT ">$i\n$read_sequences{$i}\n";
		}
	}

}



sub read_fasta{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	open SEQ_INFO, ">>../Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
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
	print "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	print SEQ_INFO "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	close SEQ_INFO;
	#print "\n";
	%sequences;
}
close REM;
close OUT;

