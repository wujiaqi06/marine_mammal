#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
#For protein
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use FindBin;

my %opts;
getopts('i:j:', \%opts);
my $genome_folder = $opts{'i'} or die "use: $0 -i genome_folder -j species_assembly_to_latin_file\n";
my $species_assembly_to_latin_file = $opts{'j'} or die "use: $0 -i genome_folder -j species_assembly_to_latin_file\n";

my $script_dir = $FindBin::RealBin;
chdir $genome_folder or die $!;
my @files = glob ("*.gz");
chdir $script_dir;

my %names_assembly;
my %names_species;
open NAME, $species_assembly_to_latin_file or die $!;
while (<NAME>){
	chomp;
	my @line = split (/\t/, $_);
	$names_assembly{$line[0]} = $line[1];
	$names_species{$line[0]} = $line[2];
}

my $output_folder1 = $genome_folder."_assembly_name";
my $output_folder2 = $genome_folder."_species_name";
my $output_folder3 = $genome_folder."_copy_number";
mkdir $output_folder1;
mkdir $output_folder2;
mkdir $output_folder3;

foreach my $file (@files){
	my %seq = &read_fasta("$genome_folder/$file");
	my @file_name = split(/\./, $file);
	my $gene_name = $file_name[1];
	my %seq_assembly;
	my %seq_assembly_count;
	my %seq_speciese_name;
	foreach my $i (keys %seq){
		my $seq_name = $i;
		$seq_name =~ s/\s.*$//;
		$seq_name =~ s/vs_//;
		if (exists ($names_assembly{$seq_name})){
			if (exists $seq_assembly{$names_assembly{$seq_name}}){
				#print "Find $seq_name\n";
				$seq_assembly_count{$names_assembly{$seq_name}} ++;
			} else {
				$seq_assembly{$names_assembly{$seq_name}} = $seq{$i};
				$seq_assembly_count{$names_assembly{$seq_name}} = 1;
			}
			$seq_speciese_name{$names_species{$seq_name}} = $seq{$i};
		} else {
			print "$seq_name\n";
		}
	}

	open OUT1, ">$output_folder1/$gene_name.fas";
	open OUT3, ">$output_folder3/$gene_name.txt";
	foreach my $i (sort keys %seq_assembly){
		print OUT1 ">$i\n$seq_assembly{$i}\n";
		print OUT3 "$i\t$seq_assembly_count{$i}\n";
	}
	close OUT1;
	close OUT3;

	open OUT2, ">$output_folder2/$gene_name.fas";
	foreach my $i (sort keys %seq_speciese_name){
		print OUT2 ">$i\n$seq_speciese_name{$i}\n";
	}
	close OUT2;
}

sub read_fasta{
	open READ_SEQ, "gzip -dc @_ |" or die "Can not read Sequence file: $!\n";
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
	close READ_SEQ;
	%sequences;
}
