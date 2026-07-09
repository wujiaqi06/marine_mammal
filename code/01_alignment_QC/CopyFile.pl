#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
#For protein
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy qw(copy);

my %opts;

getopts('i:j:o:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i input_folder -j input_geneList -o output_folder\n";
my $input_geneList = $opts{'j'} or die "use: $0 -i input_folder -j input_geneList -o output_folder\n";
my $output_folder = $opts{'o'} or die "use: $0 -i input_folder -j input_geneList -o output_folder\n";

my @gene_list;
open LIST, "./$input_geneList" or die "Cannot read $input_geneList: $!\n";
while (<LIST>){
	chomp;
	s/\s//g;
	my $file = "$_";#.".fas";
	push @gene_list, $file;
}

mkdir "./$output_folder" or die "Cannot create folder $output_folder: $!\n";
foreach my $file (@gene_list){
	copy "./$input_folder/$file", "./$output_folder/$file";
}

close LIST;