#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy qw(copy);

my %opts;
getopts('i:o:n:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i Genome_folder -o output_folder -n number_of_folder\n";
my $output_folder = $opts{'o'} or die "use: $0 -i Genome_folder -o output_folder -n number_of_folder\n";
my $number_of_folder = $opts{'n'} or die "use: $0 -i Genome_folder -o output_folder -n number_of_folder\n";


chdir "./$input_folder";
my @sequence_file = glob "*.phy";
print "@sequence_file\n";
chdir "../";

mkdir "$output_folder" or die "Can not create folder $output_folder: $!\n";
#chdir "./$output_folder" or die "Can not enter folder $output_folder: $!\n";

my $count = 1;
my $folder_number;
my $sequence_number = $#sequence_file+1;
my $per_folder = int($sequence_number/$number_of_folder) + 1;
mkdir "./$output_folder/$count" or die "Can not make folder $count: $!";
#mkdir "./$count" or die "Can not make folder ./$output_folder/$count: $!\n";
#chdir "../";
for my $i (0..$#sequence_file){
	copy "./$input_folder/$sequence_file[$i]", "./$output_folder/$count/$sequence_file[$i]";
	if (($i+1)%$per_folder == 0){
		$count++;
		#chdir "./$output_folder";
		mkdir "./$output_folder/$count" or die "Can not make folder ./$output_folder/$count: $!\n";
		#chdir "../";
	}
}

print "Made $count folders\n";
