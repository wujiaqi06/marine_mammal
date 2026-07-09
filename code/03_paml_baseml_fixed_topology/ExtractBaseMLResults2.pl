#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
#For protein
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use FindBin;
my $script_dir = $FindBin::RealBin;

my %opts;
getopts('i:o:', \%opts);
my $run_number = $opts{'i'} or die "use: perl $0 -i run_number -o output_file\n";
my $output_file = $opts{'o'} or die "use: perl $0 -i run_number -o output_file\n";

open OUT, ">$output_file" or die $!;

foreach my $i (1..$run_number){
	chdir "./run$i/BaseMLResult";
	my @files = glob ("*.out");
	chdir $script_dir;
	foreach my $j (@files){
		open IN, "./run$i/BaseMLResult/$j" or die $!;
		my $boot = 0;
		while (<IN>){
			chomp;
			if (/tree length/){
				$boot = 1;
			}

			if ((/Detailed output/)){
				$boot = 2;
			}

			if ($boot == 1){
				if ((/\d\.\d+/)&!(/tree length/)){
					print OUT "$j: $_\n";
				}
			}
		}
		close IN;
	}
}
