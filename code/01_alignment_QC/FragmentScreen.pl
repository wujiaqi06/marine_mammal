#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:j:', \%opts);
my $input_file = $opts{'i'} or die "use: perl $0 -i gap_info_file -j output_fragment_length\n";
my $length = $opts{'j'} or die "use: perl $0 -i gap_info_file -j output_fragment_length\n";

open GAP, "./$input_file" or die "Cannot open $input_file: $!\n";
open OUT, ">./${input_file}_FragLess${length}.txt";
my $count = 1;
while (<GAP>){
	chomp;
	if ($count == 1){
		print OUT "$_\n";
		$count ++;
		} else {
			my @lines = split ("\t", $_);
			if ($lines[4] <= $length){
				print OUT "$_\n";
				$count ++;
			}
		} 
}
close GAP;
close OUT;
