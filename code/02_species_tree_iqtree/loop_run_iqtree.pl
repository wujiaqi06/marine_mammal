#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
#version 3
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy qw(move);

my @all_seq = glob "./failed/*.fas";

my $seq_number = scalar @all_seq;

for (;$seq_number != 0;){
	move "autoIQTree_nofolder.pl", "./failed/autoIQTree_nofolder.pl";
	#move "iqtree", "./failed/iqtree";
	move "seq_length.txt", "./failed/seq_length.txt";
	move "_partition.txt", "./failed/_partition.txt";
	chdir "./failed";
	system "chmod 777 iqtree";
	system "perl autoIQTree_nofolder.pl -o run";
	@all_seq = glob "./failed/*.fas";
	$seq_number = scalar @all_seq;
}
