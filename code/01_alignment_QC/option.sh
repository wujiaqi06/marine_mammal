1, RemoveGap_cds.pl
Anno: remove alignment parts with coverage less than 70% of the species
RemoveGap_cds.pl del0.3
perl RemoveGap_cds.pl -i Alignment_Finished -o mammal127_gap0.7 -d 0.7 -t fas

2, RemoveFragmentSeq_cds.pl
Anno: remvoe fragment sequences in each species, with length less than 10aa, one side of its neighbour gaps less than 5 aa space.
RemoveFragmentSeq_cds.pl f10g5
perl RemoveFragmentSeq_cds.pl -i mammal127_gap0.3 -o mammal127_gap0.3f10g5 -f10 -g 5

3, RemoveFragmentSeq_cds.pl
RemoveFragmentSeq_cds.pl f3g1
Anno: remvoe fragment sequences in each species, with length less than 3aa, one side of its neighbour gaps less than 1 aa space.
perl RemoveFragmentSeq_cds.pl -i mammal127_gap0.3f10g5 -o mammal127_gap0.3f10g5_f3g1 -f 3 -g 1

4, RemoveGap_cds.pl
RemoveGap_cds.pl del0.3
Anno: after RemoveGap_cds.pl, the coverange of sequences changed. Again remove alignment parts with coverage less than 30% species
perl RemoveGap_cds.pl -i mammal127_gap0.3f10g5_f3g1 -o mammal127_gap0.3f10g5_f3g1_gap0.3 -d 0.7 -t fas

5.1 perl FragmentCount.pl -i mammal127_gap0.3f10g5_f3g1_gap0.3
Count fregment information. 
Example output:
Cov7027_gap0.3f10g5_f3g1_gap0.3.gap_info.txt
File Name	Seq Name	Fragment Start	Fragment End	Fragment Length	Gap length before Fragment	Gap length after Fragment
A2M.fas	acinonyx_jubatus	1	48	47	0	3
A2M.fas	acinonyx_jubatus	1930	2106	177	9	12
A2M.fas	acinonyx_jubatus	2119	2514	396	12	9
A2M.fas	acinonyx_jubatus	2524	4431	1908	9	0
A2M.fas	acinonyx_jubatus	304	1920	1617	12	9
A2M.fas	acinonyx_jubatus	52	291	240	3	12

5.2 PairwisePi.pl
Anno: Do summary statistics of fregment, then count the pairwise p-distance of fregment of each sequence to all other sequences.
Count the p-diseance of each fragment to all other sequences.
Actually did not used for secreening.
Cov7027_gap0.3f10g5_f3g1_gap0.3.p_distance.txt
Example output
File Name	Sequence1	Sequence2	OverlapStart	OverlapEnd	OverlapLength	Diff	Ident	P_Distance
A2M.fas	acinonyx_jubatus	ailuropoda_melanoleuca	1	48	48	22	26	0.458333333333333
A2M.fas	acinonyx_jubatus	aotus_nancymaae	4	36	33	22	11	0.666666666666667
A2M.fas	acinonyx_jubatus	aotus_nancymaae	40	48	9	2	70.222222222222222
A2M.fas	acinonyx_jubatus	balaenoptera_acutorostrata_scammoni	1	48	48	21	27	0.4375
A2M.fas	acinonyx_jubatus	bison_bison_bison	1	48	48	22	26	0.458333333333333
A2M.fas	acinonyx_jubatus	bos_indicus	1	21	21	12	90.571428571428571
A2M.fas	acinonyx_jubatus	bos_mutus	1	48	48	22	26	0.458333333333333
A2M.fas	acinonyx_jubatus	bos_taurus	1	48	48	24	24	0.5
A2M.fas	acinonyx_jubatus	bubalus_bubalis	1	48	48	27	21	0.5625

5.3 FragmentCount.pl FragmentScreen.pl
Anno: FragmentCount.pl summarize the fragments for each sequences of each species.
FragmentScreen.pl: summarize fragments with less than certain size of amino acid length, normally 3aa length.
This screen is only for the sequences in the beginning or the ends of the sequences.

6. RemoveSmallFrag.pl Remove small fragments
RemoveSmallFrag.pl -i mammal127_gap0.3f10g5_f3g1_gap0.3 -o mammal127_gap0.3f10g5_f3g1_gap0.3_rmF3 -r mammal127_gap0.3f10g5_f3g1_gap0.3.gap_info.txt_FragLess3.txt
Anno: based on the output of FragmentScreen.pl, remove these small fragments.

7. RemoveEmptySequences.pl
perl RemoveEmptySequences.pl -i mammal127_gap0.3f10g5_f3g1_gap0.3_rmF3 -o mammal127_gap0.3f10g5_f3g1_gap0.3_rmF3_rmEmpty0.7 -m 0.7
Anno: count the missing percent of each sequences, for example, Species X had 40% of the sequences is gap in geneA.
Remove sequences of the species with its gap more than gap percent.

Example:
Successfuly read ./mammal127_gap0.3f10g5_f3g1_gap0.3_rmF3/ZSCAN5B.fas, in total 94 sequences with 1512 bp.
Removed One Sequence: capra_hircus, 0.740079365079365 #74% of this sequence is gap in alignment.
Removed One Sequence: elephantulus_edwardii, 0.720238095238095

Notice: this will lose many sequences in terminal branch in gene tree.

