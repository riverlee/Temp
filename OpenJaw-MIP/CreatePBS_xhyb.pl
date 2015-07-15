#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Mon Nov 17 13:39:00 2014
###################################
use strict;
use warnings;
use Cwd;

my $gene = $ARGV[0];

my $pbs_folder = "${gene}_pbs";
my $result_folder="${gene}_results";
my $currentdir=getcwd;

mkdir $pbs_folder if (! -e $pbs_folder);
mkdir $result_folder if (! -e $result_folder);

my $script="/home/jiali/MIP_Designer/xhyb_check.pl";

#### PBS template
my $pbsDesc=<<PBS;
#!/bin/bash
#Beginning of PBS bash script
#PBS -l walltime=10:00:00
#You must specify Wall Clock time (hh:mm:ss) [Maximum allowed 30 days = 720:00:00]
PBS


while(<DATA>){
    chomp;
    open(OUT,">${pbs_folder}/$_.sh") or die $!;
    print OUT $pbsDesc;
    print OUT "#pbs -N $_\n";

    print OUT "source /nfs/share/env/perl.env\n";
    print OUT "cd $currentdir/$result_folder\n";

    print OUT "$script --hash_length 12 --query_file ../${gene}.homs.fa --database_file /nfs/gp/Hs/2009_02-hg19/fa/$_.fa --allow_mismatches --pad_alignment --use_qnames --verbose >fasta_all.$_.matches.tab\n";
    close OUT;

}

__DATA__
chr1
chr10
chr11
chr11_gl000202_random
chr12
chr13
chr14
chr15
chr16
chr17
chr17_ctg5_hap1
chr17_gl000203_random
chr17_gl000204_random
chr17_gl000205_random
chr17_gl000206_random
chr18
chr18_gl000207_random
chr19
chr19_gl000208_random
chr19_gl000209_random
chr1_gl000191_random
chr1_gl000192_random
chr2
chr20
chr21
chr21_gl000210_random
chr22
chr3
chr4
chr4_ctg9_hap1
chr4_gl000193_random
chr4_gl000194_random
chr5
chr6
chr6_apd_hap1
chr6_cox_hap2
chr6_dbb_hap3
chr6_mann_hap4
chr6_mcf_hap5
chr6_qbl_hap6
chr6_ssto_hap7
chr7
chr7_gl000195_random
chr8
chr8_gl000196_random
chr8_gl000197_random
chr9
chr9_gl000198_random
chr9_gl000199_random
chr9_gl000200_random
chr9_gl000201_random
chrM
chrMT
chrUn_gl000211
chrUn_gl000212
chrUn_gl000213
chrUn_gl000214
chrUn_gl000215
chrUn_gl000216
chrUn_gl000217
chrUn_gl000218
chrUn_gl000219
chrUn_gl000220
chrUn_gl000221
chrUn_gl000222
chrUn_gl000223
chrUn_gl000224
chrUn_gl000225
chrUn_gl000226
chrUn_gl000227
chrUn_gl000228
chrUn_gl000229
chrUn_gl000230
chrUn_gl000231
chrUn_gl000232
chrUn_gl000233
chrUn_gl000234
chrUn_gl000235
chrUn_gl000236
chrUn_gl000237
chrUn_gl000238
chrUn_gl000239
chrUn_gl000240
chrUn_gl000241
chrUn_gl000242
chrUn_gl000243
chrUn_gl000244
chrUn_gl000245
chrUn_gl000246
chrUn_gl000247
chrUn_gl000248
chrUn_gl000249
chrX
chrY
