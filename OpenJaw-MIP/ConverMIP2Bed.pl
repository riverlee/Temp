#!/bin/env perl
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Mon Dec  1 10:08:29 2014
###################################
use strict;
use warnings;
###########################
## Take the MIP result and convert it to bed 12 format for visualization
##MIP result format
#Mip_name   upstream    gap     downstream  H2  gap H1  upstreamSeq downStreamSeq H2Seq H1Seq   DesignSeq
#
my $usage="ConvertMIP2Bed.pl mip.tsv name chrom\n";
my ($in,$name,$chr) = @ARGV;

if(@ARGV !=3){
    print $usage;
    exit 1;
}

## IN FORMAT
#mipName  upstream, gap ,downstream
#abcd  1-4  5-6  7-8

open(IN,$in) or die $!;
<IN>;
print "track name=$name itemRgb='On'\n";
#my $count=0;
while(<IN>){
    s/\r|\n//g;
    my($name,$upstream,$gap,$downstream) = split "\t";
    my $strand = "+";
    my $rgb="255,0,0";
    #   $count++;
    my (undef,$count) = split /\-/,$name;
    
    if($count % 2 ==0){
        $strand="-";
        $rgb="0,0,255";
    }
    my ($a,$b) = split /\-/,$upstream;
    my ($c,$d) = split /\-/,$downstream;
    my $start = $a-1;
    my $end = $d;

    my $blocksizeUp=$b-$a+1;
    my $blocksizeDown=$d-$c+1;
    my $blocksize = "$blocksizeUp,$blocksizeDown";

    my $start2=$c-$a;
    my $blockstart="0,$start2";

    print  join "\t",($chr,$start,$end,$name,0,$strand,$start,$end,$rgb,2, $blocksize,$blockstart);
    print "\n";
}
