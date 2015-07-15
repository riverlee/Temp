#!/bin/env python
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Mon Dec  1 16:45:10 2014
###################################
from nib_frag import nibFragger
import argparse
import sys

##MIP result format
#Mip_name   upstream    gap     downstream  H2  gap H1  upstreamSeq downStreamSeq H2Seq H1Seq   DesignSeq
parser = argparse.ArgumentParser(description = "Check the MIP design output file")

parser.add_argument("--input_file",dest = "input_file", help="MIP Design input file in tsv format", required = True)
parser.add_argument("--chrom", dest="chrom",help="Chromosome name", required=True)

try:
    args = parser.parse_args()
except:
    print("")
    parser.print_help()
    sys.exit(1)

input_file = args.input_file
chrom = args.chrom

MIP_backbone = "TTCCAACCTTCGATCTGTGCUUUGCCGCTCCGAGAACTTG"
smMIP_tag = "NNNNNNNNNNNN"
Five_Phos = "/5Phos/"


def revcomp(dna, reverse=True):
    bases = 'ATGCTACG'
    complement_dict = {bases[i]:bases[i+4] for i in range(4)}
    if reverse:
        dna = reversed(dna)
    result = [complement_dict[base] for base in dna]
    return ''.join(result)


fh = open(input_file,"r")

### Read each line
for line in fh:
    if line.startswith("Mip_Name"):
        continue
    else:
        (mipName, upstream,gapfill,downstream,h2pos,gapfill2,h1pos,strand,upstreamSeq,downstreamSeq,H2Seq,H1Seq,design) = line.rstrip().split("\t")

        flag = "OK"
        msg="";
        if strand == "+":
            if upstream != h2pos:
                flag = "Fail"
                msg+="Up not equal H2 (pos) |"

            if downstream !=h1pos:
                flag = "Fail"
                msg+="Down not equal H1 (pos) |"

            if upstreamSeq != H2Seq :
                flag = "Fail"
                msg+="Up not equal H2 (Seq) |"

            if downstreamSeq !=H1Seq:
                flag = "Fail"
                msg+="Down not equal H1 (Seq) |"
        else:
            if upstream != h1pos:
                flag = "Fail"
                msg+="Up not equal H1 (- pos) |"

            if downstream !=h2pos:
                flag = "Fail"
                msg+="Down not equal H2 (- pos) |"

            if upstreamSeq != revcomp(H1Seq) :
                flag = "Fail"
                msg+="Up not equal H1 (- Seq) |"

            if downstreamSeq !=revcomp(H2Seq):
                flag = "Fail"
                msg+="Down not equal H2 (-Seq) |"

        if gapfill != gapfill2:
            flag = "Fail";
            msg+="gap not equal |"

    
    ##based on upstream coordinate to excat seq
    upStart,upStop = upstream.split("-")
    downStart,downStop = downstream.split("-")

    myUpSeq = nibFragger(chrom,upStart,int(upStop)-int(upStart)+1)
    myDownSeq = nibFragger(chrom,downStart,int(downStop)-int(downStart)+1)

    if myUpSeq != upstreamSeq:
        flag = "Fail"
        msg+="Upstream not equal to fetched |"

    if myDownSeq != downstreamSeq:
        flag = "Fail"
        msg +="Downstream not equal to fetched |"

    myDesign="{}{}{}{}{}".format(Five_Phos,H1Seq,MIP_backbone,smMIP_tag,H2Seq)
    if myDesign != design:
        flag="Fail"
        msg+="Design not match |"
    print("{}\t{}\t{}".format(flag,mipName,msg));
