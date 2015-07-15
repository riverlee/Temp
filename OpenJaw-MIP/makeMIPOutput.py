#!/bin/env python
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Mon Dec  1 16:45:10 2014
###################################
from nib_frag import nibFragger
import argparse
import sys

## MIP result format
## Mip_name   upstream    gap     downstream  strand and other
## Will only use the first 5 columns information
parser = argparse.ArgumentParser(description = "Based on upstream and downstream coordiantes to make MIPs design output")

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

## Print header
print("Mip_Name\tupstream\tgap_fill_positions\tdownstream\tH2_positions\tgap_fill_positions\tH1_positions\tstrand\tupstream_seq\tdownstream_seq\tDesign_H2_seq\tDesign_H1_seq\tDesign_MIP")

### Read each line

for line in fh:
    if line.startswith("Mip_Name"):
        continue
    else:
        #tmp = line.rstrip().split("\t")
        mipName, upstream,gapfill,downstream = line.rstrip().split("\t")[0:4]
    
    strand = "+"
    mm, count = mipName.split("-")
    count=int(count)
    if count % 2 == 0:
        strand = "-"

    ##based on upstream coordinate to excat seq
    upStart,upStop = upstream.split("-")
    downStart,downStop = downstream.split("-")

    myUpSeq = nibFragger(chrom,upStart,int(upStop)-int(upStart)+1)
    myDownSeq = nibFragger(chrom,downStart,int(downStop)-int(downStart)+1)

    h2pos=""
    h1pos=""
    
    H1Seq=""
    H2Seq=""

    if strand == "+":
        h2pos=upstream
        h1pos=downstream
        H2Seq = myUpSeq
        H1Seq = myDownSeq
    else:
        h1pos=upstream
        h2pos = downstream
        H2Seq = revcomp(myDownSeq)
        H1Seq = revcomp(myUpSeq)

    ## Rewrite gap fill
    aa = int(upStop)+1
    bb = int(downStart)-1
    gapfill = "{}-{}".format(aa,bb)

    myDesign="{}{}{}{}{}".format(Five_Phos,H1Seq,MIP_backbone,smMIP_tag,H2Seq)
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(mipName,upstream,gapfill,downstream,h2pos,gapfill,h1pos,strand,myUpSeq,myDownSeq,H2Seq,H1Seq,myDesign));
