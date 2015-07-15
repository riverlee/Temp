#!/bin/env python
import sys
import argparse

# Customerized modules
from mylog import info,error

## Take a gene symbol name and output all its exon in bed format

parser = argparse.ArgumentParser(description = "Export gene's exons in bed format from it's canonical model")

parser.add_argument("--gene_name", dest = "gene_name", help = "Gene Symbol", required = True)
parser.add_argument("--output_file", dest = "output_file", help = "Output file name", required = True)

try:
    args = parser.parse_args()
except:
    print("")
    parser.print_help()
    sys.exit(1)


gene_name = args.gene_name
output_file = args.output_file

import GenomeInformation as gi

ginfo = gi.GeneInformation()

## check whether gene exists
if ginfo.geneExist(gene_name):
    info("Your input is '{}' and output is '{}'".format(gene_name,output_file))
    ginfo.makeGeneExonBed(gene_name,output_file)
    info("Finished")
else:
    error("Gene '"+gene_name+"' not exists\n")

