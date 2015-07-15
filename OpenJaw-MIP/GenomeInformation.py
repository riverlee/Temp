#!/bin/env python
import pandas as pd
import os

class GeneInformation():
    gene_data = ""

    #def __init__(self, knownGeneFile="sources/knownGene", geneXrefFile="sources/kgXref", knownCannonicalFile="sources/knownCanonical"):
    def __init__(self):
        
        current_path = os.path.dirname(__file__)

        knownGeneFile = current_path + "/sources/knownGene"
        geneXrefFile = current_path + "/sources/kgXref"
        knownCannonicalFile = current_path + "/sources/knownCanonical"

        #Check for presence of files
        kc = pd.read_table(knownCannonicalFile)
        kg = pd.read_table(knownGeneFile)
        kgx = pd.read_table(geneXrefFile,low_memory=False)
        tmp = pd.merge(kc[["transcript","protein"]], kgx[["#kgID","mRNA","geneSymbol","refseq"]][~kgx.refseq.isnull()], left_on='transcript', right_on='#kgID', how='inner')
        gene_data = pd.merge(tmp, kg, left_on='transcript', right_on='#name', how='inner')
        gene_data["char_chr"] = gene_data["chrom"].str.replace("chr","")
        self.gene_data = gene_data

    def __call__(self):
        return self.gene_data

    def makeGeneBed(self, gene_name, out_file):
        """
        This function writes the cannonical gene boundaries in bed format
        """
        self.gene_data[["char_chr","cdsStart","cdsEnd","geneSymbol"]][self.gene_data.geneSymbol == gene_name].to_csv(out_file, sep="\t", index=False, header=False)

    def makeGeneExonBed(self, gene_name, out_file):
        """
        This function writes the cannonical gene exon boundaries in bed format
        """
        file_data = ""
        file_data += self.gene_data[self.gene_data.geneSymbol == gene_name].apply(self._make_exon_string, axis=1)
        output = open(out_file, "w")
        for file_line in file_data:
            output.write(file_line)
        output.close()
        return ""

    def _make_exon_string(self, row):
        output = []
        chromosome, gene_name, num_exons = row['chrom'], row['geneSymbol'], row["exonCount"]
        exon_set = zip(row['exonStarts'].split(","), row['exonEnds'].split(","))

        for exon_index, (exon_start, exon_end) in enumerate(exon_set):
            if exon_start != "" and exon_end != "":
                output.append("{0}\t{1}\t{2}\t{3}.{4}\n".format(chromosome[3:], exon_start, exon_end, gene_name, exon_index + 1))
        return "".join(output)
    
    # added by Jiang
    def geneExist(self, gene_name):
        """
        This function check whether the gene name exists
        """
        return gene_name in set(self.gene_data.geneSymbol)



if __name__ == "__main__":
    main()

