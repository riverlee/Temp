##############################
### MIPSEQ Design pipeline ###
##############################


## Usage

### Step 1. Make bed files of all the exons for a given gene using canonical gene model

Bed_make.py will invoke the GenomeInformation.py

```
python Bed_make.py --gene_name BRCA1 --output_file BRCA1.Exons.bed
```

### Step 2. Make all the possible homologous regions (homs) by hom_designer.py
This step will use modules from 
dust.py
local_snp_finder.py
mip_tm_clip.py
nib_frag.py


### Step 3. Run xhyd to search hits of homs genome wide
xhyb_check.pl --hash_length 12 --query_file ../yourgene.homs.fa --databbase_file /nfs/gp/Hs/2009-02-hg19/fa/chr22.fa --allow_mismatches --pad_alignment --use_qnames --verbose > fasta_all.chr22.matches.tab


### Step 4. Add homs hits by hom_hits.py
hom_hit.py 

### Step 5. Output MIP
select_MIPs.py


##############################################
# Use MIPPineline.pl yourgene to create the scripts for each step 


