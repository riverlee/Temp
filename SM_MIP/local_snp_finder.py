#!/usr/bin/python
import pandas as pd

class local_snp_id:
    snp_data = None
    cosmic_data = None

    def __init__(self, chrom):
        chrom = chrom.replace("chr","")
        self.snp_data = pd.read_table("/nfs/ss11/mghent/MIPSeq/sources/dbSNP_fill_chr{}.data".format(chrom), names=["rsID","chrom","position","strand","Alleles"])
        self.cosmic_data = pd.read_table("/nfs/ss11/mghent/MIPSeq/sources/Cosmic_chr{}_parsed.tsv".format(chrom))

    def local_snps(self, start_position, stop_position):
        # use >= or =< instead of > and < to avoid missing snps in homs
        data_rsid = list(self.snp_data.rsID[(self.snp_data.position >= start_position) & (self.snp_data.position <= stop_position)].values)
        data_positions = list(self.snp_data.position[(self.snp_data.position >= start_position) & (self.snp_data.position <= stop_position)].values)
        data_set = zip(data_rsid, data_positions)
        # sorted by position
        data_set =  sorted(list(data_set),key=lambda tmp: tmp[1])
        out_string = ""
        for rsid, pos in data_set:
            out_string += "rs{}:{}|".format(rsid,pos)
        return out_string[:len(out_string)-1]

    def local_sms(self, start_position, stop_position):
        in_range = ((self.cosmic_data.start >= start_position) & (self.cosmic_data.start <= stop_position)) | ((self.cosmic_data.stop >= start_position) & (self.cosmic_data.stop <= stop_position))
        mutation_ids = list(self.cosmic_data.MutationID[in_range])
        starts = list(self.cosmic_data.start[in_range])
        stops = list(self.cosmic_data.stop[in_range])

        data_set = zip(mutation_ids, starts, stops)
        # sorted by position
        data_set =  sorted(list(data_set),key=lambda tmp: tmp[1])
        out_string = ""
        for mut_id, start, stop in (data_set):
            out_string += "{}:{}-{}|".format(mut_id, start, stop)
        return out_string[:len(out_string)-1]

