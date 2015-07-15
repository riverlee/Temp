#!/usr/bin/python
import pandas as pd

class local_repeatMasker_id:
    repeatMasker_data = None

    def __init__(self, chrom):
        #chrom = chrom.replace("chr","")
        self.repeatMasker_data = pd.read_table("/home/jiali/reference/RepeatMasker/{}_RepeatMasker.txt".format(chrom))

    def local_repeatMaskers(self, start_position, stop_position):
        rm_data = []
        con1 = (self.repeatMasker_data.genoStart<start_position) & (self.repeatMasker_data.genoEnd>=start_position)
        con2 = (self.repeatMasker_data.genoStart>=start_position) & (self.repeatMasker_data.genoStart<stop_position)
        sub_repeatMasker_data = self.repeatMasker_data[con1 | con2]

        if not sub_repeatMasker_data.empty :
            for start, end in zip(sub_repeatMasker_data.genoStart, sub_repeatMasker_data.genoEnd):
                for pos in range(start,end):
                    if pos >=start_position and pos <=stop_position:
                        rm_data.append(pos)
        
        return rm_data
