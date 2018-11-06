'''
Author: Saideep Gona

This script takes in enhancer-gene datasets downloaded from EnhancerAtlas and
adds a column of "gene name" to the files based on the Ensemble id 
'''

import os, sys
import glob
from pyliftover import LiftOver

lo = LiftOver('hg19', 'hg38')

pwd = os.getcwd()

mapping_file = pwd + "/ensemble_mapping"
enhancer_dir = pwd
processed_dir = pwd + "/processed/"

pwd = os.getcwd()

mapping_dict = {}
with open(mapping_file, "r") as map_f:
    for line in map_f:
        line_p = line.rstrip("\n").split("\t")
        mapping_dict[line_p[0]]=line_p[2]
        mapping_dict[line_p[1]]=line_p[2]
print(mapping_dict)
enhancer_files = glob.glob(enhancer_dir+"/*.txt")
print(enhancer_files)

for e_file in enhancer_files:
    if e_file.endswith(".py") or e_file.endswith("mapping"):
        continue
    with open(e_file, "r") as e:
        e_file_name = e_file.split("/")[-1]
        with open(processed_dir+"converted_"+e_file_name, "w") as p:
            for line in e:
                p_line = line.rstrip("\n").split("\t")
                print(p_line)
                candidate = p_line[6]
                if candidate not in mapping_dict:
                    continue
                new_name = mapping_dict[candidate]
                new_line = p_line[0:6] + [new_name]
                new_line[1] = str(lo.convert_coordinate(p_line[0], int(p_line[1])))
                new_line[2] = str(lo.convert_coordinate(p_line[0], int(p_line[2])))
                new_line = "\t".join(new_line) + "\n"
                p.write(new_line)

