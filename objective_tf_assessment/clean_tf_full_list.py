'''
Author: Saideep Gona

This script is intended to filter the updated TF list from 
https://www.sciencedirect.com/science/article/pii/S0092867418301065?via%3Dihub#app2with 
to only contain proper TFs
'''

import os,sys,glob

# IO ************************************************************

pwd = os.getcwd()

if len(sys.argv) > 1:
    clean_file_path = sys.argv[1]
    output_path = sys.argv[2]
else:
    clean_file_path = "all_tfs_full_sheet_important_columns.tsv"
    output_path = "all_tfs_list.txt"

# IO END *********************************************************
with open(output_path, "w") as out:
    with open(clean_file_path, "r") as clean:
        for line in clean:
            p_line = line.rstrip("\n").split("\t")
            if p_line[1] == "Yes":
                out.write(p_line[0] + "\n")

    

