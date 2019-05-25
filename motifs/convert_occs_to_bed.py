"""
Author: Saideep Gona

Converts raw motif occurences into bed file format with absolute 
coordinates instead of relative to regulatory region
"""

import sys
import os
import argparse
import glob
import numpy as np

pwd = os.getcwd()

parser = argparse.ArgumentParser(description="Convert raw motif occurrences \
into bed format with absolute genomic corordinates")
parser.add_argument("occurence_dir", help="Directory containing motif occurences",
                    default= pwd + "/motif_occurences/", nargs='?')
parser.add_argument("all_regulatory_regions", help="Directory containing regulatory \
region occurences", default= pwd + "/all_regs/allregregionsmod.bed", nargs='?')
parser.add_argument("output_dir", help="Output directory for bed file",
default= pwd + "/../pass_motifs/", nargs='?')

args = parser.parse_args()

def sort_in_place(f):
    temp_dir = pwd + "/tmp/"
    temp_f = f+".tmp"

    os.system("sort -k 1,1 -k2,2n -T " +temp_dir+" "+f+" > "+temp_f)
    os.system("rm "+f)
    os.system("mv "+temp_f+" "+f)

# Create dictionary mapping regulatory region to genome coords

reg_map = {}
with open(args.all_regulatory_regions, "r") as arr:
    for line in arr:
        p_line = line.rstrip("\n").split("\t")
        region = p_line[-1]
        reg_map[region] = {
            "chr": p_line[0],
            "start": p_line[1],
            "end": p_line[2]
        }
# print(args.occurence_dir)
motif_occs_files = glob.glob(args.occurence_dir + "*")
# print(motif_occs_files)
for mof in motif_occs_files:
    with open(mof, "r") as m:
        file_suf = mof.split("/")[-1]
        if os.path.exists(args.output_dir + file_suf):
            print(args.output_dir + file_suf, " exists")
            continue
        with open(args.output_dir + file_suf, "w") as out:

            line_count = 0
            for line in m:
                try:
                    if line_count == 0:
                        line_count += 1
                        continue

                    p_l = line.rstrip("\n").split("\t")
                    # print(p_l)
                    # motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
                    cur_reg = p_l[2]
                    chrom = reg_map[cur_reg]["chr"]
                    new_start = str(int(reg_map[cur_reg]["start"]) + int(p_l[3]))
                    new_end = str(int(reg_map[cur_reg]["start"]) + int(p_l[4]))
                    dif = str(abs(int(new_start) - int(new_end)))    
                    summit = str(int(int(new_start) + int(dif)/2))
                    strand = p_l[5]
                    score = p_l[6]
                    p_value = str(-1 * np.log10(float(p_l[7])))
                    q_value = p_l[8]
                    bed_line = [
                        chrom,
                        new_start,
                        new_end,
                        dif,
                        summit,
                        score,
                        p_value,
                        q_value
                    ]

                    out.write("\t".join(bed_line) + '\n')
                except Exception as e:
                    print(e)
                    print(line)