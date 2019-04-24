'''
Author: Saideep Gona

This script maps the motif_ids found in the MEME format motif PWM files to the
names of their corresponding transcription factors.
'''

import os, sys
import pickle
import glob

pwd = os.getcwd()

motif_pwms_dir = pwd + "/motif_dir/"
mapping = {}

motif_files = glob.glob(motif_pwms_dir + "*")

def extract_motif_id(filename, dict):
    tf = filename.rstrip(".meme").split("/")[-1]
    with open(filename, "r") as mf:
        for line in mf:
            if line.startswith("MOTIF"):
                print(line)
                p_line = line.split(" ")
                motif_id = p_line[1]
                mapping[motif_id] = tf

for mf in motif_files:
    extract_motif_id(mf, mapping)
 
print(mapping)

with open("motif_mapping.pkl", "wb") as mmp:
    pickle.dump(mapping, mmp, pickle.HIGHEST_PROTOCOL)

            