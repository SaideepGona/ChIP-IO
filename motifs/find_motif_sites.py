'''
Author: Saideep Gona

This script searches the genome for predicted motif binding sites based on meme binding matrices. It
incorporates epigenetic priors.
'''

import os, sys
import glob 
import pickle

pwd = os.getcwd()

# reference_genome = pwd+"/../GRCh38/GRCh38.p12.genome.fa"
reference_genome = pwd + "/all_regs/allregs.fa"
occurences_dir = pwd + "/motif_occurences/"

# GENERAL ANALYSIS



# TISSUE SPECIFIC ANALYSIS

# metadata_path = pwd + "/../pass_metadata/metadata.pkl"     

# metadata_dict = pickle.load(open(metadata_path, "rb"))      # Contains metadata on all datasets
# print(metadata_dict)

motifs_dir = pwd + "/motif_dir/"        # Find all the available meme style motifs
motif_files = glob.glob(motifs_dir+"*")
all_motif_tfs = set([x.split("/")[-1].rstrip(".meme") for x in motif_files])        

# chip_tfs_file = pwd + "/../static_lists/all_tfs_cur.txt"    # Find all the tfs which already have ChIP-Seq Data
# chip_tfs = set()
# with open(chip_tfs_file, "r") as ctfs:
#     for line in ctfs:
#         chip_tfs.add(line.rstrip("\n"))

# tfs_without_chip = all_motif_tfs.difference(chip_tfs)

priors_dir = pwd + "/priors/final_priors/"
prior_files_all = glob.glob(priors_dir+"*")
prior_files_wig = []
for pf in prior_files_all:
    if pf[-3:] == "wig":
        prior_files_wig.append(pf)
print(len(prior_files_wig))
# prior_tissues = set([x.split("/")[-1].rstrip(".wig") for x in prior_files])

p_thresh="4"

def exe(command):
    joint = " ".join(command)
    print("executing: " + joint)
    os.system(joint)

def find_motifs_for(motif_file, tissue_prior, ref, occurences_dir):
    '''
    Finds a set of motifs associated with a specific tf using associated priors for a specific tissue 
    on a reference sequence.
    '''

    tissue_prior_dist = tissue_prior.rstrip(".wig") + ".dist"
    tf = motif_file.split("/")[-1].rstrip(".meme")
    tissue = tissue_prior.split("/")[-1].rstrip(".wig")
    output_dir = occurences_dir + tf+"_"+tissue+"_occurences.tsv"

    motif_find = [
        "fimo",
        # "--o",
        # output_dir,
        # "--psp",
        # tissue_prior,
        # "--prior-dist",
        # tissue_prior_dist,
        "--text",
        "--thresh",
        "1e-"+p_thresh,
        "--verbosity",
        "4",
        motif_file,
        ref,
        ">",
        output_dir
    ]

    if os.path.isfile(output_dir):
        print(output_dir)
        return

    exe(motif_find)

# Cycle through all motifs and tissues
tissue = "general"
for motif in motif_files:
    find_motifs_for(motif, tissue, reference_genome, occurences_dir)
