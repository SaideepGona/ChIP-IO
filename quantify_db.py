'''
Author: Saideep Gona

This script will go through the current database and 
output summary statistics
'''

import os, sys
import glob
import pickle

pwd = os.getcwd()
metadata_file = pwd + "/pass_metadata/metadata.pkl"

motifs_dir = pwd + "/motifs/"
jaspar_glob = glob.glob(motifs_dir + "/motif_dir/*")
hoco_glob = glob.glob(motifs_dir + "/hocomoco/individual/*")

dnase_dir = pwd + "/gtrd_footprints/"
dnase_glob = glob.glob(dnase_dir+"*")

enhancer_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
enhancer_glob = glob.glob(enhancer_dir+"*")

with open(metadata_file, 'rb') as pass_m:
    meta = pickle.load(pass_m) 
    # print(meta)
studies = set()
studies_in_house = set()
studies_GEO_in_house = set()
studies_ENCODE_in_house = set()
studies_gtrd = set()

unique_tfs = set()
unique_tfs_in_house = set()
unique_tfs_GEO_in_house = set()
unique_tfs_ENCODE_in_house = set()
unique_tfs_gtrd = set()
unique_tfs_chip = set()

unique_tfs_motifs_JASPAR = set()
unique_tfs_motifs_HOCO = set()
unique_tfs_motifs = set()

unique_tissues = set()
unique_tissues_in_house = set()
unique_tissues_GEO_in_house = set()
unique_tissues_ENCODE_in_house = set()
unique_tissues_gtrd = set()
unique_tissues_chip = set()

unique_tissues_dnase = set()
unique_tissues_enhancers = set()

for study in meta.keys():
    if study[0:2] == "EX":
        studies.add(study)
        studies_gtrd.add(study)
        unique_tfs.add(meta[study]["tf"])
        unique_tfs_gtrd.add(meta[study]["tf"])
        unique_tfs_chip.add(meta[study]["tf"])
        unique_tissues = unique_tissues.union(set(meta[study]["tissue"]))
        unique_tissues_gtrd = unique_tissues_gtrd.union(set(meta[study]["tissue"]))
        unique_tissues_chip = unique_tissues_chip.union(set(meta[study]["tissue"]))

    elif study[0:2] == "EN":
        studies.add(study)     
        studies_in_house.add(study)
        studies_ENCODE_in_house.add(study)
        unique_tfs.add(meta[study]["tf"])
        unique_tfs_in_house.add(meta[study]["tf"])
        unique_tfs_ENCODE_in_house.add(meta[study]["tf"])
        unique_tfs_chip.add(meta[study]["tf"])
        unique_tissues = unique_tissues.union(set(meta[study]["tissue"]))
        unique_tissues_in_house = unique_tissues_in_house.union(set(meta[study]["tissue"]))
        unique_tissues_ENCODE_in_house = unique_tissues_ENCODE_in_house.union(set(meta[study]["tissue"]))
        unique_tissues_chip = unique_tissues_chip.union(set(meta[study]["tissue"]))

    elif study[0:2] == "GS":
        studies.add(study)     
        studies_in_house.add(study)
        studies_GEO_in_house.add(study)
        unique_tfs.add(meta[study]["tf"])
        unique_tfs_in_house.add(meta[study]["tf"])
        unique_tfs_GEO_in_house.add(meta[study]["tf"])
        unique_tfs_chip.add(meta[study]["tf"])
        unique_tissues = unique_tissues.union(set(meta[study]["tissue"]))
        unique_tissues_in_house = unique_tissues_in_house.union(set(meta[study]["tissue"]))
        unique_tissues_GEO_in_house = unique_tissues_GEO_in_house.union(set(meta[study]["tissue"]))
        unique_tissues_chip = unique_tissues_chip.union(set(meta[study]["tissue"]))

for motif in jaspar_glob:
    unique_tfs.add(motif.split("/")[-1].rstrip(".meme"))
    unique_tfs_motifs.add(motif.split("/")[-1].rstrip(".meme"))
    unique_tfs_motifs_JASPAR.add(motif.split("/")[-1].rstrip(".meme"))

for motif in hoco_glob:
    unique_tfs.add(motif.split("/")[-1].rstrip(".meme"))
    unique_tfs_motifs.add(motif.split("/")[-1].rstrip(".meme"))
    unique_tfs_motifs_HOCO.add(motif.split("/")[-1].rstrip(".meme"))

for dnase_file in dnase_glob:
    unique_tissues.add(dnase_file.split("/")[-1].rstrip(".meme"))
    unique_tissues_dnase.add(dnase_file.split("/")[-1].rstrip(".meme"))

for enhancer_file in enhancer_glob:
    unique_tissues.add(enhancer_file.split("/")[-1].rstrip(".meme"))
    unique_tissues_enhancers.add(enhancer_file.split("/")[-1].rstrip(".meme"))

# studies = set()
print("Studies: ", len(studies))
# studies_in_house = set()
print("Studies in-house: ", len(studies_in_house))
# studies_GEO_in_house = set()
print("Studies GEO in-house: ", len(studies_GEO_in_house))
# studies_ENCODE_in_house = set()
print("Studies ENCODE in-house: ", len(studies_ENCODE_in_house))
# studies_gtrd = set()
print("Studies GTRD: ", len(studies_gtrd))
print()

# unique_tfs_chip = set()
print("Unique TFs All ChIP-Seq: ", len(unique_tfs_chip))
# unique_tfs_in_house = set()
print("Unique TFs in-house: ", len(unique_tfs_in_house))
# print(unique_tfs_in_house)
# unique_tfs_GEO_in_house = set()
print("Unique TFs GEO in-house: ", len(unique_tfs_GEO_in_house))
# unique_tfs_ENCODE_in_house = set()
print("Unique TFs ENCODE in-house: ", len(unique_tfs_ENCODE_in_house))
# unique_tfs_gtrd = set()
print("Unique TFs GTRD: ", len(unique_tfs_gtrd))
# print(unique_tfs_gtrd)
print()

# unique_tfs_motifs = set()
print("Unique TF Motifs: ", len(unique_tfs_motifs))
# unique_tfs_JASPAR = set()
print("Unique TF Motifs JASPAR: ", len(unique_tfs_motifs_JASPAR))
# unique_tfs_HOCO = set()
print("Unique TF Motifs HOCOMOCO: ", len(unique_tfs_motifs_HOCO))
print()

# unique_tfs = set()
print("Unique TFs: ", len(unique_tfs))
print()

# unique_tissues_chip = set()
print("Unique Tissues All ChIP-Seq: ", len(unique_tissues_chip))
# unique_tissues_in_house = set()
print("Unique Tissues in-house: ", len(unique_tissues_in_house))
# unique_tissues_GEO_in_house = set()
print("Unique Tissues GEO in-house: ", len(unique_tissues_GEO_in_house))
# unique_tissues_ENCODE_in_house = set()
print("Unique Tissues ENCODE in-house: ", len(unique_tissues_ENCODE_in_house))
# unique_tissues_gtrd = set()
print("Unique Tissues GTRD: ", len(unique_tissues_gtrd))
print()

# unique_tissues_dnase = set()
print("Unique Tissues DNase-Seq: ", len(unique_tissues_dnase))
print()

# unique_tissues_enhancers = set()
print("Unique Tissues Enhancers: ", len(unique_tissues_enhancers))
print()

# unique_tissues = set()
print("Unique Tissues: ", len(unique_tissues))

# print(unique_tfs)