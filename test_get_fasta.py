import os,sys

pwd = os.getcwd()
reference = pwd + "/GRCh38/GRCh38.p12.genome.fa"
bed_file = "/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/enhancer-gene/processed/tissue_beds/blood.bed"
intermediate_file = "test_getfasta.txt"
convert_bed_fasta = [
    "bedtools",
    "getfasta",
    "-fi",
    reference,
    "-bed",
    bed_file,
    "-fo",
    intermediate_file
    ]

os.system(" ".join(convert_bed_fasta))