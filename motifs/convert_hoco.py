"""
Author: Saideep Gona

Convert a whole hocomoco meme format file into individual files
"""

import sys, os, argparse




parser = argparse.ArgumentParser(description="Converts whole hocomoco meme file into individual files")
parser.add_argument("whole_meme_path", help="Path to whole meme file to be split up")
parser.add_argument("output_dir", help="Output directory for meme files")
args = parser.parse_args()


boiler = """
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25
"""

with open(args.whole_meme_path, "r") as wmp:
    line_c = 0
    blocks = []
    cur_block = ""
    for line in wmp:
        # print(line)
        if line_c < 9:
            line_c += 1
            continue

        if line == "\n":
            cur_block_s = cur_block.split("\n")
            tf = cur_block_s[0].split(" ")[1].split("_")[0]
            new_file_path = args.output_dir + tf + ".meme"
            with open(new_file_path, "w") as nfp:
                nfp.write(boiler + "\n" + cur_block)
            cur_block = ""
        else:
            cur_block += line




print(boiler)