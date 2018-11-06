'''
Author: Saideep Gona

Converts EnhancerAtlas-based enhancer-gene mapping info to a GTF.
'''

import os, sys

pwd = os.getcwd()
enhancer_file = "unionAHSIL_norepeat"
gtf_out = "enhancers.gtf"

with open(enhancer_file, "r") as ef:
    with open(gtf_out, "w") as gtf:
        for line in ef:
            p_line = line.rstrip("\n").split("\t")
            new_line = [p_line[0],
                        "EnhancerAtlasUnion",
                        "exon",
                        p_line[1],
                        p_line[2],
                        ".",
                        ".",
                        ".",
                        "."
                        ]
            new_line_join = "\t".join(new_line) + "\n"
            gtf.write(new_line_join)
