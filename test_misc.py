import os,sys

pwd = os.getcwd()

def create_promoters(filename):

    ref_genes = pwd + "/annotations/genes.bed"

    with open(ref_genes, "r") as genes:
        with open(filename, "w") as tmp_proms:
            for line in genes:
                p_lines = line.rstrip("\n").split("\t")
                chrom = p_lines[0]
                rang = [int(p_lines[1]), int(p_lines[2])]
                strand = p_lines[3]
                g_name = p_lines[4]

                upstream = 10
                downstream = 20

                if strand == "+":
                    tss = min(rang)
                    new_range = [tss-upstream, tss+downstream]
                elif strand == "-":
                    tss = max(rang)
                    new_range = [tss-downstream, tss + upstream]
                
                new_line = [chrom, str(new_range[0]), str(new_range[1]), strand, g_name, p_lines[1], p_lines[2]]

                tmp_proms.write("\t".join(new_line) + "\n")

create_promoters("test_proms.bed")