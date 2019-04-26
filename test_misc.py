import os,sys
import glob

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

# create_promoters("test_proms.bed")

def bed_intersect(file1, file2, out):
    '''
    Bedtools intersect file2 with file1 keeping the associated metadata
    
    Output Style:

    '''
    bed_command = [
    "bedtools",
    "intersect",
    "-a",
    file1,
    "-b",
    file2,
    "-wb",
    "-sorted",
    ">",
    out
    ]
    os.system(" ".join(bed_command))

# peaks = glob.glob(pwd + "/pass_peaks/*")
b1 = pwd + "/testb1.bed"
b2 = pwd + "/testb2.bed"
out = pwd + "/test_bintersect.bed"

bed_intersect(b1,b2,out)