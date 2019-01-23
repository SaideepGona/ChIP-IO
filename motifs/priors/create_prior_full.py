'''
Author:Saideep Gona

Takes in a bigwig file and processes it into a valid prior for use in fimo
'''
import os,sys

pwd = os.getcwd()

def exe(command):
    joint = " ".join(command)
    print("executing: " + joint)
    os.system(joint)


def edit_wig(valid_headers_list, edit_file, new_file):
    '''
    Remove sections from a wig file if they are not part of valid headers
    '''
    valid_headers = set(valid_headers_list)

    remove = False
    with open(edit_file, "r") as e_f:
        with open(new_file, "w+") as n_f:
            for line in e_f:
                if line[0] == "f":
                    print(line)
                    p_line = line.rstrip("\n").split(" ")
                    chrom = p_line[1].split("=")[1]
                    # print(chrom)
                    if chrom in valid_headers:
                        remove = False
                    elif chrom not in valid_headers:
                        print(chrom, " not valid and being removed")
                        remove = True
                        continue
                if not remove:
                    n_f.write(line)


def order_wig(valid_headers_list, edit_file, new_file):
    '''
    Reorder a wig file by chromosome
    '''
    headers = [("fixedStep chrom="+x+" start=1 step="+wig_step+" span="+wig_step+"\n") for x in valid_headers_list]
    remove = False
    current_header = ""
    with open(new_file, "w+") as n_f:
        first_line = read_first_line(edit_file)
        n_f.write(first_line)
        for header in headers:
            n_f.write(header)
            print(header, " written")
            with open(edit_file, "r") as e_f:
                for line in e_f:
                    if line[0] == "f":
                        print(line)
                        current_header = line
                        continue
                    
                    if current_header == header:
                        n_f.write(line)

    
def read_first_line(filename):
    with open(filename, "r") as f:
        line_num = 0
        string = ""
        for line in f:
            if line_num > 0:
                break
            string += line
            line_num += 1
        return string


valid_chroms = [("chr"+str(x)) for x in range(1,23)] + ["chrX", "chrY", "chrM"]
ref = "/home/saideep/Documents/GitHub_Repos/Saideep/ChIP-IO/GRCh38/GRCh38.p12.genome.fa" 
wig_step = "2"

# Convert starting bigwig to wig file(but it is actually outputting a bed-like file)
print("Converting bigwig to wiglike bed file")

tissue = "heart"
bigwig_file = pwd + "/" + tissue +"_dnase.bigWig"
bed1 = bigwig_file.rstrip(".bigWig") + "_raw.bed"

bigwig_to_bed = [
    "./bigWigToWig",
    bigwig_file,
    bed1
]
# exe(bigwig_to_bed)              # TOGGLE

# print("Removing invalid chromosomes")

# bed2 = bigwig_file.rstrip(".bigWig") + ".bed"

# # Convert bed-like to more proper wig file
# print("Converting bed file to wig")

# bed2wiggle = [
#     "bash",
#     "bedtowiggle.sh",
#     bed2,
#     wig
# ]
# exe(bed2wiggle)              # TOGGLE

# Bedgraph to wig

wig = bigwig_file.rstrip(".bigWig") + ".wig"
bedgraph_to_wig = [
    "perl",
    "bedgraph_to_wig.pl",
    "--bedgraph",
    bed1,
    "--wig",
    wig,
    "--step",
    wig_step
]
# exe(bedgraph_to_wig)              # TOGGLE

# Edit wig to be an ordered subset of reference

print("Removing unwanted chromosomes")
wig_edited = bigwig_file.rstrip(".bigWig") + "_edited.wig"
# edit_wig(valid_chroms, wig, wig_edited)              # TOGGLE
print("Sorting wig file to match reference")
wig_final = bigwig_file.rstrip(".bigWig") + "_final.wig"
# order_wig(valid_chroms, wig_edited, wig_final)              # TOGGLE

# Create priors
remove_dir = [
    "rm",
    "-rf",
    tissue
]
exe(remove_dir)
create_prior = [
    "create-priors",
    "-o",
    tissue,
    ref,
    wig_final
]
exe(create_prior)              # TOGGLE

priors_dir = pwd + "/priors/final_priors/"



# bed_to_wig = [
#     "awk",
#     '{if(NR>1) {if($1!=lastChrom){printf("variableStep chrom=%s\\n",$1);lastChrom=$1;}print $2,$4}}',
#     bed2,
#     ">",
#     wig
# ]
# exe(bed_to_wig)

# Remove lines with certain substrings

# substrings = [
#     "#",
#     "section",
#     "random"
# ]

# for ss in substrings:

#     remove_line = [
#         "sed",
#         "i.bak",
#         "'/" + ss + "/d'",
#         wig
#     ]
#     exe(remove_line)