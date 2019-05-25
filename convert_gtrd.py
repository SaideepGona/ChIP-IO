'''
Author: Saideep Gona

Script converts bulk download files of ChIP-Seq Peaks 
and DNase footprints to a split file form
'''

import os, sys
import pickle
import glob

pwd = os.getcwd()
gtrd_raw_peak_file = pwd + "/gtrd_raw_peaks.tsv"
gtrd_peaks_dir = pwd + "/pass_peaks/"
tf_list = pwd + "/objective_tf_assessment/all_tfs_list.txt"
metadata_path = pwd + "/pass_metadata/metadata.pkl"

gtrd_raw_footprint_dir = pwd + "/gtrd_raw_footprints/"
gtrd_raw_metadata_file = pwd + "/metadata_footprints.txt"
gtrd_footprints_dir = pwd + "/gtrd_footprints/"

motif_occ_dir = pwd + "/pass_motifs/"
pass_motif_dir = pwd + "/pass_motifs/"

all_tfs = set()

def sort_in_place(f):
    temp_dir = pwd + "/tmp/"
    temp_f = f+".tmp"
    print("sort -k 1,1 -k2,2n -T " +temp_dir+" "+f+" > "+temp_f)
    os.system("sort -k 1,1 -k2,2n -T " +temp_dir+" "+f+" > "+temp_f)
    os.system("rm "+f)
    os.system("mv "+temp_f+" "+f)

def remove_trailing_tabs(f):
    os.system("sed -i 's/[\t]*$//' "+f)

def add_to_bed(file_name, line):

    if os.path.exists(file_name):
        append_write = 'a' # append if already exists
    else:
        append_write = 'w' # make a new file if not

    with open(file_name, append_write) as out:
        out.write("\t".join(line)  + "\n")

def keep_intersect(source, filt, out_file):
    '''
    Use bedtools intersect to filter a file to include only
    intersectional regions from the source file
    '''

    intermediate = pwd + "/intermediates/cur_intermediate.txt"

    bed_command = [
    "bedtools",
    "intersect",
    "-a",
    source,
    "-b",
    filt,
    "-wb",
    "-sorted",
    ">",
    intermediate
    ]

    os.system(" ".join(bed_command))

    with open(intermediate, "r") as inter:
        for line in inter:
            p_line = line.rstrip("\n").split("\t")
            keep = p_line[0:7]
            add_to_bed(out_file, keep)

    os.remove(intermediate)


with open(tf_list, "r") as tfs:
    for line in tfs:
        tf = line.rstrip("\n")
        all_tfs.add(tf)

# Convert GTRD peaks to individual files

# Gtrd Columns:
#CHROM	START	END	-10*log10(pvalue)	FDR(%)	antibody	cellLine	experiment	fold_enrichment	summit	tags	tfClassId	tfTitle	treatment	uniprotId

# Macs Columns
#chrom  start   end length  summit  pileup  -log10(pvalue)  fold_enrichment -log10(qvalue)

# All peaks columns
#chrom  start   end peak_number length  summit  pileup  -log_p  fold_enrichment -logq   tf  tissue  accession

if False:
    metadata_dict = pickle.load(open(metadata_path, "rb"))
    os.system('rm '+gtrd_peaks_dir+"EXP*")
    with open(gtrd_raw_peak_file, "r") as raw_peaks:
        line_num = 0
        for line in raw_peaks:
            try:
                if line_num == 0:
                    line_num += 1
                    continue

                p_line = line.rstrip("\n").split("\t")
                # print(p_line)
                tf = p_line[12]
                if tf in all_tfs:
                    accession = p_line[7]
                    tissue = p_line[6]
                    print(tissue)
                    metadata_dict[accession] = {
                        "tissue": tissue,
                        "tf": tf
                    }
                    new_line = [
                        p_line[0],
                        p_line[1],
                        p_line[2],
                        str(line_num-1),
                        str(int(p_line[2])-int(p_line[1])),
                        p_line[9],
                        "1000",
                        p_line[3],
                        p_line[10],
                        "1000",
                        tf,
                        tissue,
                        accession
                    ]

                    f_name = gtrd_peaks_dir + accession + "_peaks.xls"
                    add_to_bed(f_name, new_line)

                    line_num +=1
            except Exception as e:
                print(e)
                print(line)
                continue

    with open(metadata_path, 'wb') as meta:
        pickle.dump(metadata_dict, meta, protocol=pickle.HIGHEST_PROTOCOL)

print("Done 1")

# FOOTPRINT PROCESSING

# Metadata Columns
# id	peaks_id	species	info_(treatment)	source	source_id	cell_type	cell_type_id	cellosaurus_id	cell_ontology_id	exp_factor_ontology_id	uberon_id	external references
if True:
    os.system('rm '+gtrd_footprints_dir+"*")
    raw_footprint_files = glob.glob(gtrd_raw_footprint_dir + "*")

    with open(gtrd_raw_metadata_file, "r") as meta:
        line_count = 0
        for line in meta:
            print(line_count)
            try:
                if line_count == 0:
                    line_count += 1
                    continue
                p_line = line.rstrip("\n").split("\t")
                if "Homo sapiens" not in p_line[2]:
                    continue
                if len(p_line[3]) < 1:
                    continue
                # pass_age = False
                # data = p_line[3].split(";")
                # for part in data:
                #     if "age" in part:
                #         pass_age = True
                #         if "days" in part:
                #             pass_age = False
                # if not pass_age:
                #     continue
                # print(line)
                for f in raw_footprint_files:
                    if p_line[1] in f:
                        # print(f)
                        tissue = p_line[4]
                        tissue = "_".join(tissue.split(" "))

                        # print(tissue)
                        
                        with open(f, "r") as fps:
                            for nline in fps:
                                np_line = nline.rstrip("\n").split("\t")
                                add_to_bed(gtrd_footprints_dir + tissue + ".bed", np_line[0:3])
                line_count+=1
            except Exception as e:
                print(e)
                print(line)
                continue

print("Done 2") 
#Create tissue-specific motif predictions

if True:
    motif_occ_files = glob.glob(motif_occ_dir+"*")
    print("Sorting files")
    for m_f in motif_occ_files:
        sort_in_place(m_f)
    footprint_files = glob.glob(gtrd_footprints_dir+"*")
    for f_f in footprint_files:
        sort_in_place(f_f)
    # os.system('rm '+pass_motif_dir+"*")
    print(motif_occ_files)
    print(footprint_files)
    for m_f in motif_occ_files:
        for f_f in footprint_files:
            tf = m_f.split("/")[-1].split("_")[0]
            tissue = f_f.split("/")[-1].split(".")[0]
            f_name = pass_motif_dir + "/" + tf + "_" + tissue + "_occurences.tsv"
            print(tf,f_name)
            keep_intersect(m_f, f_f, f_name)

