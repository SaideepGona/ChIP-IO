'''
Author: Saideep Gona

This script is intended to populate the given database from local files
'''

from application_for_build_db import db
from application_for_build_db import ChIP_Meta, Peaks, Presets
import glob
import os
import sys
import pickle
import numpy as np
# import matplotlib
# matplotlib.use('Agg')

# import matplotlib.pyplot as plt
import unicodedata
import multiprocessing

# plt.ion()
# plt.use('Agg')

# IO ************************************************************

pwd = os.getcwd()
num_cpus = multiprocessing.cpu_count()

if len(sys.argv) > 1:
    peak_dir = sys.argv[1]
    metadata_dir = sys.argv[2]
    metrics_dir = sys.argv[3]
    peak_file = sys.argv[4]
else:
    peak_dir = pwd + "/pass_peaks/"
    metadata_path = pwd + "/pass_metadata/metadata.pkl"     
    metrics_dir =  pwd + "/static/images/"
    peak_file = pwd + "/all_peaks.tsv"

if os.path.isfile(peak_file):
    os.remove(peak_file)
    os.system("touch "+peak_file)
else:
    os.system("touch")

# IO END *********************************************************

def find_duplicates(in_list):  
    unique = set(in_list)
    # print(unique)  
    for each in unique:  
        count = in_list.count(each)  
        if count > 1:  
            print ("duplicate: " + each + " " + count)

# def histogram(field, data, metrics_dir):
#     plt.rc('axes',edgecolor='white')
#     plt.rc('lines', color='white')
#     plt.rc('text', color='white')
#     plt.rc('xtick', color='white')
#     plt.rc('ytick', color='white')
#     np_data = np.array(data)

#     plt.hist(data, color = 'white', bins = 50)
#     plt.xlabel("Bin Ranges", color='white')
#     plt.ylabel("Frequency", color='white')
#     plt.title("Distribution of "+field+" values across all peaks in database")

#     savefile = metrics_dir + "/" + field + "_hist.png"
#     print(savefile)
#     plt.savefig(savefile, transparent=True)
#     plt.cla()

def slugify(value):
    """
    Mini version of slugify which just converts spaces to underscores
    """
    new_val = ""
    for char in value:
        if char == " ":
            new_val += "_"
        else:
            new_val += char
    return new_val

def setify(in_list):
    return list(set(in_list))

# MAIN ***********************************************************

db.create_all()

# Read in metadata
metadata_dict = pickle.load(open(metadata_path, "rb"))
metadata_dict_ref = {}

find_duplicates(list(metadata_dict.keys()))
print(len(list(metadata_dict.keys())), " NUMBER OF STUDIES")
# sys.exit()

# print(metadata_files)
for m_f in list(metadata_dict.keys()):
    # print(m_f)

    tissue_obj = metadata_dict[m_f]["tissue"]                       # If list of tissues, creates duplicate entries for each 
    if type(tissue_obj) == list:
        tissue_obj = setify(tissue_obj)
        if len(tissue_obj) == 0:
            tissue = "NA"
            meta = ChIP_Meta(
                experiment_accession = m_f,
                tissue_types = tissue,
                transcription_factors = metadata_dict[m_f]["tf"]
            )

            metadata_dict_ref[m_f] = [tissue_obj, metadata_dict[m_f]["tf"]]
            # print(meta)
            db.session.add(meta)

        else:
            tissue_obj = setify(tissue_obj)
            for tissue_p in tissue_obj:
                tissue = slugify(tissue_p)
                meta = ChIP_Meta(
                    experiment_accession = m_f,
                    tissue_types = tissue,
                    transcription_factors = metadata_dict[m_f]["tf"]
                )
                metadata_dict_ref[m_f] = [tissue_obj, metadata_dict[m_f]["tf"]]
                # print(meta)
                db.session.add(meta)


    elif type(tissue_obj) == str:                                   # Standard single tissue entry
        tissue = slugify(tissue_obj)
        # print(tissue)

        meta = ChIP_Meta(
            experiment_accession = m_f,
            tissue_types = tissue,
            transcription_factors = metadata_dict[m_f]["tf"]
        )

        metadata_dict_ref[m_f] = [tissue_obj, metadata_dict[m_f]["tf"]]
        print(meta)
        db.session.add(meta)
    
    else:
        continue

db.session.commit()

# Read in all peak files


with open(peak_file, "a") as p_file:

    peak_id_count = 0

    pileup = []
    p_values = []
    fold_enrichment = []
    q_values = []

    full_peak_array = []

    peak_files = glob.glob(peak_dir+"/*")
    p_f_count = 1

    for p_f in peak_files:
        print(p_f_count)
        p_f_count += 1
        if peak_id_count > 100:
            break
        with open(p_f, "r") as pre_f:
            f = pre_f.readlines()[24:]
            for line in f:
                p_l = line.rstrip("\n").split("\t")
                if len(p_l) != 10:
                    continue
                try:
                    int(p_l[1])
                except:
                    continue
                exp_acc = p_f.rstrip("_peaks.xls").split("/")[-1]
                if exp_acc not in metadata_dict_ref:
                    print(exp_acc, " not in metadata!")
                    continue

                if p_l[1] == "0":
                    print(p_l)
                    sys.exit()
                pileup.append(float(p_l[5]))
                p_values.append(float(p_l[6]))
                fold_enrichment.append(float(p_l[7]))
                q_values.append(float(p_l[8]))

                # Handle multiple tissue-specification entries

                tissue_obj = metadata_dict_ref[exp_acc][0] 
                if type(tissue_obj) == list:
                    # print(tissue_obj)
                    if len(tissue_obj) == 0:
                        tissue = "NA"
                    else:
                        for tissue_p in tissue_obj:
                            tissue = slugify(tissue_p)
                            write_dict = {
                                "experiment_accession": p_f.rstrip("_peaks.xls").split("/")[-1],
                                "tissue_types": tissue,
                                "transcription_factors": metadata_dict_ref[exp_acc][1],
                                "chrom": p_l[0],
                                "start": str(int(p_l[1]) - 1),
                                "end": str(int(p_l[2]) - 1),
                                "length": p_l[3],
                                "summit": str(int(p_l[4]) - 1),
                                "pileup": p_l[5],
                                "log_p": p_l[6],
                                "fold_enrichment": p_l[7],
                                "log_q": p_l[8],
                                "id": peak_id_count
                            }
                            peak_id_count+=1
                            peaks_column_list = [
                                "chrom",
                                "start",
                                "end",
                                "id",
                                "length",
                                "summit",
                                "pileup",
                                "log_p",
                                "fold_enrichment",
                                "log_q",
                                "transcription_factors",
                                "tissue_types",
                                "experiment_accession",
                                "id"
                                ]
                            write_list = [str(write_dict[x]) for x in peaks_column_list]
                            p_file.write("\t".join(write_list)+"\n")
                            # full_peak_array.append(write_list)
                            # peak = Peaks(
                            #     experiment_accession = p_f.rstrip("_peaks.xls").split("/")[-1],
                            #     tissue_types = tissue,
                            #     transcription_factors = metadata_dict_ref[exp_acc][1],
                            #     chrom = p_l[0],
                            #     start = str(int(p_l[1]) - 1),
                            #     end = str(int(p_l[2]) - 1),
                            #     length = p_l[3],
                            #     summit = str(int(p_l[4]) - 1),
                            #     pileup = p_l[5],
                            #     log_p = p_l[6],
                            #     fold_enrichment = p_l[7],
                            #     log_q = p_l[8]
                            # )  
                            # db.session.add(peak)

                elif type(tissue_obj) == str:                                   # Standard single tissue entry
                    tissue = slugify(tissue_obj)

                    # Peak columns:
                    # chr   start   end length  summit  pileup  -log_p   fold_enrichment -log_q
                    write_dict = {
                        "experiment_accession": p_f.rstrip("_peaks.xls").split("/")[-1],
                        "tissue_types": tissue,
                        "transcription_factors": metadata_dict_ref[exp_acc][1],
                        "chrom": p_l[0],
                        "start": str(int(p_l[1]) - 1),
                        "end": str(int(p_l[2]) - 1),
                        "length": p_l[3],
                        "summit": str(int(p_l[4]) - 1),
                        "pileup": p_l[5],
                        "log_p": p_l[6],
                        "fold_enrichment": p_l[7],
                        "log_q": p_l[8],
                        "id": peak_id_count
                    }
                    peak_id_count += 1
                    peaks_column_list = [
                        "chrom",
                        "start",
                        "end",
                        "id",
                        "length",
                        "summit",
                        "pileup",
                        "log_p",
                        "fold_enrichment",
                        "log_q",
                        "transcription_factors",
                        "tissue_types",
                        "experiment_accession",
                        "id"
                        ]
                    write_list = [str(write_dict[x]) for x in peaks_column_list]
                    p_file.write("\t".join(write_list)+"\n")
                    # full_peak_array.append(write_list)
                    # peak = Peaks(
                    #     experiment_accession = p_f.rstrip("_peaks.xls").split("/")[-1],
                    #     tissue_types = tissue,
                    #     transcription_factors = metadata_dict_ref[exp_acc][1],
                    #     chrom = p_l[0],
                    #     start = str(int(p_l[1]) - 1),
                    #     end = str(int(p_l[2]) - 1),
                    #     length = p_l[3],
                    #     summit = str(int(p_l[4]) - 1),
                    #     pileup = p_l[5],
                    #     log_p = p_l[6],
                    #     fold_enrichment = p_l[7],
                    #     log_q = p_l[8]
                    # )
                    # print(peak)   
                    # db.session.add(peak)
        # db.session.commit()

print(num_cpus)
os.system("split --number=l/"+str(num_cpus)+" "+peak_file+" "+peak_file+"_")

peak_fields = {
    "pileup": pileup,
    "-log_p": p_values,
    "fold_enrichment": fold_enrichment,
    "-log_q": q_values
}



# for field in peak_fields.keys():
#     histogram(field, peak_fields[field], metrics_dir)

# END MAIN ********************************************************