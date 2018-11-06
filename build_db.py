'''
Author: Saideep Gona

This script is intended to populate the given database from local files
'''

from flask_app import db
from flask_app import ChIP_Meta, Peaks, Presets
import glob
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import unicodedata

plt.ion()

# IO ************************************************************

pwd = os.getcwd()

if len(sys.argv) > 1:
    peak_dir = sys.argv[1]
    metadata_dir = sys.argv[2]
    metrics_dir = sys.argv[3]
else:
    peak_dir = pwd + "/temp_peaks/"
    metadata_path = pwd + "/pass_metadata/metadata.pkl"     
    metrics_dir =  pwd + "/static/images/"

# IO END *********************************************************

def find_duplicates(in_list):  
    unique = set(in_list)
    # print(unique)  
    for each in unique:  
        count = in_list.count(each)  
        if count > 1:  
            print ("duplicate: " + each + " " + count)

def histogram(field, data, metrics_dir):
    plt.rc('axes',edgecolor='white')
    plt.rc('lines', color='white')
    plt.rc('text', color='white')
    plt.rc('xtick', color='white')
    plt.rc('ytick', color='white')
    np_data = np.array(data)
    # print(np_data)
    # hist, bins = np.histogram(np_data, bins = 50)
    # print(hist,bins)
    plt.hist(data, color = 'white', bins = 50)
    plt.xlabel("Bin Ranges", color='white')
    plt.ylabel("Frequency", color='white')
    plt.title("Distribution of "+field+" values across all peaks in database")

    # for ax, color in zip([ax1, ax2, ax3, ax4], ['white', 'white', 'white', 'white']):
    #     for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
    #         ticks.set_color('white')
    # for pos in ['top', 'bottom', 'right', 'left']:
    #     plt.spines[pos].set_edgecolor('white')
    # plt.xlabel.set_color('white')
    # plt.ylabel.set_color('white')
    # plt.tick_params(axis='x', colors='white')
    # plt.tick_params(axis='y', colors='white')
    # plt.title.set_color('white')

    # # plt.plot(np.array([0,1,2,3]), color = 'blue')
    # plt.show()
    # print(x)
    savefile = metrics_dir + "/" + field + "_hist.png"
    print(savefile)
    plt.savefig(savefile, transparent=True)
    plt.cla()

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
            print(meta)
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
                print(meta)
                db.session.add(meta)


    elif type(tissue_obj) == str:                                   # Standard single tissue entry
        tissue = slugify(tissue_obj)
        print(tissue)

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

pileup = []
p_values = []
fold_enrichment = []
q_values = []

peak_files = glob.glob(peak_dir+"/*")
for p_f in peak_files:
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
                print(tissue_obj)
                if len(tissue_obj) == 0:
                    tissue = "NA"
                else:
                    for tissue_p in tissue_obj:
                        tissue = slugify(tissue_p)
                        peak = Peaks(
                            experiment_accession = p_f.rstrip("_peaks.xls").split("/")[-1],
                            tissue_types = tissue,
                            transcription_factors = metadata_dict_ref[exp_acc][1],
                            chrom = p_l[0],
                            start = str(int(p_l[1]) - 1),
                            end = str(int(p_l[2]) - 1),
                            length = p_l[3],
                            summit = str(int(p_l[4]) - 1),
                            pileup = p_l[5],
                            log_p = p_l[6],
                            fold_enrichment = p_l[7],
                            log_q = p_l[8]
                        )  
                        db.session.add(peak)

            elif type(tissue_obj) == str:                                   # Standard single tissue entry
                tissue = slugify(tissue_obj)

                # Peak columns:
                # chr   start   end length  summit  pileup  -log_p   fold_enrichment -log_q
                peak = Peaks(
                    experiment_accession = p_f.rstrip("_peaks.xls").split("/")[-1],
                    tissue_types = tissue,
                    transcription_factors = metadata_dict_ref[exp_acc][1],
                    chrom = p_l[0],
                    start = str(int(p_l[1]) - 1),
                    end = str(int(p_l[2]) - 1),
                    length = p_l[3],
                    summit = str(int(p_l[4]) - 1),
                    pileup = p_l[5],
                    log_p = p_l[6],
                    fold_enrichment = p_l[7],
                    log_q = p_l[8]
                )
                print(peak)
                db.session.add(peak)
    db.session.commit()

print(len(p_values), max(p_values), min(p_values))

peak_fields = {
    "pileup": pileup,
    "-log_p": p_values,
    "fold_enrichment": fold_enrichment,
    "-log_q": q_values
}

for field in peak_fields.keys():
    histogram(field, peak_fields[field], metrics_dir)

# END MAIN ********************************************************