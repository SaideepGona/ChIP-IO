'''
Author: Saideep Gona

This script is intended to do any processing related to TF Motifs. Run from one directory above target dir.
Currently supports meme format motifs from JASPAR database
'''
import os,sys
import requests,json
import glob
import pickle

pwd = os.getcwd()
target_dir = pwd + "/motif_dir/"

def get_request(motif_id):

    headers = {'accept': 'application/json'}
    url = "http://jaspar.genereg.net/api/v1/matrix/" + motif_id + "/"
    response = requests.get(url, headers = headers)
    response_dict = response.json()

    return response_dict

def check_parent(string):
    # print(string)
    start = -1
    end = -1
    for ind in range(len(string)):
        if string[ind] == "(":
            start = ind
        if string[ind] == ")":
            end = ind
    if start != -1:
        print(string, start)
        return string[0:start]
    else:
        return string

motif_files_all = glob.glob(target_dir + "/*.meme")
motifs_all = [x.split("/")[-1].rstrip(".meme") for x in motif_files_all]

for f in motif_files_all:
    f_m = f.split("/")[-1].rstrip(".meme")
    print(f_m)
    if not ((f_m[0:3] == "MA0") or (f_m[0:3] == "MA1")):
        # print(f)
        print("hm")
        os.remove(f)

motif_files = glob.glob(target_dir + "/*.meme")

motifs = [x.split("/")[-1].rstrip(".meme") for x in motif_files]
# print(motifs)

id_motif_mapping = {}
motif_count = 0

# Rename individual motif files with their gene names, and create multiple copies for aliases

for motif in motifs:

    m_file = target_dir + "/" + motif + ".meme"

    with open(m_file, "r") as m:
        lines = m.read().splitlines()
        g_names = lines[9].split(" ")[-1].split("::")
        g_names = [x.upper() for x in g_names]
        for g_name_r in g_names:
            g_name = check_parent(g_name_r)
            print(g_name)
            cp_motif = [
                "cp",
                target_dir + "/" + motif + ".meme",
                target_dir + "/" + g_name + ".meme"
            ]    
            # print(cp_motif)       
            os.system(" ".join(cp_motif))


# for motif in motifs:
#     resp = get_request(motif)
#     try:
#         symbol = resp["symbol"]
#     except:
#         print(motif)
#         continue
#     # id_motif_mapping[resp] = symbol
#     # print(symbol)
#     cp_motif = [
#         "cp",
#         pwd + "/" + motif + ".meme",
#         pwd + "/" + symbol + ".meme"
#     ]
#     os.system(" ".join(cp_motif))
#     motif_count += 1
#     print(motif + " copied, "+ str(motif_count))
