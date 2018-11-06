'''
Author: Saideep Gona

This script is intended to: 
1.) Clean study metadata
2.) Filter out non-TF associated studies
3.) Standardize TF names between studies
4.) Establish some pseudonym matching and modify metadata accordingly
5.) Output finalized files into the correct directories
This should be run as an initialization step prior to constructing database. 
'''

import os,sys,glob
from fuzzywuzzy import fuzz
import pickle
import glob

# IO ************************************************************

pwd = os.getcwd()

if len(sys.argv) > 1:
    tf_standard_path = sys.argv[1]
    metadata_dir = pwd + sys.argv[2]
    pass_metadata_dir = pwd + sys.argv[3]                   #"pass" implies that study targets a valid tf
    all_peaks_dir = pwd + sys.argv[4]
    pass_peaks_dir = pwd + sys.argv[5]
else:
    tf_standard_path = pwd + "/all_tfs_list.txt"
    metadata_dir = pwd + "/metadata/tissue_types.pkl"
    pass_metadata_dir = pwd + "/pass_metadata/metadata.pkl"                     #"pass" implies that study targets a valid tf
    all_peaks_dir = pwd + "/peaks/"
    pass_peaks_dir = pwd + "/pass_peaks/"

# IO END *********************************************************

class CleanTF():
    '''
    Holds different methods for cleaning a procured tf-name, assessing validity, and providing database-ready
    files to the appropriate folders    
    '''

    def __init__(self, tf_name_standard_file, fr):
        self.tf_standard = set()
        with open(tf_name_standard_file, "r") as tf:
            for l in tf:
                self.tf_standard.add(l.rstrip("\n"))
        self.fuzz_ratio = fr
        self.conversion_mappings = {}
        self.GEO_no_match = {}
        self.GEO_no_match_count = 0
        self.ENC_no_match = {}
        self.ENC_no_match_count = 0


    removable_prefixes = [
        "gfp",
        "-gfp",
        "gfp-",
        "-GFP",
        "GFP-",
        "-anti",
        "-eGFP",
        "eGFP-",
        "FLAG-",
        "_isoform1",
        "_isoform2",
        "scl/",
        "3xFLAG-"
    ]

    removable_suffixes = [
        "-tag",
        "-Antibody",
        "(THermo",
        "(abcam",
        "(anderson",
        ":",
        ";",
        "-",
        "_isoform1",
        "_isoform2",
        "/3",
        "(Bethyl"
    ]

    # removable_suffixes = []
    # removable_prefixes = []

    direct_conversions = {
        "C/EBPB": "CEBPB"
    }

    def within_standard_fuzzy_mapping_gen(self):

        mapping = {}

        for tf in self.tf_standard:
            match_set = set()
            for tf1 in self.tf_standard:
                if tf == tf1:
                    continue
                elif float(fuzz.ratio(tf, tf1)) > self.fuzz_ratio:
                    match_set.add(tf1)
            mapping[tf] = match_set
        self.within_standard_mapping = mapping

    def remove_pre_suf(self, intf):
        in_tf = intf[:]
        in_tf = in_tf.strip()
        for pref in self.removable_prefixes:
            in_tf = in_tf.lstrip(pref)
        for suf in self.removable_suffixes:
            in_tf = in_tf.rstrip(suf)
        
        return in_tf

    def check_if_close_match_tf(self, intf):
        
        for tf in self.tf_standard:
            if intf == tf:
                return [True, tf, 0]
        for tf in self.tf_standard:
            if float(fuzz.token_sort_ratio(intf, tf)) == 100.0:
                print(intf, tf, " token sort success")
                return [True, tf, 1]
        for tf in self.tf_standard:
            if float(fuzz.ratio(intf, tf)) > self.fuzz_ratio:
                print(intf, tf, " fuzzy ratio:"+str(fuzz.ratio(intf, tf)))                
                return [True, tf, 2]
        
        return [False, None, 3]

    def clean_meta(self, metadata_dir, pass_metadata_dir):
        '''
        Creates a cleaned version of the metadata dictionary
        '''
        with open(metadata_dir, 'rb') as pass_m:
            prev_meta = pickle.load(pass_m)        

        clean_meta = {}                          # Clean the individual data strings and output a clean metadata dict
        for key, meta in prev_meta.items():
            new_clean_meta = {}
            s_clean = self.remove_pre_suf(meta["tf"])
            new_clean_meta["tf"] = s_clean
            new_clean_meta["tissue"] = meta["tissue"]
            clean_meta[key] = new_clean_meta
        
        self.clean_meta = clean_meta

    def process_metadata(self, metadata_dir, pass_metadata_dir):

        '''
        Checks to see if the current TFs are within the tf standard list, or if they are similar enough to
        be converted to a known entry in the list. Saves a final metadata dict with these conversions
        '''

        pass_meta = {}

        for key, meta in self.clean_meta.items():
            new_meta = {}
            close_match = self.check_if_close_match_tf(meta["tf"])
            if close_match[0]:                  # If there is a close match as defined by the checker
                print(key)
                if meta["tf"] in self.within_standard_mapping[close_match[1]]:      # If the similarity can be explained as both tfs being similar standard tf names, then just use the default rather than corrected similar
                    print(close_match[1], " in standard_mapping")
                    new_meta["tf"] = meta["tf"]
                else:   # If the original is not in the standard set, then replace it with its match in the standard set
                    new_meta["tf"] = close_match[1]
                    if new_meta["tf"] != close_match[1]:    # Add the matching to stored mappings if thee match is not equivalent
                        self.conversion_mappings[meta["tf"]] = close_match[1]
            else:   # When there is no similar match in the matching set
                print("No Match", meta["tf"])

                if key[0] == "E":
                    if meta["tf"] in self.ENC_no_match:     # Add to the stored mismatches 
                        self.ENC_no_match[meta["tf"]] += 1
                        self.ENC_no_match_count += 1
                    else:
                        self.ENC_no_match[meta["tf"]] = 1
                        self.ENC_no_match_count += 1
                else:
                    if meta["tf"] in self.GEO_no_match:     # Add to the stored mismatches 
                        self.GEO_no_match[meta["tf"]] += 1
                        self.GEO_no_match_count += 1
                    else:
                        self.GEO_no_match[meta["tf"]] = 1 
                        self.GEO_no_match_count += 1                  
                
                continue                            # Skip to next study if no match

            new_meta["tissue"] = meta["tissue"]

            pass_meta[key] = new_meta

        self.pass_meta = pass_meta

        with open(pass_metadata_dir, 'wb') as pass_m:
            pickle.dump(pass_meta, pass_m, pickle.HIGHEST_PROTOCOL)

    def move_peaks(self, all_peaks_dir, pass_peaks_dir):
        '''
        Copies over only peaks which exist in the final metadata table to the "pass peaks" directory
        '''
        pass_cont = glob.glob(pass_peaks_dir + "/*")
        for peaks in pass_cont:
            os.remove(peaks)

        valid_studies = list(self.pass_meta.keys())
        valid_study_files = [(x + "_peaks.xls") for x in valid_studies]

        for study in valid_study_files:
            os.system("cp " + all_peaks_dir + study + " " + pass_peaks_dir)
        


standard_cleaner = CleanTF(tf_standard_path, 80.0)
standard_cleaner.within_standard_fuzzy_mapping_gen()
standard_cleaner.clean_meta(metadata_dir, pass_metadata_dir)
standard_cleaner.process_metadata(metadata_dir, pass_metadata_dir)
print("conversion mappings")
print(standard_cleaner.conversion_mappings)
print("no match set ENCODE")
print(standard_cleaner.ENC_no_match)
print(standard_cleaner.ENC_no_match_count)
print("no match set GEO")
print(standard_cleaner.GEO_no_match)
print(standard_cleaner.GEO_no_match_count)

standard_cleaner.move_peaks(all_peaks_dir, pass_peaks_dir)
peak_set = set()
for key, meta in standard_cleaner.pass_meta.items():
    peak_set.add(meta["tf"])
print("peak list")
print(peak_set)
print(len(peak_set))
