'''
Author: Saideep Gona

This is the core script for the ChIP-Base application which hosts large-scale ChIP-Seq data. It allows for parameter specification and generation of binary
tf-gene binding tables. Built using the python Flask framework.
'''
import os, sys
from datetime import datetime
from multiprocessing import Process, Queue
import pickle
import time
import glob

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

from flask import Flask
from flask import render_template, send_file, redirect, url_for
from flask_wtf import FlaskForm
from flask_sqlalchemy import SQLAlchemy
from wtforms import StringField, PasswordField, BooleanField, SubmitField, FloatField, IntegerField, SelectField, SelectMultipleField
from wtforms.validators import DataRequired, NumberRange, Length
from wtforms.fields.html5 import IntegerRangeField 
from wtforms_html5 import AutoAttrMeta

import multiprocessing

#TODO Replace data with GTRD data
#TODO Set up reusable aliasing for TF names, tissue names
#TODO Figure out load handling before wide adoption
#TODO Finish motif mapping!
    # Reorganize motif mapping occurences to be bed files
    # Modify the motif filtering to match the bed file
    # Treat the motifs as peak files in terms of overlap testing
    # Generate seperate and joint count tables as output
    # TEST OUT CURRENT IMPLIMENTATION

#TODO Motif-finding from mapped peaks 
#TODO Motif finding with epigenetic priors
#TODO Standardize tissue names
#TODO Calculate more complex distributions of peak p-values
#TODO Start writing formal writeup for BioArXive
#TODO Collect more motifs from hocomoco, etc.
#TODO Subset peaks based on the presence of an overlapping motif(should it occur near summit?)
#TODO Loading bar for query
#TODO Split application.py into app_modules
#TODO Reorganize python code
#TODO Get ansible deployment running
#TODO Set up a test server with ansible deployment
#TODO Parallelize IO to improve speed
#TODO Include aliases in gene target lists
#TODO Get email sending working
#TODO Incorporate TRRUST interactions
#TODO Collect data for mouse
#TODO Improve searchability via google

# DONE?
#TODO Convert hg19 studies and update input peaks + metadata

# OTHER COMMENTS
# Similar work includes TRANSFAC, iRegulon, GTRD

app = Flask(__name__)

app.config['SECRET_KEY'] = 'sai-key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///ChIP_Base.db'

pwd = os.getcwd()
num_cpus = multiprocessing.cpu_count()
# url_root = "http://ec2-54-145-225-122.compute-1.amazonaws.com"
url_root = "chipio.cs.cmu.edu"
current_stats = {}
version = "v1.0.1"

all_genes_path = pwd + "/static_lists/all_genes.txt"
all_tfs_path = pwd + "/static_lists/all_tfs_list.txt"

en_tiss_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
en_tiss_files = glob.glob(en_tiss_dir+"*")
enhancer_tissues_list = [x.split("/")[-1][:-4] for x in en_tiss_files]
enhancer_tissues = ", ".join(enhancer_tissues_list)
enhancer_tissues = enhancer_tissues.replace("_", " ")

print("Enhancer-supported tissues: " + enhancer_tissues)

all_genes = []
with open(all_genes_path, "r") as ag:
    for line in ag:
        all_genes.append(line.rstrip("\n"))
current_stats["all_genes"]=str(len(all_genes))
print("NUMBER OF GENES: ", len(all_genes))

all_tfs = []
with open(all_tfs_path, "r") as at:
    for line in at:
        all_tfs.append(line.rstrip("\n"))
current_stats["all_tfs"]=str(len(all_tfs))
print("NUMBER OF TFS: ", len(all_tfs))

metadata_path = pwd + "/pass_metadata/metadata.pkl"     
metadata_dict = pickle.load(open(metadata_path, "rb"))      # Contains metadata on all datasets

all_peaks = pwd+"/all_peaks_small.tsv"
all_motif_occs = pwd + "/all_motif_occs.tsv"

# New params

# Database Specs

db = SQLAlchemy(app)

meta_column_list = ["experiment_accession",
            "transcription_factors",
            "tissue_types"
            ]

class ChIP_Meta(db.Model):
    '''
    Stores metadata on ChIP-Seq Studies
    '''
    id = db.Column(db.Integer, primary_key=True)
    experiment_accession = db.Column(db.String(50), unique=False, nullable=False)
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return ('<Experiment Accession {}>'.format(self.experiment_accession)
        + '<Transcription Factor {}>'.format(self.transcription_factors)
        + '<Tissue Type {}>'.format(self.tissue_types) 
        )

# The order of this list determines write order 
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
            "experiment_accession"
            ]

motif_occs_column_list = [
            "chrom",
            "start",
            "end",
            "id",
            "length",
            "score",
            "log_p",
            "transcription_factors",
            "tissue_types"
]

class Peaks(db.Model):
    '''
    This stores all the peaks that have been collected, and uses the experiment accession # to relate them to each other
    '''

    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String, nullable=False)
    start = db.Column(db.Integer, nullable=False)
    end = db.Column(db.Integer, nullable=False)
    length = db.Column(db.Integer, nullable=False)
    summit = db.Column(db.Integer, nullable=False)
    pileup = db.Column(db.Float, nullable=False) 
    log_p = db.Column(db.Float, nullable=False)
    fold_enrichment = db.Column(db.Float, nullable=False)
    log_q = db.Column(db.Float, nullable=False)
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)
    experiment_accession = db.Column(db.String, nullable=False)

    def __repr__(self):
        return ('<Experiment Accession {}>'.format(self.experiment_accession)
        + '<Transcription Factor {}>'.format(self.transcription_factors)
        + '<Tissue Type {}>'.format(self.tissue_types) 
        + '<Chrom {}>'.format(self.chrom) 
        + '<Start {}>'.format(self.start) 
        + '<end {}>'.format(self.end) 
        )

class DNAaseFootprints(db.Model):
    '''
    This stores tissue-specific footprints for tighter validation of tf-binding
    '''
    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String, nullable=False)
    start = db.Column(db.Integer, nullable=False)
    end = db.Column(db.Integer, nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)

class Footprints(db.Model):

    id = db.Column(db.Integer, primary_key=True)
    chrom = db.Column(db.String, nullable=False)
    start = db.Column(db.Integer, nullable=False)
    end = db.Column(db.Integer, nullable=False) 
    transcription_factors = db.Column(db.String(30), nullable=False)
    tissue_types = db.Column(db.String(50), nullable=False)

class Query_History(db.Model):
    '''
    Stores a history of queries and query data
    '''
    id = db.Column(db.Integer, primary_key=True)
    promoter = db.Column(db.Boolean)
    enrichment = db.Column(db.Boolean)

    transcription_factors = db.Column(db.String)
    tissue_types = db.Column(db.String)

    log_p = db.Column(db.Float)
    fold_enrichment = db.Column(db.Float)

    distance_from_TSS_upstream = db.Column(db.Integer)
    distance_from_TSS_downstream = db.Column(db.Integer)

    email = db.Column(db.String)

    date_posted = db.Column(db.DateTime, nullable = False, default = datetime.utcnow)

class Presets(db.Model):
    '''
    Store a set of preset TF matrices for rapid download
    '''
    id = db.Column(db.Integer, primary_key=True)
    log_p = db.Column(db.Float)
    fold_enrichment = db.Column(db.Float)
    promoter = db.Column(db.Boolean)
    enrichment = db.Column(db.Boolean)
    table_file_path = db.Column(db.String) 

class DownloadFiles():

    def __init__(self):
        self.make_time = datetime.utcnow


    def collect_presets(self, preset_dir):

        # Format for preset filenames should be: field1:value,field2:value2,... e.g. tissues=All,tfs=All  
        def parse_f_name(fname):
            all_fields = fname.split(',')
            mapping = {}
            for field in all_fields:
                keyval = field.split("=")
                mapping[keyval[0]] = keyval[1]
            return mapping

        preset_files = glob.glob(pwd + "/" + preset_dir + "/*")
        unpath_preset_files = ["#".join(x.split("/")) for x in preset_files]
        preset_files_strip = [x.split("/")[-1] for x in preset_files]
        self.preset_files = preset_files
        self.unpath_preset_files = unpath_preset_files
        self.preset_files_strip = preset_files_strip
        self.num_preset_files = len(preset_files)
        self.preset_mappings = [parse_f_name(x) for x in preset_files_strip]
        # print(self.preset_files_strip)
        # print(self.preset_mappings)

    def collect_peaks(self, peak_dir):
        peak_files = glob.glob(pwd + "/" + peak_dir + "/*")
        unpath_peak_files = ["#".join(x.split("/")) for x in peak_files]
        peak_files_strip = [x.split("/")[-1] for x in peak_files]
        accessions = [x.rstrip("_peaks.xls") for x in peak_files_strip]
        self.peak_files = peak_files
        self.unpath_peak_files = unpath_peak_files
        self.peak_files_strip = peak_files_strip
        self.num_peak_files = len(peak_files)
        self.tfs = [metadata_dict[x]["tf"] for x in accessions]
        self.tissues = [metadata_dict[x]["tissue"] for x in accessions]


all_possible = {
"transcription_factors": list(set([x.transcription_factors for x in ChIP_Meta.query.all()])),
"tissue_types": list(set([x.tissue_types for x in ChIP_Meta.query.all()]))
}

cur_tfs_path = pwd + "/static_lists/all_tfs_cur.txt"
cur_tiss_path = pwd + "/static_lists/all_tissues_cur.txt"

print("NUMBER OF TFS FROM ChIP-SEQ STUDIES: " + str(len(all_possible["transcription_factors"])))
current_stats["study_tfs"]=str(len(all_possible["transcription_factors"]))
with open(cur_tfs_path, "w") as cur_tfs:
    for tf in all_possible["transcription_factors"]:
        cur_tfs.write(tf + "\n")

with open(cur_tiss_path, "w") as cur_ts:
    for tiss in all_possible["tissue_types"]:
        cur_ts.write(tf + "\n")

# Forms

class ParameterForm(FlaskForm):
    # log_p = StringField('-log(p) Value', validators=[DataRequired()])
    # fold_enrichment = StringField('Fold Enrichment', validators=[DataRequired()])
    # promoter_bool = BooleanField('Promoters')
    # enhancer_bool = BooleanField('Enhancers')

    promoter_bool = BooleanField('Promoters')
    enhancer_bool = BooleanField('Enhancers')

    tfs = [(x,x) for x in all_possible["transcription_factors"]]
    tissues = [(x,x) for x in all_possible["tissue_types"]]

    transcription_factors = StringField('Transcription Factors')                  # Blank form format
    tissue_types = StringField('Tissue Types')
    
    # transcription_factors = SelectField('Transcription Factors', choices=tfs)                  # Dropdown format
    # tissue_types = SelectField('Tissue Types', choices=tissues)

    pileup = FloatField('Pileup', validators=[DataRequired(message="pileup not right")])
    log_p = FloatField('-log(p) Value', validators=[DataRequired(message="logp not right")])
    fold_enrichment = FloatField('Fold Enrichment', validators=[DataRequired(message="enrichment not right")])
    log_q = FloatField('-log(q) Value', validators=[DataRequired(message="logq not right")])


    # distance_from_TSS = IntegerField('Distance from TSS', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_upstream = IntegerField('Distance from TSS Upstream', validators=[NumberRange(0, 10000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_downstream = IntegerField('Distance from TSS Downstream', validators=[NumberRange(0, 1000, message="Must be an integer in range [0,100000]")])
    peak_count = IntegerField('Regulatory Peak Count', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])

    include_motif_sites = BooleanField("Motif Inclusion")

    motif_discovery = BooleanField('Motif Discovery')

    motif_p_val = FloatField('-log(p) Value', validators=[DataRequired(message="logp for motifs not right")])
    motif_score = FloatField('Motif Score', validators=[DataRequired(message="Motif score not right")])

    email = StringField("Email", validators=[DataRequired(message="email not right")])

    submit = SubmitField('Submit')

# General Functions

def run_pipeline(user_params):

    step_num = 1
    time_string = user_params["time"]
    print("QUERY TIME STRING: "+time_string)

    print("BEGIN PIPELINE *************************************************************************************************************")
    print(user_params)

    tf_set = set(user_params["transcription_factors"])
    tissue_set = set(user_params["tissue_types"])

    removable_junk = []                         # Store disposable file paths here

    # All regulatory regions columns: chr   start   end    gene    region_type
    all_reg_regions = pwd + "/intermediates/" + time_string + "_allregregions.bed"         # Place all regulatory regions for current query here
    removable_junk.append(all_reg_regions)
    with open(all_reg_regions, "w") as all_reg:
        print(all_reg_regions, " all regulatory regions stored here")

    # Perform filtering of ChIP-Seq peaks
    temp_peaks_file = pwd + "/intermediates/" + "temppeaks_" + time_string + ".bed"
    removable_junk.append(temp_peaks_file) 
    filter_peaks(peaks_column_list, user_params, all_peaks, temp_peaks_file, time_string, removable_junk)  
    sort_in_place(temp_peaks_file)

    # Perform similar filtering of motif occurences
    temp_motif_occs_file = pwd + "/intermediates/tempmotifs_" + time_string + ".bed"
    removable_junk.append(temp_motif_occs_file)
    filter_motif_occs(motif_occs_column_list, user_params, all_motif_occs, temp_motif_occs_file, time_string, removable_junk)
    sort_in_place(temp_motif_occs_file)

    step_num = make_check(step_num, time_string, removable_junk)
    print("||||||||||||||||||Peaks Queried and Filtered")

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only

        # Creates promoter regions bed file
        inter_promoter = pwd + "/intermediates/promoters_" + time_string + ".bed"
        removable_junk.append(inter_promoter)
        create_promoters(user_params, inter_promoter, all_reg_regions)
        os.system("sort -k 1,1 -k2,2n -o " + inter_promoter + " " +inter_promoter)
        print("promoters created")     

        # Calculates peak-promoter intersection
        promoter_intersect = pwd + "/intermediates/" + "promintersect_" + time_string + ".bed"       # Computes intersect of peaks with promoter regions
        removable_junk.append(promoter_intersect)
        bed_intersect(inter_promoter, temp_peaks_file, promoter_intersect)
        sort_in_place(promoter_intersect)

        # Calculates motif-promoter intersection
        if user_params["include_motif_sites"]:
            promoter_motif_intersect = pwd + "/intermediates/" + "prommotifintersect_" + time_string + ".bed"
            removable_junk.append(promoter_motif_intersect)
            bed_intersect(inter_promoter, temp_motif_occs_file, promoter_motif_intersect)
            sort_in_place(promoter_motif_intersect)       

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only

        # Calculates enhancer-peak intersection
        enhancer_bed_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
        intersect_output = pwd + "/intermediates/" + "annotation_" + time_string + ".bed"
        removable_junk.append(intersect_output)
        process_enhancers(user_params, intersect_output, enhancer_bed_dir, all_reg_regions, time_string)   # Computes tissue-specific intersect of peaks with enhancers

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer

        # Creates promoter regions bed file
        inter_promoter = pwd + "/intermediates/promoters_" + time_string + ".bed"
        removable_junk.append(inter_promoter)
        create_promoters(user_params, inter_promoter, all_reg_regions)
        os.system("sort -k 1,1 -k2,2n -o " + inter_promoter + " " +inter_promoter)
        print("promoters created")     

        # Calculates peak-promoter intersection
        promoter_intersect = pwd + "/intermediates/" + "promintersect_" + time_string + ".bed"       #Computes intersect of peaks with promoter regions
        removable_junk.append(promoter_intersect)
        bed_intersect(inter_promoter, temp_peaks_file, promoter_intersect)
        sort_in_place(promoter_intersect)

        # Calculates motif-promoter intersection
        if user_params["include_motif_sites"]:
            promoter_motif_intersect = pwd + "/intermediates/" + "prommotifintersect_" + time_string + ".bed"
            removable_junk.append(promoter_motif_intersect)
            bed_intersect(inter_promoter, temp_motif_occs_file, promoter_motif_intersect)
            sort_in_place(promoter_motif_intersect)

        # Calculates enhancer-peak intersection
        enhancer_bed_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
        enhancer_intersect = pwd + "/intermediates/" + "enhancerintersect_" + time_string + ".bed"
        removable_junk.append(enhancer_intersect)
        process_enhancers(user_params, enhancer_intersect, enhancer_bed_dir, all_reg_regions, time_string)      # Computes tissue-specific intersect of peaks with enhancers
        sort_in_place(enhancer_intersect)
        print("enhancers processed")

        # Calculates enhancer-motif intersection
        if user_params["include_motif_sites"]:
            enhancer_motif_intersect = pwd + "/intermediates/" + "enhancermotifintersect_" + time_string + ".bed"
            removable_junk.append(enhancer_motif_intersect)
            process_enhancers(user_params, enhancer_motif_intersect, enhancer_bed_dir, all_reg_regions, time_string)      # Computes tissue-specific intersect of peaks with enhancers
            sort_in_place(enhancer_motif_intersect)


    step_num = make_check(step_num, time_string, removable_junk)
    print("||||||||||||||||||Peaks and Motifs Intersected")

    final_peak_file = pwd + "/intermediates/" + time_string + "_finalpeaks.bed"    # Place all the peaks which successfully match regulatory regions here.
    removable_junk.append(final_peak_file)
    with open(final_peak_file, "w") as final:
        print(final_peak_file, "contains all final peaks")
    tg_table = create_empty_table(all_genes, all_tfs, tf_set)           # Create empty tf-gene table data structure
    motif_tg_table = create_empty_table(all_genes, all_tfs, tf_set)

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only
        promoter_tgtable_update(user_params, promoter_intersect, tg_table, final_peak_file)      # Adds promoter-mapped peaks to tf-gene table
        if user_params["include_motif_sites"]:
            promoter_motif_tgtable_update(user_params, promoter_motif_intersect, motif_tg_table)

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only
        enhancer_tgtable_update(user_params, enhancer_intersect, tg_table, final_peak_file)        # Adds enhancer-mapped peaks to tf-gene table
        if user_params["include_motif_sites"]:
            enhancer_motif_tgtable_update(user_params, enhancer_motif_intersect, motif_tg_table)

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer
        promoter_tgtable_update(user_params, promoter_intersect, tg_table, final_peak_file)      # Adds promoter-mapped peaks to tf-gene table
        if user_params["include_motif_sites"]:
            promoter_motif_tgtable_update(user_params, promoter_motif_intersect, motif_tg_table)
        enhancer_tgtable_update(user_params, enhancer_intersect, tg_table, final_peak_file)        # Adds enhancer-mapped peaks to tf-gene table
        if user_params["include_motif_sites"]:
            enhancer_motif_tgtable_update(user_params, enhancer_motif_intersect, motif_tg_table)

    step_num = make_check(step_num, time_string, removable_junk)
    print("||||||||||||||||||Annotations Parsed")

    # Output file directory for current pipeline completed files
    output_files_dir = pwd + "/intermediates/" + time_string + "_output/"
    os.mkdir(output_files_dir)

    out_readme = output_files_dir+"Readme.txt"
    os.system("cp "+pwd+"/readme_files/ChIP_IO_tg_README.txt "+out_readme)         # Copy over readme file
    build_readme(time_string, user_params, out_readme)

    # print("Starting Motif Matching")
    # motif_sites_file = output_files_dir + time_string + "_motif_matches/"                       # Perform motif finding on regulatory regions
    # removable_junk.append(motif_sites_file)
    # output_files.append(motif_sites_file)
    # motif_site_find(motif_sites_file)

    print("Writing tg table")
    tg_write_file = output_files_dir + "tgtable_" + time_string + ".tgtable"                       # Output file for tg-table
    removable_junk.append(tg_write_file)
    write_dict_tsv(tg_table, all_genes, all_tfs, tg_write_file, tf_set)

    bin_tg_write_file = output_files_dir + "tgtable_" + time_string + ".bintgtable"                       # Output file for tg-table
    removable_junk.append(bin_tg_write_file)
    binarize_tsv(tg_table, all_genes, all_tfs, bin_tg_write_file, tf_set, user_params)

    print("Writing motif tg table")
    motiftg_write_file = output_files_dir + "motiftgtable_" + time_string + ".tgtable"                       # Output file for tg-table
    removable_junk.append(motiftg_write_file)
    write_dict_tsv(motif_tg_table, all_genes, all_tfs, motiftg_write_file, tf_set)

    bin_motiftg_write_file = output_files_dir + "motiftgtable_" + time_string + ".bintgtable"                       # Output file for tg-table
    removable_junk.append(bin_motiftg_write_file)
    binarize_tsv(motif_tg_table, all_genes, all_tfs, bin_motiftg_write_file, tf_set, user_params)

    # pickle_tg_table = output_files_dir + "tgtable_" + time_string + ".pkl"                       # Pickle file for tg-table
    # removable_junk.append(pickle_tg_table)
    # with open(pickle_tg_table, "w")  as ptable:
    #     pickle.save(pickle.highest_protocol)

    print("Starting Motif Discovery")
    motifs_disc_file = output_files_dir + time_string + "_motif_discovery.txt"                     # Perform motif discovery on mapped peaks
    removable_junk.append(motifs_disc_file)
    # output_files.append(motifs_disc_file)
    motif_discovery(final_peak_file, time_string, motifs_disc_file)

    print("||||||||||||||||||Motif Matching, Discovery and TG Table Complete")

    print("Zipping Output Files")
    zipped_contents = pwd + "/intermediates/" + time_string + ".tar.gz" 
    removable_junk.append(zipped_contents)
    # print("tar -cvzf "+zipped_contents+ " -C " + pwd+"/intermediates/ " + " ".join(output_files))
    os.system("tar -cvzf "+zipped_contents+ " -C " + output_files_dir + " .")
    solo_zipped_filename = zipped_contents.split("/")[-1]
    # print("cp "+zipped_contents+" "+pwd+"/results/"+solo_zipped_filename)
    os.system("cp "+zipped_contents+" "+pwd+"/results/"+solo_zipped_filename)

    step_num = make_check(step_num, time_string, removable_junk)
    print("||||||||||||||||||Results Ready")

    # print(removable_junk)
    
    # os.system("rm -rf "+output_files_dir)
    # for f in removable_junk:
    #     try:
    #         os.remove(f)
    #     except:
    #         print("Can't dispose of junk: " + f)
    #         continue

    print("END PIPELINE *************************************************************************************************************")


def make_check(num, time_string, removal):
    check_file = pwd+"/intermediates/" + time_string + "_check_" + str(num)
    print(check_file, "Checkfile", num)
    os.system("touch " + check_file)
    removal.append(check_file)
    return num+1

def sort_in_place(f):
    temp_dir = pwd + "/tmp/"
    temp_f = f+".tmp"

    os.system("sort -k 1,1 -k2,2n -T " +temp_dir+" "+f+" > "+temp_f)
    os.system("rm "+f)
    os.system("mv "+temp_f+" "+f)
 
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

def build_readme(time_string, user_params, readme_file):
    write_string = "\n\n" + "time stamp: "+time_string + "\n\n"
    
    write_string += "\n\n".join([
        "promoter: " + str(user_params["promoter"]), 
        "enhancer: " + str(user_params["enhancer"]),
        "transcription_factors: " + ",".join([str(x) for x in user_params["transcription_factors"]]),
        "tissue_types: " + str(",".join([str(x) for x in user_params["tissue_types"]])),
        "pileup: "+str(user_params["pileup"]),
        "log_p: "+str(user_params["log_p"]),
        "fold_enrichment: "+str(user_params["fold_enrichment"]),
        "log_q: "+str(user_params["log_q"]),
        "dist_tss_upstream: "+str(user_params["dist_tss_upstream"]),
        "dist_tss_downstream: "+str(user_params["dist_tss_downstream"]),
        "peak_count: "+str(user_params["peak_count"]),
        "include_motif_sites: "+ str(user_params["include_motif_sites"]),
        "motif_score: "+ str(user_params["motif_score"]),
        "motif_log_p: "+ str(user_params["motif_p_val"]),
        "motif_discovery: "+ str(user_params["motif_discovery"]),
        "email: "+ str(user_params["email"])
    ])

    with open(readme_file, "a") as rf:
        rf.write(write_string)

    sub_log = pwd + "/logs/submission_logs.txt"
    with open(sub_log, "a") as sl:
        sl.write(write_string +
        "\n\n******************************************************************************")

def create_promoters(user_params, filename, all_reg_regions):
    '''
    Takes in user-defined promoter definition and constructs promoters from known
    TSS sites
    '''
    chrom_lengths = {}
    with open(pwd + "/GRCh38/chrom_lengths.txt") as cl:
        for line in cl:
            p_line = line.rstrip("\n").split("\t")
            chrom_lengths[p_line[0]] = int(p_line[1])

    ref_genes = pwd + "/annotations/genes.bed"

    with open(ref_genes, "r") as genes:
        with open(filename, "w") as tmp_proms:
            with open(all_reg_regions, "a") as reg:
    
                for line in genes:
                    p_lines = line.rstrip("\n").split("\t")
                    chrom = p_lines[0]
                    rang = [int(p_lines[1]), int(p_lines[2])]
                    strand = p_lines[3]
                    g_name = p_lines[4]

                    upstream = user_params["dist_tss_upstream"]
                    downstream = user_params["dist_tss_downstream"]

                    if strand == "+":
                        tss = min(rang)
                        new_range = [tss-upstream, tss+downstream]
                        if new_range[0] < 0:
                            new_range[0] = 0
                        elif new_range[1] < 0:
                            new_range[1] = 0
                        if new_range[0] > chrom_lengths[chrom]:
                            new_range[0] = chrom_lengths[chrom]
                        elif new_range[1] > chrom_lengths[chrom]:
                            new_range[1] = chrom_lengths[chrom]

                    elif strand == "-":
                        tss = max(rang)
                        new_range = [tss-downstream, tss + upstream]
                        if new_range[0] < 0:
                            new_range[0] = 0
                        elif new_range[1] < 0:
                            new_range[1] = 0
                        if new_range[0] > chrom_lengths[chrom]:
                            new_range[0] = chrom_lengths[chrom]
                        elif new_range[1] > chrom_lengths[chrom]:
                            new_range[1] = chrom_lengths[chrom]

                    new_line = [chrom, str(new_range[0]), str(new_range[1]), strand, g_name]
                    reg_region_line = [chrom, str(new_range[0]), str(new_range[1]), g_name, "promoter"]
                    tmp_proms.write("\t".join(new_line) + "\n")
                    reg.write("\t".join(reg_region_line) + "\n")


def motif_discovery(bed_file, time_string, output_file):
    '''
    Converts a bed file into a fasta, and then performs motif discovery 
    on it using MEME.
    '''
    reference = pwd + "/GRCh38/GRCh38.p12.genome.fa"

    intermediate_file = pwd + "/intermediates/" + time_string + "_motifinter.fa"

    print("Converting bed to fasta")
    convert_bed_fasta = [
        "bedtools",
        "getfasta",
        "-fi",
        reference,
        "-bed",
        bed_file,
        "-fo",
        intermediate_file
    ]
    
    find_motifs = [
        "meme",
        intermediate_file,
        "-p",
        str(num_cpus),
        "-o",
        output_file 
    ]

    os.system(" ".join(convert_bed_fasta))

    os.system(" ".join(find_motifs))

    os.remove(intermediate_file)

def process_enhancers(user_params, intersect_out, enh_bed_dir, all_reg_regions, time_string):
    '''
    Processes enhancers on a tissue-by-tissue basis. 
    Annotation is performed individually
    by annotating peaks for a given tissue type to their 
    respective tissue-specific enhancers
    and the combined results 
    concatenated into a single file
    '''

    file_prefix = intersect_out.rstrip(".bed")
    tissue_peak_beds = []
    tissue_intersect_beds = []
    for tissue in user_params["tissue_types"]:
        # print(tissue)
        enhancer_bed = enh_bed_dir + "/" + tissue + ".bed"   # Write enhancers to regulatory region file
        if not os.path.isfile(enhancer_bed):
            print(enhancer_bed, " is not an enhancer bed file")
            continue

        temp_peaks_file = pwd + "/intermediates/annotation_"+ tissue +"_" + time_string + ".bed"   # Temporary, tissue-specific peak file
        print(temp_peaks_file)
        tissue_peak_beds.append(temp_peaks_file)                                                            
        tissue_intersect = file_prefix+"_tempintersect_"+tissue+".bed"
        tissue_intersect_beds.append(tissue_intersect)                                                      

        try:    # write enhancer regions to all regulatory regions
            with open(enhancer_bed, "r") as eb:
                with open(all_reg_regions, "a") as reg:
                    for line in eb:
                        p_line = line.rstrip("\n").split("\t")
                        reg_line = [p_line[0],p_line[1],p_line[2],p_line[6], "enhancer"]
                        reg.write("\t".join(reg_line) + "\n")
        except Exception as e:
            print("Error in process_enhancers",e)

        # enhancer_bed: 7 cols, index 6 is gene name
        # temp_peaks_file: 13 cols, index 10 is tf name (17 when intersected)
        sort_in_place(temp_peaks_file)
        bed_intersect(enhancer_bed, temp_peaks_file, tissue_intersect)
        sort_in_place(tissue_intersect)
    
    intersect_files = " ".join(tissue_intersect_beds)
    if len(intersect_files) == 0:
        os.system("touch "+intersect_out)
    else:
        os.system("cat "+ intersect_files+ " > "+intersect_out)

    for f in tissue_intersect_beds:
        try:
            os.remove(f)
        except:
            print("Can't dispose of intersect file: " + f)
            continue

def write_dict_tsv(tg_table, all_genes, all_tfs, table_write, tf_set):
    '''
    Writes dictionary table to a tsv file
    '''
    with open(table_write, "a") as tw:

        cur_tfs = [x for x in all_tfs if x in tf_set]

        tw.write("Genes\t" + "\t".join(cur_tfs) + "\n")

        for gene in all_genes:
            # if gene not in tg_table:
            #     continue
            cur_gene_string = gene + "\t"
            for tf in all_tfs:
                if tf in tf_set:
                    cur_gene_string += (str(tg_table[gene][tf]) + "\t")
            tw.write(cur_gene_string.rstrip("\t")+"\n")

def binarize_tsv(tg_table, all_genes, all_tfs, table_write, tf_set, user_params):
    '''
    Writes dictionary table to a tsv file but in binary form based on supplied peak count
    '''
    with open(table_write, "a") as tw:

        cur_tfs = [x for x in all_tfs if x in tf_set]

        tw.write("Genes\t" + "\t".join(cur_tfs) + "\n")

        for gene in all_genes:
            # if gene not in tg_table:
            #     continue
            cur_gene_string = gene + "\t"
            for tf in all_tfs:
                if tf in tf_set:
                    if tg_table[gene][tf] >= user_params["peak_count"] :
                        cur_gene_string += ("1" + "\t")
                    else:
                        cur_gene_string += ("0" + "\t")

            tw.write(cur_gene_string.rstrip("\t")+"\n")

def create_empty_table(gene_list, tf_list, tf_set):
    '''
    Creates an empty table (dictionary[tfs]->dictionary[genes]) initializing all values
    to 0 based upon a list of genes and tfs
    '''
    table = {}
    for gene in gene_list:
        tf_dict = {}
        for tf in tf_list:
            if tf in tf_set: 
                tf_dict[tf] = 0
        table[gene] = tf_dict
    
    return table

def promoter_tgtable_update(user_params, anno_file, tf_gene_table, full_peak_bed):
    '''
    Updates tf-gene table with promoter annotation results which pass constraints
    '''
    entry_count = 0
    line_count = 0
    with open(full_peak_bed, "a") as full:
        with open(anno_file, "r") as anno:
            for line in anno:
                p_l= line.rstrip("\n").split("\t")
                if line_count == 0:
                    line_count += 1
                    continue

                peak_bedline_l = p_l[5:8]                           # Write regualtory peaks to to cumulative regulatory peak bedfile   
                peak_bedline = "\t".join(peak_bedline_l) + "\n"
                full.write(peak_bedline)

                peak_id = int(p_l[8])
                gene_id = p_l[4]
                ori_peak_tf = p_l[15]    

                if gene_id in tf_gene_table:
                    if ori_peak_tf in tf_gene_table[gene_id]:
                        entry_count += 1
                        tf_gene_table[gene_id][ori_peak_tf] += 1
                else:
                    print("Did not pass Check 1 ", gene_id)
                line_count += 1
    print("Sum of promoter entries in output table: " + str(entry_count))

def enhancer_tgtable_update(user_params, anno_file, tf_gene_table, full_peak_bed):
    '''
    Updates tf-gene table with enhancer annotation results which pass constraints
    '''

    # Example line from anno_file: chrX	12974102	12974429	88841	chrX	12809489	PRPS2	chrX	12974102	12974429	70528	328	12974406	22.35	6.99766	3.33623	3.13225	eGFP-PYGO2	blood	ENCSR410DWC

    print("parsing enhancer")
    print("anno", anno_file)
    line_count = 0
    enh_add_count = 0
    with open(full_peak_bed, "a") as full:
        with open(anno_file, "r") as anno:
            for line in anno:
                p_l= line.rstrip("\n").split("\t")
                if line_count == 0:
                    line_count += 1
                    continue

                peak_bedline_l = p_l[7:10]                          # Write regualtory peaks to to cumulative regulatory peak bedfile       
                peak_bedline = "\t".join(peak_bedline_l) + "\n"
                full.write(peak_bedline)

                peak_id = int(p_l[10])
                gene_id = p_l[6]
                ori_peak_tf = p_l[17]

                if gene_id in tf_gene_table:
                    if ori_peak_tf in tf_gene_table[gene_id]:
                        enh_add_count += 1
                        tf_gene_table[gene_id][ori_peak_tf] += 1
                else:
                    print("Did not pass Check 1 ", gene_id)
                line_count += 1
    print("enhancer add count: ", enh_add_count)

def promoter_motif_tgtable_update(user_params, anno_file, tf_gene_table):
    '''
    Updates tf-gene table with promoter annotation results which pass constraints
    '''

    # Example line from anno_file: chrX	12974102	12974429	88841	chrX	12809489	PRPS2	chrX	12974102	12974429	70528	328	12974406	22.35	6.99766	3.33623	3.13225	eGFP-PYGO2	blood	ENCSR410DWC

    entry_count = 0
    line_count = 0
    with open(anno_file, "r") as anno:
        for line in anno:
            p_l= line.rstrip("\n").split("\t")
            if line_count == 0:
                line_count += 1
                continue

            peak_id = int(p_l[8])
            gene_id = p_l[4]
            ori_peak_tf = p_l[15]    

            if gene_id in tf_gene_table:
                if ori_peak_tf in tf_gene_table[gene_id]:
                    entry_count += 1
                    tf_gene_table[gene_id][ori_peak_tf] += 1
            else:
                print("Did not pass Check 1 ", gene_id)
            line_count += 1
    print("Sum of motif promoter entries in output table: " + str(entry_count))

def enhancer_motif_tgtable_update(user_params, anno_file, tf_gene_table):
    '''
    Updates tf-gene table with enhancer annotation results which pass constraints
    '''

    # Example line from anno_file: chrX	12974102	12974429	88841	chrX	12809489	PRPS2	chrX	12974102	12974429	70528	328	12974406	22.35	6.99766	3.33623	3.13225	eGFP-PYGO2	blood	ENCSR410DWC

    print("parsing enhancer")
    print("anno", anno_file)
    line_count = 0
    enh_add_count = 0
    with open(anno_file, "r") as anno:
        for line in anno:
            p_l= line.rstrip("\n").split("\t")
            if line_count == 0:
                line_count += 1
                continue

            peak_id = int(p_l[10])
            gene_id = p_l[6]
            ori_peak_tf = p_l[17]

            if gene_id in tf_gene_table:
                if ori_peak_tf in tf_gene_table[gene_id]:
                    enh_add_count += 1
                    tf_gene_table[gene_id][ori_peak_tf] += 1
            else:
                print("Did not pass Check 1 ", gene_id)
            line_count += 1
    print("motif enhancer add count: ", enh_add_count)

def constraints_met(data, user_params, constraints_type):
    '''
    Given some data in list form(a line in a bed file), the user input parameters, and the "type"
    of constraints corresponding to the given data, assesses whether the data "passes"
    '''
    if constraints_type == "peaks":
        # print(data)
        # print(user_params)
        if (float(data[6]) > user_params["pileup"] and
            float(data[7]) > user_params["log_p"] and
            float(data[8]) > user_params["fold_enrichment"] and
            float(data[9]) > user_params["log_q"]
        ):
            return True
        else:
            return False
    
    elif constraints_type == "motif_occs":
        print(data, "to be filtered in -constraints met-")
        print(user_params["motif_score"], " ", user_params["motif_p_val"])
        if (float(data[6]) > user_params["motif_score"] and
            float(data[5]) > user_params["motif_p_val"]
        ):
            print("Passed")
            print(data[6], " ", user_params["motif_score"])
            print(data[5], " ", user_params["motif_p_val"])
            return True
        else:
            return False

    elif constraints_type == "annotations":
        try:
            int(data[9])
        except:
            print("TSS DISTANCE NOT INTABLE: ", data[9])
            return False
        if ((int(data[9]) > (-1)*user_params["dist_tss_upstream"]) and (int(data[9]) < user_params["dist_tss_downstream"])
        ):
            print("TSS DISTANCE: " + data[9] + " < " + str(user_params["dist_tss"]))
            return True
        else:
            return False
    

def filter_peaks(columns, user_params, in_file_path, out_file_path, time_string, removal_bin):
    '''
    filters the full set of peaks based on the user provided parameters and 
    pre-computed motif occurences
    '''

    def write_tissue(tissue_file, line):
        if os.path.exists(tissue_file):
            with open(tissue_file, "a") as tf:
                tf.write(line)
        else:
            os.system("touch "+tissue_file)
            removal_bin.append(tissue_file)
            to_be_sorted.append(tissue_file)
            with open(tissue_file, "a") as tf:
                tf.write(line)

    to_be_sorted = []

    if not os.path.exists(out_file_path):
        os.system("touch "+out_file_path)
        removal_bin.append(out_file_path)
        to_be_sorted.append(out_file_path)

    # write_dict = {
    #     "experiment_accession": p_f.rstrip("_peaks.xls").split("/")[-1],
    #     "tissue_types": tissue,
    #     "transcription_factors": metadata_dict_ref[exp_acc][1],
    #     "chrom": p_l[0],
    #     "start": str(int(p_l[1]) - 1),
    #     "end": str(int(p_l[2]) - 1),
    #     "length": p_l[3],
    #     "summit": str(int(p_l[4]) - 1),
    #     "pileup": p_l[5],
    #     "log_p": p_l[6],
    #     "fold_enrichment": p_l[7],
    #     "log_q": p_l[8]
    # }

    num_pass_peaks = 0
    with open(in_file_path, "r") as infile:
        with open(out_file_path, "a") as out:
            for line in infile:
                p_line = line.rstrip("\n").split("\t")
                # print(p_line)
                if constraints_met(p_line, user_params, "peaks"):
                    write_tissue(pwd + "/intermediates/annotation_"+ p_line[11] +"_" + time_string + ".bed", "\t".join(p_line)+"\n")
                    out.write("\t".join(p_line)+"\n")
                    num_pass_peaks += 1

    print("number of passed peaks: "+str(num_pass_peaks))
    for f in to_be_sorted:
        sort_in_place(f)

def filter_motif_occs(columns, user_params, in_file_path, out_file_path, time_string, removal_bin):
    '''
    filters the full set of motif occurences based on the user provided parameters and 
    pre-computed motif occurences
    '''

    def write_tissue(tissue_file, line):
        if os.path.exists(tissue_file):
            with open(tissue_file, "a") as tf:
                tf.write(line)
        else:
            os.system("touch "+tissue_file)
            removal_bin.append(tissue_file)
            to_be_sorted.append(tissue_file)
            with open(tissue_file, "a") as tf:
                tf.write(line)

    to_be_sorted = []

    if not os.path.exists(out_file_path):
        os.system("touch "+out_file_path)
        removal_bin.append(out_file_path)
        to_be_sorted.append(out_file_path)

        # write_dict = {
        #     "tissue_types": tissue,
        #     "transcription_factors": tf,
        #     "chrom": p_l[0],
        #     "start": str(int(p_l[1]) - 1),
        #     "end": str(int(p_l[2]) - 1),
        #     "length": p_l[3],
        #     "score": p_l[5]
        #     "log_p": p_l[6],
        #     "id": motif_id_count
        # }

    num_pass_motifs = 0
    with open(in_file_path, "r") as infile:
        with open(out_file_path, "a") as out:
            for line in infile:
                p_line = line.rstrip("\n").split("\t")
                # print(p_line)   
                if constraints_met(p_line, user_params, "motif_occs"):
                    write_tissue(pwd + "/intermediates/annotation_"+ p_line[8] +"_" + time_string + ".bed", "\t".join(p_line)+"\n")
                    out.write("\t".join(p_line)+"\n")
                    num_pass_motifs += 1

    print("number of passed motifs: "+str(num_pass_motifs))
    for f in to_be_sorted:
        sort_in_place(f)

def convert_query_to_file(columns, query_result, user_params, file_path):
    '''
    Converts the output of a db query to a tsv file and saves the file00000
    '''

    peaks_list = []
    with open(file_path, "a") as f:
        pagesize = 5000000
        offset = 0
        row_count = Peaks.query.count()
        print("row count: ",row_count)
        while offset < row_count:


            peak_subset = (
                            Peaks.query.filter(Peaks.tissue_types.in_(user_params["tissue_types"]))
                            .filter(Peaks.transcription_factors.in_(user_params["transcription_factors"]))
                            .slice(offset, offset+pagesize).limit(pagesize)
                            )

            # new_query = query_result.slice(offset, offset+pagesize).limit(pagesize)   
            offset += pagesize
            print(offset)
            for item in new_query:
                # print(item, "item")
                outputs = [str(getattr(item,x)) for x in columns]
                if constraints_met(outputs, user_params, "peaks"):
                    peaks_list.append(outputs)
                    writeable = "\t".join(outputs) + "\n"
                    # print(writeable, "WRITE LINE")
                    f.write(writeable)

    # return peaks_list
    
def parse_input(in_string, delim, field, all_possible):
    '''
    Takes in a string and delimeter of choice and parses into a list. Handles "ALL" condition.
    '''

    if in_string.strip() == "ALL":
        return all_possible[field]
    else:
        return(in_string.strip().split(delim))

def build_query_hist(form):

    query_data = Query_History(
        log_p = form.log_p.data, 
        fold_enrichment = form.fold_enrichment.data, 
        promoter = False, 
        enrichment = True)

    return query_data

def send_mail(send_from, send_to, subject, text, file_path, file_name, server, email_user, email_password):

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Subject'] = subject

    body = text
    msg.attach(MIMEText(body,'plain'))

    # msg.attach(part)
    text = msg.as_string()  
    server = smtplib.SMTP(server,587)
    server.starttls()
    server.login(email_user,email_password)


    server.sendmail(email_user,send_to,text)
    server.quit()

# View Functions

@app.route('/')
def home():
    return render_template('home.html', current_stats=current_stats, version=version)

@app.route('/construct_query')
def construct_query():
    return render_template('construct_query.html', enhancer_tissues = enhancer_tissues)

@app.route('/promoter_form', methods=['GET', 'POST'])
def promoter_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 3, fold_enrichment = 1, log_q = 3,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
                        include_motif_sites=True,
                        motif_discovery=False,
                        motif_p_val = 3,
                        motif_score = 0.1,
                        email="send.results.here@peaks.com")
    if form.validate_on_submit():
        # print("valid")
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": True, 
            "enhancer": False,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "pileup": form.pileup.data,
            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,
            "log_q": form.log_q.data,

            "dist_tss_upstream": form.distance_from_TSS_upstream.data,
            "dist_tss_downstream": form.distance_from_TSS_downstream.data,
            "peak_count": form.peak_count.data,

            "include_motif_sites": form.include_motif_sites.data,
            "motif_discovery": form.motif_discovery.data,

            "motif_score": form.motif_score.data,
            "motif_p_val": form.motif_p_val.data,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html', time=query_data_dict["time"])
    else:
        print(form.errors)

    return render_template('promoter_form.html', form = form, tissues=all_possible["tissue_types"], tfs=all_possible["transcription_factors"])

@app.route('/enhancer_form', methods=['GET', 'POST'])
def enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 3, fold_enrichment = 1, log_q = 3,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
                        include_motif_sites=True,
                        motif_discovery=False,
                        motif_p_val = 3,
                        motif_score = 0.1,
                        email="send.results.here@peaks.com")
    if form.validate_on_submit():
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": False, 
            "enhancer": True,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "pileup": form.pileup.data,
            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,
            "log_q": form.log_q.data,

            "dist_tss_upstream": 1,
            "dist_tss_downstream": 1,
            "peak_count": form.peak_count.data,

            "include_motif_sites": form.include_motif_sites.data,
            "motif_discovery": form.motif_discovery.data,

            "motif_score": form.motif_score.data,
            "motif_p_val": form.motif_p_val.data,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html', time=query_data_dict["time"])

    else:
        print(form.errors)
    return render_template('enhancer_form.html', form = form, tissues=all_possible["tissue_types"], tfs=all_possible["transcription_factors"])

@app.route('/promoter_enhancer_form', methods=['GET', 'POST'])
def promoter_enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 3, fold_enrichment = 1, log_q = 3,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
                        include_motif_sites=True,
                        motif_discovery=False,
                        motif_p_val = 3,
                        motif_score = 0.1,
                        email="send.results.here@peaks.com")
    if form.validate_on_submit():
        query_data = build_query_hist(form)
        db.session.add(query_data)
        db.session.commit()
        query_data_dict = {                         # Input data for the pipeline
            "promoter": True, 
            "enhancer": True,

            "transcription_factors": parse_input(form.transcription_factors.data, ",", "transcription_factors", all_possible),
            "tissue_types": parse_input(form.tissue_types.data, ",", "tissue_types", all_possible),

            "pileup": form.pileup.data,
            "log_p": form.log_p.data, 
            "fold_enrichment": form.fold_enrichment.data,
            "log_q": form.log_q.data,

            "dist_tss_upstream": form.distance_from_TSS_upstream.data,
            "dist_tss_downstream": form.distance_from_TSS_downstream.data,
            "peak_count": form.peak_count.data,

            "include_motif_sites": form.include_motif_sites.data,
            "motif_discovery": form.motif_discovery.data,

            "motif_score": form.motif_score.data,
            "motif_p_val": form.motif_p_val.data,

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html', time=query_data_dict["time"])

    else:
        print(form.errors)
    return render_template('promoter_enhancer_form.html', form = form)

@app.route('/output')
def output():
    return render_template('output.html', query_data = data)

@app.route('/complete')
def complete():
    return render_template('complete.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/downloads')
def downloads():
    # send_file("/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/peaks/ENCSR480LIS_peaks.xls", as_attachment=True)
    return render_template('downloads.html', download_files = download_files)

@app.route('/tools')
def tools():
    # send_file("/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/peaks/ENCSR480LIS_peaks.xls", as_attachment=True)
    return render_template('tools.html', tissues=all_possible["tissue_types"], tfs=all_possible["transcription_factors"])

@app.route("/download_file/<file_path>", methods=['GET', 'POST'])
def download_file(file_path):
    print("downloading,", file_path)
    file_path_proper = "/".join(file_path.split("#"))
    if file_path is None:
        print("path is None")
    try:
        time = "_".join(str(datetime.utcnow()).split(" "))
        down_log = pwd + "/logs/download_logs.txt"
        with open(down_log, "a") as dl:
            dl.write(time + " | " + file_path_proper +
            "\n\n******************************************************************************")
        return send_file(file_path_proper, as_attachment=True)
    except Exception as e:
        print("problem with path")
        return 
    # return render_template('downloads.html', download_files = download_files)

@app.route("/download_results/<query_id>", methods=['GET', 'POST'])
def download_results(query_id):
    pwd = os.getcwd()
    results_file = pwd + "/results/" + query_id + ".tar.gz"
    print("download page refresh")
    return render_template(
                            'results.html', 
                            results_url = url_root+url_for('download_results', query_id=query_id),
                            query_id = query_id, 
                            results_file = results_file, 
                            results_file_unpath = "#".join(results_file.split("/")), 
                            file_ready = os.path.isfile(results_file))

@app.route('/contact')
def contact():
    return render_template('contact.html')

download_files = DownloadFiles()
download_files.collect_peaks("pass_peaks")
download_files.collect_presets("presets")

# print(download_files.preset_mappings)

