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

# NEXT STEP IS TO HAVE EVERYTHING WORKING

#TODO PEAK COUNT PER REGION constraint
#TODO Epigenetic priors
#TODO Convert hg19 studies and update input peaks + metadata
#TODO Searchable download table
#TODO Write up general text, Home, About, and Contact
#TODO Collect more motifs from hocomoco, etc.


# Implemented, but still testing
#TODO Motif-mapping for non-ChIP-Seq TFs, build new database
#TODO Motif-finding from mapped peaks

app = Flask(__name__)

# os.system("python build_db.py")

app.config['SECRET_KEY'] = 'sai-key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///ChIP_Base.db'
# app.config['SQLALCHEMY_DATABASE_URI'] = os.environ['DATABASE_URL']      #Deployment database
# print(os.environ['DATABASE_URL'])

pwd = os.getcwd()

all_genes_path = pwd + "/static_lists/all_genes.txt"
all_tfs_path = pwd + "/static_lists/all_tfs_list.txt"

all_genes = []
with open(all_genes_path, "r") as ag:
    for line in ag:
        all_genes.append(line.rstrip("\n"))
print("NUMBER OF GENES: ", len(all_genes))

all_tfs = []
with open(all_tfs_path, "r") as at:
    for line in at:
        all_tfs.append(line.rstrip("\n"))
print("NUMBER OF TFS: ", len(all_tfs))

metadata_path = pwd + "/pass_metadata/metadata.pkl"     
metadata_dict = pickle.load(open(metadata_path, "rb"))      # Contains metadata on all datasets
# print(metadata_dict.keys())

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



# IF BUILDING TABLE, COMMENT OUT THE NEXT FEW SETS OF LINES


all_possible = {
"transcription_factors": list(set([x.transcription_factors for x in ChIP_Meta.query.all()])),
"tissue_types": list(set([x.tissue_types for x in ChIP_Meta.query.all()]))
}

cur_tfs_path = pwd + "/static_lists/all_tfs_cur.txt"
cur_tiss_path = pwd + "/static_lists/all_tissues_cur.txt"

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

    transcription_factors = StringField('Transcription Factors')
    tissue_types = StringField('Tissue Types')
    
    pileup = FloatField('Pileup', validators=[DataRequired(message="pileup not right")])
    log_p = FloatField('-log(p) Value', validators=[DataRequired(message="logp not right")])
    fold_enrichment = FloatField('Fold Enrichment', validators=[DataRequired(message="enrichment not right")])
    log_q = FloatField('-log(q) Value', validators=[DataRequired(message="logq not right")])


    # distance_from_TSS = IntegerField('Distance from TSS', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_upstream = IntegerField('Distance from TSS Upstream', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    distance_from_TSS_downstream = IntegerField('Distance from TSS Downstream', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])
    peak_count = IntegerField('Regulatory Peak Count', validators=[NumberRange(0, 100000, message="Must be an integer in range [0,100000]")])

    email = StringField("Email", validators=[DataRequired(message="email not right")])

    submit = SubmitField('Submit')

# General Functions

def run_pipeline(user_params):

    time_string = user_params["time"]

    print("BEGIN PIPELINE *************************************************************************************************************")
    # print(user_params)

    removable_junk = []

    # All regulatory regions columns: chr   start   end    gene    region_type
    all_reg_regions = pwd + "/intermediates/" + time_string + "_allregregions.bed"         # Place all regulatory regions for current query here
    removable_junk.append(all_reg_regions)
    with open(all_reg_regions, "w") as all_reg:
        print(all_reg_regions, " all regulatory regions stored here")

    peak_subset = (
                    Peaks.query.filter(Peaks.tissue_types.in_(user_params["tissue_types"]))
                    .filter(Peaks.transcription_factors.in_(user_params["transcription_factors"]))
                    )
            
    study_list = [x.experiment_accession for x in peak_subset]

    inter_promoter = pwd + "/intermediates/promoters_" + time_string + ".bed"
    removable_junk.append(inter_promoter)
    create_promoters(user_params, inter_promoter, all_reg_regions)                          # Creates custom promoter set with query vals

    temp_peaks_file = pwd + "/intermediates/" + "temppeaks_" + time_string + ".bed"
    removable_junk.append(temp_peaks_file)
    convert_query_to_file(peaks_column_list, peak_subset, user_params,  temp_peaks_file)        # Creates peak bed file    

    # os.system("sort -k 1,1 -k2,2n " + temp_peaks_file)

    print("||||||||||||||||||Peaks Queried and Promoters Created")

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only

        promoter_intersect = pwd + "/intermediates/" + "promintersect_" + time_string + ".bed"       # Computes intersect of peaks with promoter regions
        removable_junk.append(promoter_intersect)
        bed_command = [
        "bedtools",
        "intersect",
        "-a",
        inter_promoter,
        "-b",
        temp_peaks_file,
        "-wb",
        ">",
        promoter_intersect
        ]
        os.system(" ".join(bed_command))

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only

        enhancer_bed_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
        intersect_output = pwd + "/intermediates/" + "annotation_" + time_string + ".bed"
        removable_junk.append(intersect_output)
        process_enhancers_peaks(user_params, intersect_output, enhancer_bed_dir, all_reg_regions)   # Computes tissue-specific intersect of peaks with enhancers

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer

        enhancer_bed_dir = pwd + "/enhancer-gene/processed/tissue_beds/"
        intersect_output = pwd + "/intermediates/" + "annotation_" + time_string + ".bed"
        removable_junk.append(intersect_output)
        process_enhancers(user_params, intersect_output, enhancer_bed_dir, all_reg_regions)      # Computes tissue-specific intersect of peaks with enhancers

        promoter_intersect = pwd + "/intermediates/" + "promintersect_" + time_string + ".bed"       #Computes intersect of peaks with promoter regions
        removable_junk.append(promoter_intersect)
        bed_command = [
        "bedtools",
        "intersect",
        "-a",
        inter_promoter,
        "-b",
        temp_peaks_file,
        "-wb",
        ">",
        promoter_intersect
        ]
        os.system(" ".join(bed_command))

    print("||||||||||||||||||Peaks Annotated")

    final_peak_file = pwd + "/intermediates/" + time_string + "_finalpeaks.bed"    # Place all the peaks which successfully match regulatory regions here.
    removable_junk.append(final_peak_file)
    with open(final_peak_file, "w") as final:
        print(final_peak_file, "contains all final peaks")

    tg_table = create_empty_table(all_genes, all_tfs)           # Create empty tf-gene table data structure

    if user_params["promoter"] and not user_params["enhancer"]: # Promoter only
        parse_promoter(user_params, promoter_intersect, tg_table, final_peak_file)      # Adds promoter-mapped peaks to tf-gene table

    elif user_params["enhancer"] and not user_params["promoter"]:   # Enhancer only
        parse_enhancer(user_params, intersect_output, tg_table, final_peak_file)        # Adds enhancer-mapped peaks to tf-gene table

    elif user_params["promoter"] and user_params["enhancer"]:   # Promoter and Enhancer
        parse_promoter(user_params, promoter_intersect, tg_table, final_peak_file)      # Adds promoter-mapped peaks to tf-gene table
        parse_enhancer(user_params, intersect_output, tg_table, final_peak_file)        # Adds enhancer-mapped peaks to tf-gene table

    print("||||||||||||||||||Annotations Parsed")

    output_files = []

    # print("Starting Motif Matching")
    # motif_sites_file = pwd + "/intermediates/" + time_string + "_motif_matches/"                       # Perform motif finding on regulatory regions
    # removable_junk.append(motif_sites_file)
    # output_files.append(motif_sites_file)
    # motif_site_find(motif_sites_file)

    print("Writing tg table")
    tg_write_file = pwd + "/intermediates/" + "tgtable_" + time_string + ".tgtable"                       # Output file for tg-table
    removable_junk.append(tg_write_file)
    output_files.append(tg_write_file)
    write_dict_tsv(tg_table, all_genes, all_tfs, tg_write_file)

    print("Starting Motif Discovery")
    motifs_disc_file = pwd + "/intermediates/" + time_string + "_motif_discovery.txt"                     # Perform motif discovery on mapped peaks
    removable_junk.append(motifs_disc_file)
    output_files.append(motifs_disc_file)
    motif_discovery(all_reg_regions, time_string, motifs_disc_file)

    print("||||||||||||||||||Motif Matching, Discovery and TG Table Complete")

    print("Zipping Output Files")
    zipped_contents = pwd + "/intermediates/" + time_string + "_output.tar.gz" 
    os.system("tar -cvzf "+zipped_contents+ " " + " ".join(output_files))
    solo_zipped_filename = zipped_contents.split("/")[-1]

    send_from = "ChIPBaseApp@gmail.com"
    password = "chipbase"
    send_to = user_params["email"]
    subject = "Your output from ChIP-IO"
    text = "Here is your transcription factor - gene interaction table"
    server = 'smtp.gmail.com'
    send_mail(send_from, send_to, subject, text, zipped_tg_table, solo_zipped_filename, server, send_from, password)
    print(user_params["email"])

    print("||||||||||||||||||Email Sent")

    for f in removable_junk:
        try:
            os.remove(f)
        except:
            continue

    print("END PIPELINE *************************************************************************************************************")

def create_promoters(user_params, filename, all_reg_regions):

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

def motif_site_find(bed_file, time_string, output_dir, tf_set):
    '''
    Converts a bed file of regulatory regions into a fasta, and then performs motif finding 
    for each tf of interest. Outputs a file of locations and corresponding genes/tfs
    '''

    motif_mapping = [
        "fimo",
        motif_file,

    ]

def motif_discovery(bed_file, time_string, output_file):
    '''
    Converts a bed file into a fasta, and then performs motif discovery on it using MEME.
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
        "-o",
        output_file
    ]

    os.system(" ".join(convert_bed_fasta))

    os.system(" ".join(find_motifs))

    os.remove(intermediate_file)

def process_enhancers(user_params, intersect_out, bed_dir, all_reg_regions):
    '''
    Processes enhancers on a tissue-by-tissue basis. Takes in an already filtered set of
    peaks and partitions them by tissue-type. Annotation is then performed individually
    and the combined results concatenated into a single file
    '''



    file_prefix = intersect_out.rstrip(".bed")
    tissue_peak_beds = []
    tissue_intersect_beds = []
    for tissue in user_params["tissue_types"]:
        # print(tissue)
        enhancer_bed = bed_dir + "/" + tissue + ".bed"                              # Write enhancers to regulatory region file
        if not os.path.isfile(enhancer_bed):
            continue

        temp_peaks_file = file_prefix+"_temppeaks_"+tissue+".bed"
        tissue_peak_beds.append(temp_peaks_file)
        tissue_intersect = file_prefix+"_tempintersect_"+tissue+".bed"
        tissue_intersect_beds.append(tissue_intersect)

        peak_subset = (
                        Peaks.query.filter(Peaks.tissue_types.in_([tissue]))
                        .filter(Peaks.transcription_factors.in_(user_params["transcription_factors"]))
                        )

        convert_query_to_file(peaks_column_list, peak_subset, user_params, temp_peaks_file)


        try:
            with open(enhancer_bed, "r") as eb:
                with open(all_reg_regions, "a") as reg:
                    for line in eb:
                        p_line = line.rstrip("\n").split("\t")
                        reg_line = [p_line[0],p_line[1],p_line[2],p_line[6], "enhancer"]
                        reg.write("\t".join(reg_line) + "\n")
        except Exception as e:
            print(e)

        bed_command = [
        "bedtools",
        "intersect",
        "-a",
        enhancer_bed,
        "-b",
        temp_peaks_file,
        "-wb",
        ">",
        tissue_intersect
        ]
        os.system(" ".join(bed_command))
    
    intersect_files = " ".join(tissue_intersect_beds)
    os.system("cat "+ intersect_files+ " > "+intersect_out)
    
    for f in (tissue_intersect_beds + tissue_peak_beds):
        os.remove(f)


def write_dict_tsv(tg_table, all_genes, all_tfs, table_write):
    '''
    Writes dictionary table to a tsv file
    '''
    with open(table_write, "a") as tw:

        tw.write("Genes\t" + "\t".join(all_tfs) + "\n")

        for gene in all_genes:
            # if gene not in tg_table:
            #     continue
            cur_gene_string = gene + "\t"
            for tf in all_tfs:
                # if tf not in tg_table[gene]:
                #     continue
                cur_gene_string += (str(tg_table[gene][tf]) + "\t")
            tw.write(cur_gene_string.rstrip("\t")+"\n")

def create_empty_table(gene_list, tf_list):
    '''
    Creates an empty table (dictionary[tfs]->dictionary[genes]) initializing all values
    to 0 based upon a list of genes and tfs
    '''
    table = {}
    for gene in gene_list:
        tf_dict = {}
        for tf in tf_list:
            tf_dict[tf] = 0
        table[gene] = tf_dict
    
    return table

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
            # print("Peaks constraints satisfied")
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

def parse_enhancer(user_params, anno_file, tf_gene_table, full_peak_bed):
    '''
    Updates tf-gene table with enhancer annotation results which pass constraints
    '''

    # Example line from anno_file: chrX	12974102	12974429	88841	chrX	12809489	PRPS2	chrX	12974102	12974429	70528	328	12974406	22.35	6.99766	3.33623	3.13225	eGFP-PYGO2	blood	ENCSR410DWC

    print("parsing enhancer")
    line_count = 0
    with open(full_peak_bed, "a") as full:
        with open(anno_file, "r") as anno:
            # print(tf_gene_table["PAX7"].keys())
            for line in anno:
                # print(line)
                p_l= line.rstrip("\n").split("\t")
                # print(p_l, "annotation line")
                # print(line_count)
                if line_count == 0:
                    line_count += 1
                    continue

                peak_bedline_l = p_l[7:10]                          # Write regualtory peaks to to cumulative regulatory peak bedfile       
                peak_bedline = "\t".join(peak_bedline_l) + "\n"
                full.write(peak_bedline)

                peak_id = int(p_l[10])
                gene_id = p_l[6]
                ori_peak_tf = Peaks.query.filter(Peaks.id==peak_id)[0].transcription_factors

                if gene_id in tf_gene_table:

                    # print("CHECK gene")
                    if ori_peak_tf in tf_gene_table[gene_id]:
                        # print("CHECK 2, added")
                        tf_gene_table[gene_id][ori_peak_tf] += 1
                else:
                    print("Did not pass Check 1 ", gene_id)
                line_count += 1


def parse_promoter(user_params, anno_file, tf_gene_table, full_peak_bed):
    '''
    Updates tf-gene table with promoter annotation results which pass constraints
    '''
    line_count = 0
    with open(full_peak_bed, "a") as full:
        with open(anno_file, "r") as anno:
            # print(tf_gene_table["PAX7"].keys())
            for line in anno:
                p_l= line.rstrip("\n").split("\t")
                # print(p_l, "annotation line")
                # print(line_count)
                if line_count == 0:
                    line_count += 1
                    continue

                peak_bedline_l = p_l[5:8]                           # Write regualtory peaks to to cumulative regulatory peak bedfile   
                peak_bedline = "\t".join(peak_bedline_l) + "\n"
                full.write(peak_bedline)

                peak_id = int(p_l[8])
                gene_id = p_l[4]
                ori_peak_tf = Peaks.query.filter(Peaks.id==peak_id)[0].transcription_factors
                # print(ori_peak_tf, "ori peak")
                if gene_id in tf_gene_table:

                    # print("CHECK gene")
                    if ori_peak_tf in tf_gene_table[gene_id]:
                        # print("CHECK 2, added")
                        tf_gene_table[gene_id][ori_peak_tf] += 1
                else:
                    print("Did not pass Check 1 ", gene_id)
                line_count += 1


def send_mail(send_from, send_to, subject, text, file_path, file_name, server, email_user, email_password):

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = send_to
    msg['Subject'] = subject

    body = 'Hi there, sending this email from Python!'
    msg.attach(MIMEText(body,'plain'))

    filename=file_path
    attachment=open(filename,'rb')

    part = MIMEBase('application','octet-stream')
    part.set_payload((attachment).read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition',"attachment; filename= "+file_name)

    msg.attach(part)
    text = msg.as_string()
    server = smtplib.SMTP(server,587)
    server.starttls()
    server.login(email_user,email_password)


    server.sendmail(email_user,send_to,text)
    server.quit()

def convert_query_to_file(columns, query_result, user_params, file_path):
    '''
    Converts the output of a db query to a tsv file and saves the file
    '''
    peaks_list = []
    with open(file_path, "a") as f:
        x = 0
        for item in query_result:
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

    return(in_string.strip().split(delim))

def build_query_hist(form):

    query_data = Query_History(
        log_p = form.log_p.data, 
        fold_enrichment = form.fold_enrichment.data, 
        promoter = False, 
        enrichment = True)

    return query_data

# View Functions

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/construct_query')
def construct_query():
    return render_template('construct_query.html')

@app.route('/promoter_form', methods=['GET', 'POST'])
def promoter_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 1, fold_enrichment = 1, log_q = 1,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
                        email="send.results.here@peaks.com")
    if form.validate_on_submit():
        print("valid")
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

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
    return render_template('promoter_form.html', form = form, tissues=all_possible["tissue_types"], tfs=all_possible["transcription_factors"])

@app.route('/enhancer_form', methods=['GET', 'POST'])
def enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 1, fold_enrichment = 1, log_q = 1,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
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

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
    return render_template('enhancer_form.html', form = form, tissues=all_possible["tissue_types"], tfs=all_possible["transcription_factors"])

@app.route('/promoter_enhancer_form', methods=['GET', 'POST'])
def promoter_enhancer_form():
    form = ParameterForm(transcription_factors = "ALL", tissue_types = "ALL",
                        pileup = 1, log_p = 1, fold_enrichment = 1, log_q = 1,
                        distance_from_TSS_upstream=1000,
                        distance_from_TSS_downstream=100,
                        peak_count=1,
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

            "email": form.email.data,

            "time": "_".join(str(datetime.utcnow()).split(" "))
        }
        # print(query_data_dict)
        pid=os.fork()
        if pid==0:
            run_pipeline(query_data_dict)

        return render_template('complete.html')
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
        return send_file(file_path_proper, as_attachment=True)
    except Exception as e:
        print("problem with path")
    # return render_template('downloads.html', download_files = download_files)

@app.route('/contact')
def contact():
    return render_template('contact.html')

download_files = DownloadFiles()
download_files.collect_peaks("pass_peaks")
download_files.collect_presets("presets")
print(download_files.preset_mappings)

