from flask_app import db
from flask_app import Peaks, ChIP_Meta, Query_History

# x = Peaks.query.filter(Peaks.chrom=="chrX").all() 
# x = Query_History.query.all() 
# x = ChIP_Meta.query.filter(ChIP_Meta.tissue_types=="brain").all() 

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

def convert_query_to_file(columns, query_result, file_path):
    '''
    Converts the output of a db query to a tsv file and saves the file
    '''
    peaks_list = []
    with open(file_path, "w") as f:
        x = 0
        for item in query_result:
            outputs = [str(getattr(item,x)) for x in columns]
            peaks_list.append(outputs)
            writeable = "\t".join(outputs) + "\n"
            # print(writeable, "WRITE LINE")
            f.write(writeable)

    return peaks_list

# x = (Peaks.query.filter(Peaks.tissue_types.in_("esophagus")))

x = Peaks.query.all()
convert_query_to_file(peaks_column_list, x, "test_out.tsv")