PLEASE READ THROUGH THIS DOCUMENT BEFORE USING THE PROVIDED data

-Thank you for using ChIP-IO-

Current Reference Genome: GRCh38
Current Annotation: gencode.v29

------------ Output Description ------------

*** .tgtable file ***

Your output folder contains a .tgtable file. This is actually a tab-seperated value(tsv) file. It contains tabular data 
with the metadata:

Rows: Gene + Gene Regions
Columns: Transcription Factors

Table entries correspond to the number of ChIP-Seq peaks for a given transcription factor mapping to regulatory regions of
a given gene. Currently these peaks come only from raw ChIP-Seq peak data, but future versions will include motif/epigenetic
prior binding site prediction to reduce false positive rates where applicable.

------------ Important Notes! ------------

- In general, "0" table entries are not intended to provide confidence for lack of TF-gene association. ChIP-Seq studies still 
cover only a subset of transcription factors, and there may always be variation across cell states which is not well
captured. Instead, focus on interpreting non-zero entries as providing experimental confidence of association over neutral.

- If multi-tissue selection was done, the final output table will contain a union of data across tissue types specified in the
ChIP-IO query. If you are interested in tf-gene association for a single tissue type, submit a new query inputting only that 
tissue. 

------------ Current Submission ------------