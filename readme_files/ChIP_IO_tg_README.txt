PLEASE READ THROUGH THIS DOCUMENT BEFORE USING THE PROVIDED DATA

-Thank you for using ChIP-IO-

Current Reference Genome: GRCh38
Current Annotation: gencode.v29
Updated README: 6.11.19 

------------ Output Description ------------

*** .tgtable file ***

Your output folder contains one or more .tgtable files. Currently two are provided, one with peaks drawn ChIP-Seq studies and 
one containing peaks drawn from motif predictions. These are actually tab-seperated value(tsv) file. They contain tabular data 
with the metadata:

Rows: Gene + Gene Regions
Columns: Transcription Factors

Table entries correspond to the number of ChIP-Seq peaks for a given transcription factor mapping to regulatory regions of
a given gene. 

*** .bintgtable file ***

Same as above, except entries are binarized as 0 or 1 based on the "peak count" threshold.

*** .pkl file ***

Python readable binarytable object. Holds the same information as the above .tsv files but can be read into python as a nested
dictionary object. Once read in, entries can be pulled as table[gene][tf]

------------ Important Notes! ------------

- In general, "0" table entries are not intended to provide confidence for lack of TF-gene association (True Negatives). ChIP-Seq 
studies still cover only a subset of transcription factors, and there may always be variation across cell states which is not well
captured. Instead, focus on interpreting non-zero entries as providing experimental confidence of association over neutral.

- If multi-tissue selection was done, the final output table will contain a union of data across tissue types specified in the
ChIP-IO query. If you are interested in tf-gene association for a single tissue type, submit a new query inputting only that 
tissue. 

------------ Current Submission ------------
