# SPARSEQ_QUASISPECIES
Scripts and instructions for SPARSEQ Quasispecies pipeline. 



Requirements
-python 3.9 with openpyxl, pandas
-clustal omega [http://www.clustal.org/omega/]
-R 4.2.2 with dplyr, Biostrings, openxlsx 
-BC row match and column match reference files [same as in regular BC processing pipeline]

Scripts:
1. quasispecies_wrapper_standardized.sh
Wrapper script for the QS analysis pipeline. Unzips fastq files, writes list of filenames, executes sequence pileup script [sparseq_bc_srbd_quasispecies_standardized.py], executes cluster omega locally for alignment and conversion to aligned fasta, executes fasta cut/trim script, then re-zips fastq files and reports run time.  This script runs steps 2-3-4 of this workflow.

2. sparseq_bc_srbd_quasispecies_standardized.py
Initial processing script for raw fastq files. It processes pairs of fastq files by checking for matching read IDs, checking that the R1 and R2 BCs match the sample's well in the 384w plate, then selecting, binning, and counting Srbd sequences. Next it uses a preset cutoff to filter the binned sequences based on percentage of total sequence counts. The cutoff is calculated previously by running the pipeline at a set of standard cutoffs [1%, 0.5%, 0.2%, 0.1%, 0.05%, 0.02%, 0.01%] and then using them to form a linear model to get the x-intercept, which is used as the finalized cutoff for the dataset.   Sequences passing the cutoffs are written to a .fa file. 

3. Clustal Omega is installed used locally to align sequences that passed the cutoffs and convert to an aligned .fasta file.

4. fasta_cutter_srbd_standardized.py
This script takes the aligned .fasta file and organizes the sequences, and outputs a file that has the base for each position for each sequence. This output is used later to determine the contents of each position in the Srbd sequence. 

5. qs_aggregator.R
This script takes one run's aligned sequences and aggregates, compares, and translates them, outputting the sequences with a reference sequence for easy viewing in snapgene. 



