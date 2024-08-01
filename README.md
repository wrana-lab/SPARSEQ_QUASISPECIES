SPARSEQ QUASISPECIES

This is the pipeline for quasispecies analysis.


Requirements: 
-python 3.9 with openpyxl,pandas 
-clustal omega [http://www.clustal.org/omega/] 
-R 4.2.2 with dplyr, Biostrings, openxlsx
-BC row match and column match reference files [included]
-optional: Snapgene Viewer for visualization of sequences [https://www.snapgene.com/snapgene-viewer]
-tested on MacOS 12.7.5


Initial folder structure:
-/example_run
	-/bc1
		-sample fastq files
		-quasispecies_wrapper_standardized.sh
		-sparseq_bc_srbd_quasispecies_standardized.py
		-fasta_cutter_srbd_standardized.py
		-BCcolmatch.csv
		-BCrowmatch.csv
	-/bc2
		-sample fastq files
		-quasispecies_wrapper_standardized.sh
		-sparseq_bc_srbd_quasispecies_standardized.py
		-fasta_cutter_srbd_standardized.py
		-BCcolmatch.csv
		-BCrowmatch.csv
	-/qs_aggregator.R

Script descriptions:

quasispecies_wrapper_standardized.sh Wrapper script for the QS analysis pipeline. Unzips fastq files, writes list of filenames, executes sequence pileup script [sparseq_bc_srbd_quasispecies_standardized.py], executes cluster omega locally for alignment and conversion to aligned fasta, executes fasta cut/trim script, then re-zips fastq files and reports run time. 

sparseq_bc_srbd_quasispecies_standardized.py Initial processing script for raw fastq files. It processes pairs of fastq files by checking for matching read IDs, checking that the R1 and R2 BCs match the sample's well in the 384w plate, then selecting, binning, and counting Srbd sequences. Next it uses a preset cutoff to filter the binned sequences based on percentage of total sequence counts. The cutoff is calculated previously by running the pipeline at a set of standarizedd cutoffs [1%, 0.5%, 0.2%, 0.1%, 0.05%, 0.02%, 0.01%] and then using them to form a linear model to get the x-intercept, which is used as the finalized cutoff for the dataset (please see methods of our manuscript for details). Sequences passing the cutoffs are written to a .fa file.

Clustal Omega [installed and used locally] aligns sequences that passed the cutoffs in the .fa file and convert to an aligned .fasta file, as an intermediate step in the wrapper script and as an optional step at the end of the pipeline for visualization.

fasta_cutter_srbd_standardized.py This script takes the aligned .fasta file and organizes the sequences, and outputs a file that has the base for each position for each sequence. This output is used later to determine the contents of each position in the Srbd sequence.

qs_aggregator.R This script takes one run's aligned sequences and trims, aggregates, compares between barcode copies, and translates them, outputting the sequences with a reference sequence.

As detailed in our manuscript we did extensive studies to develop customized read percent cutoffs for each unique run. In this case the cutoffs based on our studies are 0.1707 for BC1 and 0.1841 for BC2 and these are already hard coded in this example dataset. 

Instructions:
```
1. Load fastq files for each barcoded run into folders bc1 and bc2. Exclude controls. Leave the files gzipped (ie .fastq.gz).
2. Also load scripts into the folders if needed. Edit customized cutoffs if needed. 
3. On command line (terminal on Mac OS) navigate to the Runfolder containing bc1 and bc2 folders.
4. For each of the bc1 and bc2 folders, execute the wrapper script in the respective folder:
	bash quasispecies_wrapper_standardized.sh
5. When the analysis is complete for both bc1 and bc2, open the qs_aggregator script in R and set the working directory to the example_run folder.
6. Run the R script to load in the outputs from the analysis and compare between bc1 and bc2 files. This outputs any final sequences that pass filters to .fa files that also include reference sequences. 
7. Optionally the resulting fasta files can be viewed in snap gene viewer after alignment with clustalo:
	 clustalo -i srbd_topsequences.fa -o srbd_topsequences.fasta
	 clustalo -i srbd_topAAsequences.fa -o srbd_topAAsequences.fasta 
```
