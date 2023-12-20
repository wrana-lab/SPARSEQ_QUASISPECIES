#/bin/bash
#quick wrapper script for quasispecies Analysis

start=`date +%s`

#get SE_fastq_files
echo "Unzipping fastq files..."
gunzip *.gz

ls *.fastq > SE_fastq_files.txt

##run first python scripts
echo "Running sequence pileup script..."
python3 sparseq_bc_srbd_quasispecies_standardized.py

#take those outputs and align them in clustal
echo "Aligning SWF26..."
clustalo -i filtered_Srbd_sequences.fa  -o filtered_Srbd_sequences.fasta --force

##finally run second python script with the aligned fasta files
echo "Running fasta cutter script..."
python3 fasta_cutter_srbd_standardized.py

#Zip
echo "Zipping fastq files..."
gzip *.fastq

end=`date +%s`

runtime=$((end-start))
echo "Runtime was $runtime seconds"
