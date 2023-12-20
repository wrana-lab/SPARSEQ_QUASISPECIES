#!/usr/bin/python

##script for splicing fasta files from the quasispecies script
#after they have been aligned in clustal.

import re

fileout_Srbd = open("srbd_aligned_list.txt", "w")

#write headers
fileout_Srbd.write("sample, sequence\n")
srbd_sequences = list()

#input files are expected to have lengths in a reasonable range but script may need to be edited
#to account for more rows due to gappy alignments

#Srbd
firstcount = 0
secondcount = 1
thirdcount = 2
#fourthcount = 3

print("Processing Srbd...")
with open("filtered_Srbd_sequences.fasta") as fasta:
    alllines = fasta.readlines()
    #check num of lines in file
    number_of_lines = len(alllines)
    while firstcount < number_of_lines :
        firstline = alllines[firstcount].strip()
        secondline = alllines[secondcount].strip()
        thirdline = alllines[thirdcount].strip()
        #fourthline = alllines[fourthcount].strip()

        #remove sequence and keep sample ID and count
        firstline = re.sub(">", "", firstline)
        firstline_trimmed = re.sub("\..*$", "", firstline)
        sequence = secondline + thirdline #+ fourthline
        srbd_sequences.append(sequence)
        fileout_Srbd.write(firstline_trimmed + "," + sequence + "\n")
        firstcount = firstcount+3
        secondcount = secondcount+3
        thirdcount = thirdcount+3
        #fourthcount = fourthcount+4

fasta.close()

print("Done.")
fasta.close()
fileout_Srbd.close()
