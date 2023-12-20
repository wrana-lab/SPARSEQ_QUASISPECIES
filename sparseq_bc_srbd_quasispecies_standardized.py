#!/usr/bin/python

import os, sys
import argparse
import re
import openpyxl
from openpyxl import Workbook
import pandas as pd

###This version is for R1s and R2s with BCs --
#The purpose is to get paired reads and get exact match for row and column BC from the matching pair of R1 and R2
#and run the initial QS analysis on them - sorting reads.
srbd = "ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTT" #

def main():
    srbd_var_dic = dict()
    srbd_sum = open("srbd_summary.csv", "w")
    srbd_sum.write("sample" + "," + "sequence" + "," + "count" + "\n")
    fileout_srbd = open("srbd_sequence_list.fa", "w")
    fileout_srbd.write(">Refseq_Srbd_Omicron" + "\n" + srbd + "\n")
    filtered_sequences = open("filtered_Srbd_sequences.fa", "w")
    filtered_sequences.write(">Refseq_Srbd_Omicron" + "\n" + srbd + "\n")
    seqpercenttablefile = open("seq_percents.csv", "w")
    countlist = open("countoutputlist.csv", "w")
    countlist.write("sample" + "," + "total_SRBD_count" + "\n")


    sampleID_listR1 = [] #used to check for duplicate sampleIDs
    input_file_setR1 = []
    sampleID_listR2 = [] #used to check for duplicate sampleIDs
    input_file_setR2 = []

#get R1 AND R2 as original fastq files
    for fastq in os.listdir("./"):
        #edited this to use the trimmed fastq
        if "R1_001.fastq" in fastq:
            matches = re.match(r'.*?_(.*?)_.*?', fastq)
            sampleID_listR1.append(matches.group(1))
            input_file_setR1.append(fastq)
        elif "R2_001.fastq" in fastq:
            matches = re.match(r'.*?_(.*?)_.*?', fastq)
            sampleID_listR2.append(matches.group(1))
            input_file_setR2.append(fastq)

    input_file_setR1.sort()
    input_file_setR2.sort()

    #create list of duplicate samples
    sIDR1 = set()
    sampleIDR1_dups = {x for x in sampleID_listR1 if x in sIDR1 or (sIDR1.add(x) or False)}

    #get the row and column dicts
    row_dic = {}
    BCrowtable = './BCrowmatch.csv'
    with open(BCrowtable) as rowlines:
        for line in rowlines:
            info = line.strip()
            infol = info.split(",")
            row_dic[(infol[0])] = infol[1]
    rowlines.close()

    col_dic = {}
    coltable = './BCcolmatch.csv'
    with open(coltable) as collines:
        for line in collines:
            info = line.strip()
            infol = info.split(",")
            col_dic[(infol[0])] = infol[1]
    collines.close()

    counterone = 0
    while counterone < len(input_file_setR1):
        print("Processing file ", counterone+1)

        r1file = input_file_setR1[counterone]
        r2file = input_file_setR2[counterone]
        countfilepathR1 = "./"+r1file
        countfilepathR2 = "./"+r2file

        #Using regex to parse sampleID; include well to account for duplicates e.g. 6062021_H2O-multi_S41_R2C3_23M_S41_R1_001.fastq
        matches = re.match(r'(.*?)_(.*?)_S.{1,4}_.*?_(.{2,3})_.*', r1file)
        date = matches.group(1)
        sampleID = matches.group(2)
        sampleID = re.sub("-V1-2", "", sampleID)
        well = matches.group(3)

        if sampleID in sampleIDR1_dups:
            sampleID = sampleID+"_"+well
        count_dataR1 = [sampleID]

        #open pair of fastqs
        it=3 #start at it = 3
        with open(countfilepathR1) as f1, open(countfilepathR2) as f2:

            long_read_dict_srbd = dict()
            for r1line,r2line in zip(f1,f2):
                it+=1

                if it%4 == 0:
                    #get ID line (first line, 5th line...) check if ID lines match, and remove " 1" and " 2"
                    #you do this to make sure it's a proper pair of reads
                    idlineR1 = re.sub(" 1", "", r1line)
                    idlineR2 = re.sub(" 2", "", r2line)

                #if ID lines match AND it's a 2nd, 6th etc line  get sequences on next line
                if it%4 == 1 and idlineR1 == idlineR2:
                    read1 = r1line[:5] #the 5 is not included, it is 0:4
                    read2 = r2line[:5]
                    #match reads to barcode list
                    #read 1 is row (letter) and r2 is cols (num)
                    #what if val is not in dict
                    if read1 in row_dic and read2 in col_dic:
                        givenrow = row_dic[read1]
                        givencol = col_dic[read2]
                        givenwell = givencol + givenrow

                        if givenwell == well:
                            seq = r1line.strip()

                            if len(seq)>100 :
                                if "ACCTTTTGAGA" in seq : ##srbd
                                #need to trim primer and BC
                                    seq = seq[32:] #
                                    #remove leading AAA/AA for alignment later
                                    seq = re.sub("^AAA", "A", seq)
                                    seq = re.sub("^AA", "A", seq)
                                    #hard trim
                                    seq = seq[:96]
                                    if seq not in long_read_dict_srbd : long_read_dict_srbd[seq] = 1
                                    else : long_read_dict_srbd[seq] = 1+long_read_dict_srbd[seq]

        srbd_ = list()
        srbd_c_ = 0
        long_read_t_srbd = long_read_dict_srbd.items()
        for key,val in long_read_t_srbd :
            #append to dict

            if len(key)>85:
                srbd_.append((val, key))
                srbd_c_ = srbd_c_+ val
                #srbd_c_ is total counts of sample srbd, including dropped top seq
        countlist.write(sampleID + "," +  str(srbd_c_) + "\n")
        if srbd_c_ > 32000:
            srbd_.sort(reverse=True)
            if srbd_: #drop top count sequence
                srbd_.pop(0)

                if srbd_:
                    #if there is anything in srbd_ after dropping top
                    #convert dict to pandas frame
                    seqs_merged = pd.DataFrame.from_records(srbd_)
                    #fix colnames
                    seqs_merged.columns = ['count', 'sequence']
                    #get percent within sample
                    seqs_merged['percent_count'] = seqs_merged['count'] / srbd_c_ * 100
                    rangelength = len(seqs_merged)
                    #add dummy variable
                    seqs_merged['seq_number'] = range(1,rangelength+1)
                    seqs_merged['sample'] = sampleID

                    if len(seqs_merged) > 2:
                        filtered_seqs_merged = seqs_merged[seqs_merged["percent_count"] > 0.1581]
                        filt_dict = pd.Series(filtered_seqs_merged['count'].values,index=filtered_seqs_merged['sequence']).to_dict()
                        filt_dict_items = filt_dict.items()
                        for key,val in filt_dict_items :
                            filtered_sequences.write(">" + sampleID + "_" + str(val) + "." + key + "." + "\n" + key + "\n")
                            filtered_seqs_merged.to_csv('seq_percents.csv',mode='a', header=False)


        #this second out is the complete set of sequences.
        for val,item in srbd_ :
            fileout_srbd.write(">" + sampleID + "-" + str(val) + "." + item + "." + "\n" + item + "\n")
            srbd_sum.write(sampleID + "," + item + "," + str(val) + "\n")
        counterone+=1

    #close files.
    countlist.close()
    filtered_sequences.close()
    srbd_sum.close()
    fileout_srbd.close()
    seqpercenttablefile.close()
    print("\t\t> R1 and R2 BC QS analysis now complete...")

main()
