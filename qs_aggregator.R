#script to process quasispecies pipeline outputs and aggregate sequences 
library(openxlsx)
library(dplyr)
library(Biostrings)
### Translating DNA/RNA:
#translate(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

#prepare reference sequences
srbd_omicron<-"ATCTATCAGGCCGGTAACAAACCTTGTAATGGTGTTGCAGGTTTTAATTGTTACTTTCCTTTACGATCATATAGTTTCCGACCCACTTATGGTGTT"
srbd_omicron_AA<-"IYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGV"
srbd_delta<-"ATCTATCAGGCCGGTAGCAAACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
srbd_delta_AA<-"IYQAGSKPCNGVEGFNCYFPLQSYGFQPTNGV"
srbd_wt<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTT"
srbd_wt_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGV"
srbd_alpha<-"ATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTTATGGTGTT"
srbd_alpha_AA<-"IYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGV"

#import result seqs
bc1_aligned_list<-read.csv("bc1/srbd_aligned_list.txt")
bc2_aligned_list<-read.csv("bc2/srbd_aligned_list.txt")

#drop refseq - the first row 
bc1_aligned_list<-bc1_aligned_list[2:nrow(bc1_aligned_list),]
bc2_aligned_list<-bc2_aligned_list[2:nrow(bc2_aligned_list),]

#remove gaps, remove any leading or trailing bases that will throw off alignment 
bc1_aligned_list$sampleID<-substr(bc1_aligned_list$sample, 1, 9)
bc1_aligned_list$seq_chars<-gsub("-", "", bc1_aligned_list$sequence)
bc1_aligned_list$seq_chars<-gsub("C$", "", bc1_aligned_list$seq_chars)
mean(nchar(bc1_aligned_list$seq_chars)) #check length

bc2_aligned_list$sampleID<-substr(bc2_aligned_list$sample, 1, 9)
bc2_aligned_list$seq_chars<-gsub("-", "", bc2_aligned_list$sequence)
bc2_aligned_list$seq_chars<-gsub("C$", "", bc2_aligned_list$seq_chars)
mean(nchar(bc2_aligned_list$seq_chars))

#keep prepared columns and check overlap 
bc1_aligned_list_sm<-bc1_aligned_list[,3:4]
bc2_aligned_list_sm<-bc2_aligned_list[,3:4]
table(bc1_aligned_list_sm$seq_chars %in% bc2_aligned_list_sm$seq_chars) #24 false; 102 true
table(bc2_aligned_list_sm$seq_chars %in% bc1_aligned_list_sm$seq_chars) #80 false; 111 true

#need to keep only those sequences found in both copies of the sample
##separate vers where I merge sample ID+seq to get unique total in more intuitive way 
bc1_aligned_list_merged<-bc1_aligned_list; bc1_aligned_list_merged$exp<-NULL
bc2_aligned_list_merged<-bc2_aligned_list; bc2_aligned_list_merged$exp<-NULL
bc1_aligned_list_merged$mergedsamplesequence<-paste(bc1_aligned_list_merged$sampleID, bc1_aligned_list_merged$seq_chars, sep = "_")
bc2_aligned_list_merged$mergedsamplesequence<-paste(bc2_aligned_list_merged$sampleID, bc2_aligned_list_merged$seq_chars, sep = "_")

####compile sets of common sequences#####
bc1_final<-bc1_aligned_list_merged[bc1_aligned_list_merged$mergedsamplesequence %in% bc2_aligned_list_merged$mergedsamplesequence,]
bc2_final<-bc2_aligned_list_merged[bc2_aligned_list_merged$mergedsamplesequence %in% bc1_aligned_list_merged$mergedsamplesequence,]

#initialize table of common sample / seq
compile_common_qs <- function(x,y){
  common_aligned_list_sm <- x[0,]
  for(i in unique(y$sampleID)){
    print(i)
    sub_bc1<-subset(x, x$sampleID == i)
    sub_bc2<-subset(y, y$sampleID == i)
    for(j in unique(sub_bc2$seq_chars)){
      if (any(grepl(j, sub_bc1$seq_chars, fixed=TRUE))){
        #push sample id and sequence to output table 
        common_aligned_list_sm[nrow(common_aligned_list_sm)+1,] <- c(i, j)
        print(j)
      }
    }
  }
  return(common_aligned_list_sm)
}

bc1_aligned_list_sm<-bc1_final[,3:4]
bc2_aligned_list_sm<-bc2_final[,3:4]
bcmin1<-compile_common_qs(bc1_aligned_list_sm, bc2_aligned_list_sm)

#summarize
bcmin1_agg<-bcmin1 %>% 
  group_by(seq_chars) %>% 
  summarise(num_samples = n())
bcmin1_agg<-bcmin1_agg[order(bcmin1_agg$num_samples, decreasing = T),]
bcmin1_agg<-as.data.frame(bcmin1_agg)
nrow(bcmin1_agg) 
common_aligned_list_sm<-bcmin1_agg

write.xlsx(bcmin1, file = "run27_final_sample_sequence_list.xlsx") #individual sample and sequence pairs
write.xlsx(common_aligned_list_sm, "run27_common_sampleIDsandsequences_finalized.xlsx") #sequence and counts

#convert to AAs
common_dnastring<-DNAStringSet(common_aligned_list_sm$seq_chars)
common_aas<-as.data.frame(translate(common_dnastring))
common_translated<-cbind(common_aas,common_aligned_list_sm)
colnames(common_translated)<-c("aas", "sequence", "numSamples")

###save aggregated sequences to files 
#sequences
write.table(">Refseq_srbd_omicron", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
write.table(srbd_omicron, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_delta", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_delta, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_WT", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_alpha", file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_alpha, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)

for(i in 1:nrow(common_translated)){
  seq<-common_translated[i,2]
  exp_samp<-common_translated[i,3]
  lines<-paste0(">",exp_samp,"_samples")
  write.table(lines, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "delta_srbd_topsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#AAs
write.table(">Refseq_srbd_omicron", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = FALSE, col.names = FALSE, quote = F)
write.table(srbd_omicron_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_delta", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_delta_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_WT", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(">Refseq_srbd_alpha", file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
write.table(srbd_wt_AA, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)

for(i in 1:nrow(common_translated)){
  seq<-common_translated[i,1]
  exp_samp<-common_translated[i,3]
  lines<-paste0(">",exp_samp,"_samples")
  write.table(lines, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
  write.table(seq, file = "delta_srbd_topAAsequences.fa", row.names = FALSE, append = TRUE, col.names = FALSE, quote = F)
}

#these output files would now be aligned with clustal omega and then viewed in snapgene
##end
