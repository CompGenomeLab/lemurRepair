# MUMMER - nucmer alignment of lemur chromosomes to human chromosomes
# 
# mummer outputs .delta file
#
# delta-filter filters the alignments
# For each query, leave only the alignments which form the longest 
# consistent set for the query
# 
# show-coords summarizes the .delta alignment file
# this file can be read as following: 

setwd("/Lemur/01data/")
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)

hsap_mmur_filter_global <- read_delim("./hsap_mmur.filter.coords.global", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)
# rearrange and name columns: 
# qID : Lemur chromosomes
# rID: Human chromosomes
hsap_mmur_filter_global <- hsap_mmur_filter_global[, c(1,3,6:13,18,19)]
colnames(hsap_mmur_filter_global) <- c( "qID", "qLen", "rID", 
                                        "qStart", "qEnd", "rStart", "rEnd",
                                        "percentIdent", "percentSim", "alnLen",
                                        "qStrand", "rLen")
hsap_mmur_filter_global$qID <- paste("chr", hsap_mmur_filter_global$qID, sep="")

#subset for chr names
chrNames <- paste("chr",as.character(c(1:32,"X")),sep="") 

hsap_mmur_filter_global <- hsap_mmur_filter_global[hsap_mmur_filter_global$qID %in% chrNames ,]
hsap_mmur_filter_global <- hsap_mmur_filter_global[hsap_mmur_filter_global$rID %in% chrNames ,]


# sort by qID then rID then start position on the ref:

global <- hsap_mmur_filter_global[order(hsap_mmur_filter_global$qID,hsap_mmur_filter_global$rID,hsap_mmur_filter_global$qStart),]


# generate .bed files:

# filter alignments for, 

# 1) alignment lengths
# 2) percent similarity 

# plot alignment lengths--------------------------------------------------
ggplot(global[global$alnLen <= 5000,], aes(x=alnLen)) + 
  geom_histogram(binwidth = 100) + 
  geom_vline(xintercept = 300, color="red") + 
  theme_bw() +
  labs(title="Alignment lenghts, binwidth = 100bp", x="Length (bp)", y="Count")
#looking at the histogram, remove alignments smaller than 400bp.
global <- global[global$alnLen >= 400, ]

# plot percent similarity-------------------------------------------------
ggplot(global,aes(x=percentSim)) + geom_histogram(binwidth = 5)
global <- global[global$percentSim >= 80, ]

# change the strand info to +/-
global$qStrand <- ifelse(global$qStrand == "Plus", "+", "-")

# save this object to further use in plotting: 
saveRDS(global, file="./globalObj.rds")

# bed file format: 
# 1st column: chromosome
# 2nd column: start
# 3rd column: end
# 4th column: feature name
# 5th column: strand

# lemur chromosome positions overlapping human genome, 

lemur.bed <- global[,c(1,4,5)]
lemur.bed$qID <- substring(lemur.bed$qID,5)
lemur.bed$feature <- paste(global$rID,global$rStart,global$rEnd,sep = ":")
lemur.bed$score <- 0 
lemur.bed$qStrand <- global$qStrand
# bed file requires start position < end position, so for - alignments, need to change start-end positions:
lemur.bed <- lemur.bed %>% dplyr::mutate(qStart2=ifelse(qStart>qEnd, qEnd, qStart), 
                                             qEnd2=ifelse(qEnd<qStart, qStart, qEnd)) %>% 
  dplyr::select(qID, qStart2, qEnd2, feature, score, qStrand)

write.table(lemur.bed, "./lemurOverlapsHuman_short_noFilter.bed", sep = "\t", quote = F ,row.names = F, col.names = F)

# human chromosome positions overlapping lemur genome,
# 4th column: feature name (corresponding position in lemur genome)
# 5th column: strand (always + for human, due to the way of reporting by nucmer/mummer)
human.bed <- global[,c(3,6,7)]
human.bed$rID <- substring(human.bed$rID,5)
human.bed$feature <- paste(global$qID,global$qStart,global$qEnd,sep = ":")
human.bed$score <- 0
human.bed$rStrand <- "+"
write.table(human.bed, "./humanOverlapsLemur_short_noFilter.bed", sep = "\t", quote = F ,row.names = F, col.names = F)
