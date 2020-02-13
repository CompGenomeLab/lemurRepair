setwd("/Lemur/01data/")
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)

# finding alignments that overlap with one another: 
# use bedtools intersect:
# 
# 1) for Lemur: 
# bedtools intersect -wo -s -a lemurOverlapsHuman_short_noFilter.bed -b lemurOverlapsHuman_short_noFilter.bed > lemur_intersect

lemur_intersect <- read_delim("./lemur_intersect", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# read in the bed file generated using generateBED.R: 
lemurOverlapsHuman <- read_delim("./lemurOverlapsHuman_short_noFilter.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

lemur_intersect$percent <-  lemur_intersect$X13 / (lemur_intersect$X3- lemur_intersect$X2)

duplicates <- lemur_intersect[duplicated(lemur_intersect$X4),4]
duplicates <- duplicates[["X4"]]
lemur_instersect_filterDup <- lemur_intersect[!(lemur_intersect$X4 %in% duplicates), ]
lemurOverlapsHuman_noDup <- lemurOverlapsHuman[!(lemurOverlapsHuman$X4 %in% duplicates), ]

# new bed file for overlaps between lemur-human, only unique alignments
write.table(lemurOverlapsHuman_noDup, "./lemurOverlapsHuman_short_noFilter_noDup.bed", sep = "\t", quote = F ,row.names = F, col.names = F)

# 2) for Human: 
# bedtools intersect -wo -s -a humanOverlapsLemur_short_noFilter.bed -b humanOverlapsLemur_short_noFilter.bed > human_intersect


human_intersect <- read_delim("./human_intersect", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

humanOverlapsLemur <- read_delim("./humanOverlapsLemur_short_noFilter.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


human_intersect$percent <-  human_intersect$X13 / (human_intersect$X3- human_intersect$X2)

duplicates_hum <- human_intersect[duplicated(human_intersect$X4),4]
duplicates_hum <- duplicates_hum[["X4"]]
human_instersect_filterDup <- human_intersect[!(human_intersect$X4 %in% duplicates_hum), ]
humanOverlapsLemur_noDup <-humanOverlapsLemur[!(humanOverlapsLemur$X4 %in% duplicates_hum), ]

# new bed file for overlaps between lemur-human, only unique alignments
write.table(humanOverlapsLemur_noDup, "./humanOverlapsLemur_short_noFilter_noDup.bed", sep = "\t", quote = F ,row.names = F, col.names = F)
