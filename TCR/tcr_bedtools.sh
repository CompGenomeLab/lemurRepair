#!/bin/bash


SAMPLE=""
TSS=""
TES=""
GENOME=""


bedtools intersect -wa -wb -a ${SAMPLE}_cutadapt_sorted.bed -b ${TSS}.bed ${TES}.bed -names tss tes -f 0.50 > ${SAMPLE}_tcr.bed

bedtools shuffle -i ${TSS}.bed -g ${GENOME} > ${TSS}_shuffled.bed
bedtools shuffle -i ${TES}.bed -g ${GENOME} > ${TES}_shuffled.bed

bedtools intersect -wa -wb -a ${SAMPLE}_cutadapt_sorted.bed -b ${TSS}_shuffled.bed ${TES}_shuffled.bed -names tss tes -f 0.50 > ${SAMPLE}_shuffled_tcr.bed
