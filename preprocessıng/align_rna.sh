#!/bin/bash


STARGENOMEDIR=""
Pair1=""
Pair2=""
SAMPLE=""

echo "Align with the reference genome"
STAR --genomeDir ${STARGENOMEDIR} --readFilesIn ${Pair1} ${Pair2} --readFilesCommand zcat --outFileNamePrefix ${SAMPLE} --runThreadN 16 outBAMsortingThreadN 16 --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory

echo "Convert to bed"
bedtools bamtobed -i ${SAMPLE}Aligned.out.bam >${SAMPLE}.bed

echo "Sort bed file. Use -u to remove duplicates"
sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE}.bed >${SAMPLE}_sorted.bed

