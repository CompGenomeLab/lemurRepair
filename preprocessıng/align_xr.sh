#!/bin/bash


GENOME_DIR=""
BOWTIE2_IND=""
SAMPLE=""

echo "Cut adapter"
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq

echo "Align with the reference genome"
bowtie2 -p 16 -x $BOWTIE2_IND -U ${SAMPLE}_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam

echo "Convert to bam"
samtools view -q 20 -b -o ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.sam

echo "Convert to bed"
bedtools bamtobed -i ${SAMPLE}_cutadapt.bam >${SAMPLE}_cutadapt.bed

echo "Sort bed file. Use -u to remove duplicates"
sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE}_cutadapt.bed >${SAMPLE}_cutadapt_sorted.bed

echo "Count number of mapped reads"
grep -c "^" ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_readCount.txt

echo "Get the read length distribution of the aligned and deduplicated reads"
awk '{print $3-$2}' ${SAMPLE}_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}'

echo "Split into strands"
awk '{if($6=="+"){print}}' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_plus.bed
awk '{if($6=="-"){print}}' ${SAMPLE}_cutadapt_sorted.bed >${SAMPLE}_cutadapt_sorted_minus.bed

echo "Convert bed to bedgraph"
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted_plus.bed -g ${GENOME_DIR}/genome.fa.fai -bg -scale $(cat ${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print 1000000/$1}') >${SAMPLE}_cutadapt_sorted_plus.bdg
bedtools genomecov -i ${SAMPLE}_cutadapt_sorted_minus.bed -g ${GENOME_DIR}/genome.fa.fai -bg -scale $(cat ${SAMPLE}_cutadapt_sorted_readCount.txt | awk '{print -1000000/$1}') >${SAMPLE}_cutadapt_sorted_minus.bdg

echo "Convert bedgraph to bigwig"
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted_plus.bdg ${GENOME_DIR}/genome.fa.fai ${SAMPLE}_cutadapt_sorted_plus.bw
bedGraphToBigWig ${SAMPLE}_cutadapt_sorted_minus.bdg ${GENOME_DIR}/genome.fa.fai ${SAMPLE}_cutadapt_sorted_minus.bw

echo "Count read values for transcribed (TS) and nontranscribed (NTS) strands"
bedtools intersect -sorted -a ${GENOME_DIR}/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -S -F 0.5 > ${SAMPLE}_cutadapt_sorted_TScount.txt
bedtools intersect -sorted -a ${GENOME_DIR}/genes.bed -b ${SAMPLE}_cutadapt_sorted.bed -wa -c -s -F 0.5 > ${SAMPLE}_cutadapt_sorted_NTScount.txt
