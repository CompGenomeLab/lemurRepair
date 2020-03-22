#!/bin/bash


RefGenome=""
Annotation=""
STARGENOMEDIR=""

mkdir ${STARGENOMEDIR}
STAR --runMode genomeGenerate --runThreadN 16 --sjdbGTFfile $Annotation --genomeDir ${STARGENOMEDIR} --genomeFastaFiles $RefGenome

