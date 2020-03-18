#!/bin/bash

hsap_ref_path="/reference_genomes/human/gencode/19/genome.fa"
mmur_ref_path="/reference_genomes/mouse_lemur/ensembl/genome.fa"

/MUMmer3.23/nucmer --prefix=hsap_mmur $hsap_ref_path $mmur_ref_path
# results in hsap_mmur.delta file
# proceed with filtering the alignments for:
# Global alignment using length*identity weighted LIS (longest increasing subset). 
# For every reference-query pair, 
# leave only the alignments which form the longest mutually consistent set
/MUMmer3.23/delta-filter -g hsap_mmur.delta > hsap_mmur.filter.global

# use show-coors to display summary information such as position, 
# percent identity and so on, of each alignment, in Btab format (-B)
/MUMmer3.23/show-coords -B -rclo hsap_mmur.filter.global > hsap_mmur.filter.coords.gloabl
