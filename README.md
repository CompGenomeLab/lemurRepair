# lemurRepair
Lemur Repair Study


## Directories to include

* XR-seq_pipeline
* RNA-seq_pipeline
* TCR
* Orthology
* Circos
* XR-seq_simulation


## Preprocessing RNA-seq and XR-seq reads

### XR-seq

*Required software*

    bowtie2-2.3.4.1
    samtools-1.9
    bedtools2-2.27.1
    cutadapt-1.9.1

1. Define parameters in `align_xr.sh`

    - `GENOME_DIR`: reference fasta file
    - `BOWTIE2_IND`: bowtie2 index directory
    - `SAMPLE`: full path of fastq file (without `.fastq` exstension)

2. Run `align_xr.sh`

Output files will be located at where input files are

### RNA-seq

*Required software*

    star-2.6.1
    bedtools2-2.27.1

1. Define parameters in `make_star_index.sh`

    - `RefGenome`: reference fasta file
    - `Annotation`: reference gtf file
    - `STARGENOMEDIR`: output path for star index

2. Run `make_star_index.sh`

3. Define parameters in `align_rna.sh`

    *If your samples are SingleEnd simply leave blank `Pair2`*

    - `STARGENOMEDIR`: star index from previous step
    - `Pair1`
    - `Pair2`
    - `SAMPLE`: output prefix

## XR-seq_simulation

*TODO:ART simulation. @vogulcan*

*Required software*

    go-1.14

1. Compile `filter_syn.go`

```bash
GOOS=linux go build main.go
```

2. Run `filter_syn -h` for all command line parameters

## TCR

*Required software*

    python-3.7.4
    bedtools2-2.27.1
    plotly-4.1.0

1. Run `parse_biomart.py` for both human and mouse lemur

Simple usage:
```bash
python parse_biomart.py --biomart mart_lemur.txt --out ./lemur
```

This will create two bed files for TSS and TES, if above script sample is used filenames will be `lemur_tss.bed` and `lemur_tes.bed`

2. Define parameters in `tcr_bedtools.sh`

    - `SAMPLE`: full path of aligned bed file (without `_cutadapt_sorted.bed` exstension, by default it's same as `SAMPLE` parameter in XR-seq alignment step)
    - `TSS`: TSS bed file from previous step (without `.bed` exstension)
    - `TES`: TES bed file from previous step (without `.bed` exstension)
    - `GENOME`: tab seperated file fith chromosome names and lengths

This step will produce two bed files `${SAMPLE}_tcr.bed` and `${SAMPLE}_shuffled_tcr.bed`

3. To plot results run `plot_tcr.py` after defining file paths for all human and lemur, tss and tes tcr bed files inside `plot_tcr` function