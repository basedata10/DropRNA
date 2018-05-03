# DropRNA

## Prepare configure file: config_drops.ini
The pipeline need samtools, STAR software.
The barcode whitelist files for indrop and 10X should be placed under whitelistDir.
The key process of read alignment and tagging to genes are inspired and borrowed from the open source cellranger pipeline(https://github.com/10XGenomics/cellranger). The refernces of genome index and transcriptome can be downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest.
In the config file, the directory of cellrange references is named as cellranger_<genome>.

[Drops]
samtools = /path/to/samtools
star = /path/to/STAR
whitelistDir = /path/to/whitelist_file_directory
cellranger_hg38 = /path/to/reference/refdata-cellranger-GRCh38-1.2.0/

## Extract the Cell Barcode
Counting the number of each kinds of barcode; this will genrate a barcode_count.<sample>.csv;

## Cell Barcode correction and filtering
Correcting the cell barcode with 1bp mismatch, filtering the barcode with min number of reads;

## Split the reads of valid Cell Barcodes
The raw pair-end raw reads are splitted to 16 single end, according to the 2bp prefix of barcode;
For example, we will get: split.<sample>.<AA|AT|AC|AG...|GG>.fq
This enables multiprocessing of the dataset;

## Star Alignment
4 fastq files runs at the same time;
The bam file sorted by sequence header is generated;

## Reads tagging
Tagging the reads alignment position to the corresponding gene name

## Genrating UMI table
