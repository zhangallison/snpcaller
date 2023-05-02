# About
This script simulates a SNP caller that returns a putative genotype given a BAM file of data and putative SNPS to analyze.

# How to run the code
This code loads the libraries `Rsamtools`, `GenomicRanges`, `GenomicAlignments`, and `seqinr`, so before running the script it may be necessary to run the following commands in R console to install the packages:
```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("Rsamtools"))
install.packages("seqinr") 
```
If there are any errors regarding loading libraries, install the necessary library.

To run the code with example files `align_sort.bam` and `metadata.tsv`, run the following command in terminal:
```
Rscript snpcaller.R align_sort.bam metadata.tsv
```

## Format for `metadata.tsv`
An example `metadata.tsv` file exists in the project folder, but to create your own test file, follow the following format:

| chr | pos | ref | alt | maf |
| :--- | :--- | :--- | :--- | :--- |
| chr1:1000000-2000000 | 172741 | A | T | 0.05 |
| chr1:1000000-2000000 | 325026 | C | G | 0.17 |
| chr1:1000000-2000000 | 375797 | A | T | 0.10 |
| chr1:1000000-2000000 | 423797 | T | A | 0.04 |
| chr1:1000000-2000000 | 518726 | C | G | 0.04 |

where `chr` is the chromosome, `pos` is the position of the SNP, `ref` is the major allele, `alt` is the minor allele, and `maf` is the minor allele frequency.
