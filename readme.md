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
