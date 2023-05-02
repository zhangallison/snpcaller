#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
bam <- args[1]
mtdata <- args[2]

# if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("Rsamtools"))
# install.packages("seqinr") 

suppressMessages(library(Rsamtools, quietly = TRUE))
suppressMessages(library(GenomicRanges, quietly = TRUE))
suppressMessages(library(GenomicAlignments, quietly = TRUE))
suppressMessages(library(seqinr, quietly = TRUE))


bam_file <- bam

putative_snps <- read.csv(mtdata, sep = '\t')

gal <- readGAlignments(bam_file, param=ScanBamParam(what=c("seq","qual")))


prior_prob <- function(snp) {
  # INPUT: snp with ref, alt, maf values
  # OUTPUT: vector of prior probabilities of (RR, AA, RA)
  
  ref <- snp$ref
  alt <- snp$alt
  maf <- snp$maf
  prior_rr <- (1-maf) * (1-maf)
  prior_aa <- maf * maf
  prior_ra <- 1 - (prior_rr + prior_aa)
  c(prior_rr, prior_aa, prior_ra)
}


readsAtSnp <- function(galign, snp) {
  # INPUTS
  # galign : a GAlignments object scanned with params "seq" and qual"
  # snp : a putative snp with chr and pos values
  
  # OUTPUT: data frame with column of observations and column of P(E=1)
  
  chr <- snp$chr
  pos <- snp$pos
  grange <- GRanges(chr, ranges = IRanges(start = pos, width=1))
  readsAtPos <- as.character(stackStringsFromGAlignments(galign, region=grange))
  phredq <- PhredQuality(stackStringsFromGAlignments(galign, region=grange, what="qual"))
  qual_score <- unlist(as(phredq, "IntegerList"))
  p_error <- 10^(-qual_score/10)
  res <- as.data.frame(cbind(readsAtPos, p_error))
  res$p_error <- as.numeric(res$p_error)
  colnames(res) <- c("observations", "p_error")
  res
}

post_prob <- function(snp, readsAtPos) { 
  # INPUTS
  # snp : a row in a data frame with pos, ref, alt, maf values
  # readsAtPos : a data frame with observations and corresponding error
  
  # OUTPUT: a vector of posterior probabilities of (RR, AA, RA)
  #         where RR is major allele homozygous genotype,
  #         AA is minor allele homozygous genotype,
  #         RA is heterozygous genotype
  
  chr <- snp$chr
  pos <- snp$pos
  ref <- snp$ref
  alt <- snp$alt
  maf <- snp$maf
  
  x <- readsAtPos$observations
  p_error <- readsAtPos$p_error
  
  prob_rr <- 0 # probability of RR ("RefRef")
  prob_aa <- 0
  prob_ra <- 0
  priors <- prior_prob(snp)
  prior_rr <- priors[1]
  prior_aa <- priors[2]
  prior_ra <- priors[3]
  
  for (i in seq_len(length(x))) {
    if (x[i] == ref)
      prob_rr <- prob_rr + log(1 - p_error[i])
    else
      prob_rr <- prob_rr + log(p_error[i])
  }
  
  for (i in seq_len(length(x))) {
    if (x[i] == alt)
      prob_aa <- prob_aa + log(1 - p_error[i])
    else
      prob_aa <- prob_aa + log(p_error[i])
  }
  
  for (i in seq_len(length(x))) {
    prob_ra <- prob_ra + log((0.5 * (1-p_error[i])) + (0.5 * p_error[i]))
  }
  prob_rr <- exp(prob_rr)
  prob_aa <- exp(prob_aa)
  prob_ra <- exp(prob_ra)
  
  cond_rr <- prob_rr * prior_rr
  cond_aa <- prob_aa * prior_aa
  cond_ra <- prob_ra * prior_ra
  
  post_rr <- cond_rr/(cond_rr + cond_aa + cond_ra)
  post_aa <- cond_aa/(cond_rr + cond_aa + cond_ra)
  post_ra <- cond_ra/(cond_rr + cond_aa + cond_ra)

  c(post_rr, post_aa, post_ra)
}

snps_post_prob <- function(snps, galign) {
  # INPUTS
  # snps: table of putative snps
  # galign: a GAlignments object scanned with params "seq" and qual"
  
  N <- dim(snps)[1] # number of snps
  df <- data.frame(chr = character(N),
                   pos = integer(N),
                   ref = character(N),
                   alt = character(N),
                   maj_homo_prob = numeric(N),
                   min_homo_prob = numeric(N),
                   het_prob = numeric(N),
                   n_reads = integer(N),
                   reads = character(N))
  for (i in seq_len(N)) {
    snp <- snps[i, ]
    readsAtPos <- readsAtSnp(galign, snp)
    post_probs <- post_prob(snp, readsAtPos)
    
    df[i,]$chr <- snp$chr
    df[i,]$pos <- snp$pos
    df[i,]$alt <- snp$alt
    df[i,]$ref <- snp$ref
    df[i,]$maj_homo_prob <- round(post_probs[1], digits = 4)
    df[i,]$min_homo_prob <- round(post_probs[2], digits = 4)
    df[i,]$het_prob <- round(post_probs[3], digits = 4)
    df[i,]$n_reads <- length(readsAtPos$observations)
    df[i,]$reads <- paste(readsAtPos$observations, collapse=",")
  }
  df
}

snps_post_prob(putative_snps, gal)