#!/usr/bin/env Rscript

## USAGE: Rscript inputdataformat.R file1.vcf


#### ENVIRONMENT SETUP ####
library(vcfR)
library(MASS)
library(dplyr)
library(tidyr)
library(stringr)


## helper function to convert each Phred-scaled likelihood 
#--------------------------------------------------------------------#
to_numvec <- compiler::cmpfun(function(item) {
  return(
    item %>%
      strsplit(split = ",") %>%
      unlist() %>%
      trimws %>%
      as.numeric()
  )
})
#--------------------------------------------------------------------#


## vcf to mpgl function
#--------------------------------------------------------------------#
vcf_to_mpgl <- compiler::cmpfun(function(vcf_file){
  
  print("Creating input genotype likelihood data for entropy.")
  
  # parse out the simple name for the input VCF
  simple_name <- vcf_file %>%
    basename() %>%
    str_remove_all(., ".vcf")
  
  # bring the variant data into memory as a vcfR object
  vcf_data <- read.vcfR(vcf_file)
  
  # parse out the sample IDs and the varying loci
  sample_ids <- colnames(vcf_data@gt)[-1]
  sample_size <- length(sample_ids)
  loci <- paste(vcf_data@fix[, 'CHROM'],
                vcf_data@fix[, 'POS'],
                sep="_")
  n_loci <- length(loci)
  
  # Extract phred-scaled likelihoods, equivalent to:
  # PL = -10 * log{P(Genotype | Data)}
  sample_pl <- extract.gt(vcf_data, element='PL')
  
  # Extract sample size and variant count to get the dimensionality of the data
  sample_size <- ncol(sample_pl)
  variant_count <- nrow(sample_pl)
  print(
    paste(sample_size, "samples among", variant_count, "variants.", sep = " ")
  )
  dimensionality <- sample_size * variant_count
  
  # Calculate the percentage of possible data points that are missing
  percent_missing <- sum(is.na(sample_pl)) * 100
  print(
    paste(
      round(dimensionality / dimensionality, 1),
      "% missing data in file",
      sep = ""
    )
  )
  
  # convert phred-likelihoods into numeric vectors for every item in the matrix
  sample_gl <- apply(sample_pl, c(1, 2), function(x) { to_numvec(x) })
  
  # run some checks to make sure conversions worked
  stopifnot(sample_size == dim(sample_pl)[2])
  stopifnot(n_loci == dim(sample_pl)[1])
  
  # create a vector with a ploidy value for each sample. Here we assume samples
  # are all diploid.
  ploidy <- rep(2, sample_size)
  
  # 
  for (i in 1:nrow(sample_pl)) {
    for (j in 1:ncol(sample_pl)) {
      if (is.na(sample_gl[1, i, j])) {
        sample_gl[, i, j] <- rep(0, ploidy[[1]]+1)
      }
    }
  }
  
  # prepare output file
  output_name <- paste(simple_name, "mpgl", sep = ".")
  if (file.exists(output_name)) {
    unlink(paste(simple_name, ".mpgl"))
  }
  
  # prepare the header for the MPGL output file
  cat(paste(sample_size, n_loci, "\n"), file = output_name)
  cat(paste(sample_ids), file = output_name, append = T)
  
  # append the rows of likelihoods for each locus to complete the file
  for (i in 1:n_loci) {
    cat(paste0("\n", vcf_data@fix[i, 'CHROM'], " "),
        file = output_name, append = TRUE)
    write.table(unlist(t(sample_gl[, i, ])), file = output_name,
                append = TRUE, col.names = FALSE, row.names = FALSE, eol=" ")
  }
  
})
#--------------------------------------------------------------------#


## main function
#--------------------------------------------------------------------#
main <- compiler::cmpfun(function(){
  
  # retrieve positional command line arguments
  args <- commandArgs(TRUE)
  vcf_file <- args[1]
  stopifnot(file.exists(vcf_file))
  
  # run the conversion, which involves a nesting of the above functions
  vcf_to_mpgl(vcf_file)
  
})
#--------------------------------------------------------------------#


#### EXECUTE THE PROGRAM #### 
#--------------------------------------------------------------------#
main()
#--------------------------------------------------------------------#

