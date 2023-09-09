#!/usr/bin/env Rscript

## Script to extract relevant information from vcf files (PL field)
## and prepare input data for entropy in the right format 
# March 2019 -- Vivaswat Shastry

## USAGE: Rscript inputdataformat.R file1.vcf


#### ENVIRONMENT SETUP ####
require(vcfR)
require(MASS)


#### FUNCTION SETUP #### 
## PCA function
#--------------------------------------------------------------------#
do_pca <- compiler::cmpfun(function(gmat){
  print("1. centering genotype matrix")
  gmn <- apply(gmat, 1, mean, na.rm=T)
  gmnmat <- matrix(gmn, nrow=nrow(gmat), ncol=ncol(gmat))
  gprime <- gmat-gmnmat ## remove mean
  
  print("2. finding pairwise covariance")
  gcovarmat <- matrix(NA, nrow=ncol(gmat), ncol=ncol(gmat))
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i==j){
        gcovarmat[i, j] <- cov(gprime[, i], gprime[, j],
                               use="pairwise.complete.obs")
      }
      else{
        gcovarmat[i, j] <- cov(gprime[, i], gprime[, j],
                               use="pairwise.complete.obs")
        gcovarmat[j, i] <- gcovarmat[i, j]
      }
    }
  }
  missing.inds <- unique(which(is.na(gcovarmat), arr.ind=T)[, 2])
  print("3. final principal comp analysis")
  return(list(prcomp(x=na.omit(gcovarmat),
                     center = TRUE, scale = FALSE), missing.inds))
})
#--------------------------------------------------------------------#


# point estimates to ancestry ldak estimates
#--------------------------------------------------------------------#
point_to_ldak <- compiler::cmpfun(function(genotype_likelihoods){
  
  print("Obtaining admixture proportions for chain initialization...")
  
  return.val <- do_pca(genotype_likelihoods)
  pcout <- return.val[1][[1]]
  miss.inds <- return.val[2][[1]]
  pcSummary <- summary(pcout)
  
  print("4. k-means, lda and writing initial admixture proportions to file")
  
  for(k in 2:3){
    init.admix <- matrix(0, nrow=nind, ncol=k)
    init.admix[miss.inds, ] <- rep(1/k, k)
    
    kn <- kmeans(pcout$x[, 1:5], k, iter.max=10, nstart=10,
                 algorithm="Hartigan-Wong")
    ldakn <- lda(x=pcout$x[, 1:5], grouping=kn$cluster, CV=TRUE)
    init.admix <- ldakn$posterior
    write.table(round(init.admix, 5),
                quote = FALSE, row.names = FALSE, col.names = FALSE,
                file = paste0("qk", k, "inds.txt"))
  }
  
})
#--------------------------------------------------------------------#


## mpgl to point estimate function
#--------------------------------------------------------------------#
mpgl_to_point <- compiler::cmpfun(function(genotype_likelihoods, mult_files){
  
  print("Converting genotype likelihoods into point estimates...")
  
  if(mult_files){
    total.mean.gl <- matrix(NA, nrow=nrow(genotype_likelihoods),
                            ncol=ncol(genotype_likelihoods))
    for(i in 1:nrow(genotype_likelihoods)){
      for(j in 1:ncol(genotype_likelihoods)){
        if(length(unique(unlist(genotype_likelihoods[i, j])))!=1){
          temp <- 10^(genotype_likelihoods[i, j][[1]]/-10)/(sum(10^(genotype_likelihoods[i, j][[1]]/-10)))
          if(length(temp)==2){
            total.mean.gl[i, j] <- sum(temp*(0:ploidy[j]))
          }
          else{
            total.mean.gl[i, j] <- sum(temp*(0:ploidy[j]))
          }
        }
      }
    }
    a.mean.gl <- total.mean.gl
  }else{
    a.mean.gl <- matrix(0, nrow=nloci, ncol=nind)
    a.mean.gl <- apply(genotype_likelihoods, c(1, 2),
                       function(x){temp <- sum(10^(-0.1*unlist(x)));sum((10^(-0.1*unlist(x)))*(0:max(ploidy)))/temp})
  }
  write.table(a.mean.gl, file="pntest_mean_gl.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep=" ")
  
  # convert to LDAK point estimates
  point_to_ldak(a.mean.gl)
  
})
#--------------------------------------------------------------------#


## multi-file vcf-to-mpgl conversion function
#--------------------------------------------------------------------#
multifile_convert <- compiler::cmpfun(function(args, a.pl, a.loci, a.ids){
  
  ploidy1 <- as.character(args[2])
  vcffile2 <- as.character(args[3])
  ploidy2 <- as.character(args[4])
  
  fname2 <- strsplit(vcffile2, ".vcf")
  a.vcf2 <- read.vcfR(vcffile2)
  
  a.ids2 <- colnames(a.vcf2@gt)[-1]
  write.table(a.ids2, file=paste0("inds_", fname2, ".txt"),
              row.names = FALSE, col.names = FALSE)
  
  a.pl2 <- extract.gt(a.vcf2, element='PL')
  a.pl2[is.na(a.pl2)] <- paste(rep("0", ploidy2+1), collapse=", ")
  
  a.loci2 <- paste0(a.vcf2@fix[, 'CHROM'], ';', a.vcf2@fix[, 'POS'])
  
  nind <- dim(a.pl)[2]+dim(a.pl2)[2]
  commonloci <- a.loci %in% a.loci2
  nloci <- dim(a.pl)[1]+dim(a.pl2)[1]-length(commonloci.idx)
  
  total.pl <- matrix("0", nrow=nloci, ncol=nind)
  total.pl[1:length(commonloci), 1:dim(a.pl)[2]] <- a.pl[commonloci, ]
  total.pl[(length(commonloci)+1):dim(a.pl)[1],
           1:dim(a.pl)[2]] <- a.pl[!commonloci, ]
  
  commonloci <- a.loci2 %in% a.loci
  total.pl[1:length(commonloci),
           (dim(a.pl)[2]+1):nind] <- a.pl2[commonloci, ]
  total.pl[(dim(a.pl)[1]+1):nloci,
           (dim(a.pl)[2]+1):nind] <- a.pl2[!commonloci, ]
  
  total.pl[(length(commonloci)+1):dim(a.pl)[1],
           (dim(a.pl)[2]+1):nind] <- paste(rep("0", ploidy+1), collapse=", ")
  total.pl[(dim(a.pl)[1]+1):nloci,
           1:dim(a.pl)[2]] <- paste(rep("0", ploidy2+1), collapse=", ")
  
  total.pl <- apply(total.pl, c(1, 2),
                    function(x){as.numeric(unlist(strsplit(x, split=", ")))})
  
  print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 
  
  loci.names <- a.loci[a.loci %in% a.loci2]
  loci.names <- append(a.loci[!(a.loci %in% a.loci2)], loci.names)
  loci.names <- append(a.loci2[!(a.loci2 %in% a.loci)], loci.names)
  
  unlink("combined.mpgl")
  cat(paste(ncol(total.pl), nrow(total.pl), "\n"), file="combined.mpgl")
  cat(paste(a.ids), file="combined.mpgl", append=T)
  for(locus in 1:nloci){
    cat(paste0("\n", loci.names[locus], " "), file="combined.mpgl",
        append = TRUE)
    write.table(unlist(total.pl[locus, ]), file="combined.mpgl", 
                append = TRUE, col.names = FALSE, row.names = FALSE, eol=" ")
  }
  
  # convert genotype likeligoods to points
  mpgl_to_point(total.gl, mult_files)
  
})
#--------------------------------------------------------------------#


## vcf to mpgl function
#--------------------------------------------------------------------#
vcf_to_mpgl <- compiler::cmpfun(function(vcffile, mult_files, args){
  
  print("Creating input genotype likelihood data for entropy...")
  
  fname <- strsplit(vcffile, ".vcf")[[1]]
  a.vcf <- read.vcfR(vcffile)
  
  a.ids <- colnames(a.vcf@gt)[-1]
  write.table(a.ids, file=paste0("inds_", fname, ".txt"),
              row.names = FALSE, col.names = FALSE)
  
  a.loci <- paste0(a.vcf@fix[, 'CHROM'], ';', a.vcf@fix[, 'POS'])
  
  a.pl <- extract.gt(a.vcf, element='PL')
  
  print(paste0(round(sum(is.na(a.pl))*100/(nrow(a.pl)*ncol(a.pl)), 1),
               "% missing data in file"))
  
  if (single_file){
    
    a.gl <- apply(a.pl, c(1, 2),
                  function(x){as.numeric(unlist(strsplit(x, split=", ")))})
    
    nind <- dim(a.pl)[2]
    nloci <- dim(a.pl)[1]
    ploidy <- rep(2, nind)
    
    print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 
    
    for(i in 1:nrow(a.pl)){
      for(j in 1:ncol(a.pl)){
        if(is.na(a.gl[1, i, j])){
          a.gl[, i, j] <- rep(0, ploidy[[1]]+1)
        }
      }
    }
    unlink(paste0(fname, ".mpgl"))
    cat(paste(nind, nloci, "\n"), file=paste0(fname, ".mpgl"))
    cat(paste(a.ids), file=paste0(fname, ".mpgl"), append=T)
    for(locus in 1:nloci){
      cat(paste0("\n", a.vcf@fix[locus, 'CHROM'], " "),
          file=paste0(fname, ".mpgl"), append = TRUE)
      write.table(unlist(t(a.gl[, locus, ])), file=paste0(fname, ".mpgl"),
                  append = TRUE, col.names = FALSE, row.names = FALSE, eol=" ")
    }
    
    # convert genotype likeligoods to points
    mpgl_to_point(a.gl, mult_files)
    
  } else {
    
    # run multifile function
    multifile_convert(args, a.pl, a.loci, a.ids)
    
  }
  
})
#--------------------------------------------------------------------#


## main function
#--------------------------------------------------------------------#
main <- compiler::cmpfun(function(){
  
  # retrieve positional command line arguments
  args <- commandArgs(TRUE)
  
  ## flag to check if we have mixed ploidy 
  mult_files <- F
  vcffile <- as.character(args[1])
  if (!file.exists("ploidy_inds.txt")) {
    cat("Assuming diploidy for each indv. in vcf \n")
    #quit()
  }
  if (length(args)>2) {
    single_file <- FALSE
  }
  
  # run the conversion, which involves a nesting of the above functions
  vcf_to_mpgl(vcffile, mult_files, args)
  
})
#--------------------------------------------------------------------#


#### EXECUTE THE PROGRAM #### 
#--------------------------------------------------------------------#
main()
#--------------------------------------------------------------------#
