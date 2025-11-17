library(tidyverse)
library(polyRAD)

get_var_explained <- function(object){
  # Copied and modified slightly from polyRAD addPCA function here: https://github.com/lvclark/polyRAD/blob/main/R/classes_methods.R
  nPcsInit = 10
  maxR2changeratio = 0.05
  minPcsOut = 1
  
  if(minPcsOut > nPcsInit){
    stop("minPcsOut can not be greater than nPcsInit.")
  }
  # matrix for input to PCA; depth ratios or posterior probs
  if(!CanDoGetWeightedMeanGeno(object)){
    genmat <- object$depthRatio[,-OneAllelePerMarker(object)]
  } else {
    genmat <- GetWeightedMeanGenotypes(object, omit1allelePerLocus = TRUE,
                                       naIfZeroReads = FALSE)
  }
  
  # replace NaN with NA
  genmat[is.na(genmat)] <- NA
  # check that no columns are completely NA
  if(any(colMeans(is.na(genmat)) == 1)){
    message("Alleles with completely missing data:")
    cat(colnames(genmat)[colMeans(is.na(genmat)) == 1], sep = "\n")
    stop("Alleles found with completely missing data.")
  }
  
  genfreq <- colMeans(genmat, na.rm = TRUE)
  # if any individuals are completely missing, fill them in with mean.
  # This can happen with blanks or very low quality samples when looping
  # through the genome in chunks.
  missind <- which(rowMeans(is.na(genmat)) == 1)
  for(i in missind){
    genmat[i,] <- genfreq
  }
  
  # remove non-variable sites
  genmat <- genmat[, which(genfreq > 0 & genfreq < 1)]
  
  # adjust number of PC axes if necessary
  if(nPcsInit > dim(genmat)[2]){
    nPcsInit <- dim(genmat)[2]
  }
  
  # run principal components analysis
  pc <- pcaMethods::pca(genmat, method = "ppca", nPcs = nPcsInit)
  # get rate of change in R2 values
  roc <- pc@R2[1:(nPcsInit - 1)] - pc@R2[2:nPcsInit]
  cutoff <- which(roc < roc[1] * maxR2changeratio)
  if(length(cutoff) == 0){
    cutoff <- nPcsInit
  }
  # make sure number of PCs meets minimum threshold
  if(cutoff[1] < minPcsOut){
    cutoff <- minPcsOut
  }
  
  return(pc@R2[1:2])
}

CD_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/CD_estimatedgeno_RADdata.rds") %>% get_var_explained()
DS_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/DS_estimatedgeno_RADdata.rds") %>% get_var_explained()
EC_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/EC_estimatedgeno_RADdata.rds") %>% get_var_explained()
LS_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/LS_estimatedgeno_RADdata.rds") %>% get_var_explained()
PA_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/PA_estimatedgeno_RADdata.rds") %>% get_var_explained()
TO_varexplained <- readRDS("../urban-weed-genomics-local/polyRAD_output/TO_estimatedgeno_RADdata.rds") %>% get_var_explained()

