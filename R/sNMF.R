# sNMF alternative to Structure for validation

# BiocManager::install("LEA")
library(LEA)

convert_from_structure_for_sNMF <- function(){
  base_dir <- paste0(here::here(), "/SNP_data/")
  do_conversion <- function(input.file) {
    struct2geno(input.file = input.file, ploidy = 2, FORMAT = 2, extra.row = 0, extra.col = 2)
  }
  
  # EC - polyRAD
  input.file <- paste0(base_dir, "EC/EC_estimatedgeno_noheader.structure")
  do_conversion(input.file)
  # EC - SNPs from Stacks
  input.file <- paste0(base_dir, "EC/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # LS - polyRAD
  input.file <- paste0(base_dir, "LS/LS_estimatedgeno_noheader.structure")
  do_conversion(input.file)
  # LS - SNPs from Stacks
  input.file <- paste0(base_dir, "LS/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # DS - polyRAD
  # Note that I needed to use GenoDive to convert from n-ploidy lines per individual 
  # in the structure file to 2 lines per individual..
  input.file <- paste0(base_dir, "DS/DS_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # DS - SNPs from Stacks
  # Output from Stacks is already 2 lines per individual..
  input.file <- paste0(base_dir, "DS/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
}


do_sNMF <- function(){
  # Do sNMF
  spp_ <- "DS"
  K_ <- 4
  #cols_ <- c("orange","violet", "lightgreen")
  cols_ <- c("orange","violet","lightgreen","lightblue")
  cols_ <- cols_[1:K_]
  
  # infile <- paste0(
  #     here::here(),
  #     "/SNP_data/",
  #     spp_,
  #     "/populations_20_pct/populations_noheader.structure.geno")
  # infile <- paste0(
  #     here::here(),
  #     "/SNP_data/",
  #     spp_,
  #     "/",
  #     spp_,
  #     "_estimatedgeno_noheader.structure.geno")
  infile <- paste0(
    here::here(),
    "/SNP_data/",
    spp_,
    "/",
    spp_,
    "_estimatedgeno_noheader_genodivefix.structure.geno")
  obj.snmf <- snmf(
    infile,
    K = K_,
    alpha = 100,
    project = "new",
    entropy = TRUE
  )
  qmatrix <- Q(obj.snmf, K = K_)
  barplot(t(qmatrix), col = cols_, border = NA, space = 0, ylab = "Admixture coefficients", xlab = "Individuals")
}