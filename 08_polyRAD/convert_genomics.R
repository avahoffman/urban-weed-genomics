# Convert polyRAD output to common genetics formats
library(polyRAD)
library(readr)
library(stringr)

convert_spp <- function(spp_){
  # Read in polyRAD output (RADdata object)
  wd_ <- paste0("SNP_data/", spp_, "/", spp_)
  polyrad_dat <- read_rds(paste0(wd_,"_estimatedgeno_RADdata.rds"))
  
  # Save genind format file if file isn't there already
  if (!file.exists(paste0(wd_, "_estimatedgeno_genind.rds"))) {
    genind_format_dat <-
      Export_adegenet_genind(polyrad_dat)
    # Extract population information (two character abbreviation)
    popmap <- stringr::str_extract(names(polyrad_dat$taxaPloidy), "(?<=...)[:graph:]{2,4}(?=\\..)")
    pop(genind_format_dat) <- popmap
    write_rds(genind_format_dat, paste0(wd_, "_estimatedgeno_genind.rds"))
  }
  
  # Save structure format file
  if (!file.exists(paste0(wd_, "_estimatedgeno.structure"))) {
    # Create structure file and re-import it
    Export_Structure(polyrad_dat, file = paste0(wd_, "_estimatedgeno.structure"))
    struct_ <- 
      read.delim(paste0(wd_, "_estimatedgeno.structure"), sep = "\t", check.names = F)
    colnames(struct_)[1] <- "V1"
    
    # Import population information
    popmap <- stringr::str_extract(names(polyrad_dat$taxaPloidy), "(?<=...)[:graph:]{2,4}(?=\\..)")
    
    # Recode (Structure wants integers)
    if(any(popmap == "BA")) popmap[popmap == "BA"] <- 1
    if(any(popmap == "BO")) popmap[popmap == "BO"] <- 2
    if(any(popmap == "LA")) popmap[popmap == "LA"] <- 3
    if(any(popmap == "MN")) popmap[popmap == "MN"] <- 4
    if(any(popmap == "PX")) popmap[popmap == "PX"] <- 5
    
    struct_w_pop <- merge(popmap, struct_, all = TRUE)
    
    # Structure doesn't want column names for non-locus info
    colnames(struct_w_pop)[1:2] <- c("", "")
    
    # Write structure-compatible file
    write.table(struct_w_pop, paste0(wd_, "_estimatedgeno.structure"), sep = "\t", quote = F, row.names = F)
  }
  
}


convert_all <- function(){
  convert_spp(spp_ = "CD")
  convert_spp(spp_ = "DS")
  convert_spp(spp_ = "EC")
  convert_spp(spp_ = "LS")
  convert_spp(spp_ = "PA")
  convert_spp(spp_ = "TO")
}
