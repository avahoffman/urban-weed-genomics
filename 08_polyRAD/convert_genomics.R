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
    popmap <- read.table(paste0("SNP_data/", spp_, "/popmap_", spp_, "_polyrad.txt"), sep = "\t")
    
    # Recode (Structure wants integers)
    if(any(popmap$V2 == "Baltimore")) popmap[popmap$V2 == "Baltimore",]$V2 <- 1
    if(any(popmap$V2 == "Boston")) popmap[popmap$V2 == "Boston",]$V2 <- 2
    if(any(popmap$V2 == "Los Angeles")) popmap[popmap$V2 == "Los Angeles",]$V2 <- 3
    if(any(popmap$V2 == "Minneapolis")) popmap[popmap$V2 == "Minneapolis",]$V2 <- 4
    if(any(popmap$V2 == "Phoenix")) popmap[popmap$V2 == "Phoenix",]$V2 <- 5
    
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
