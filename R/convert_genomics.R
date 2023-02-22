# Convert polyRAD output to common genetics formats
library(polyRAD)
library(readr)

convert_spp <- function(spp_){
  spp_ <- "CD"
  wd_ <- paste0("SNP_data/", spp_, "/", spp_)
  polyrad_dat <- read_rds(paste0(wd_,"_estimatedgeno_RADdata.rds"))
  
  # Save genind format file
  if (!file.exists(paste0(wd_, "_estimatedgeno_genind.rds"))) {
    genind_format_dat <-
      Export_adegenet_genind(polyrad_dat)
    write_rds(genind_format_dat, paste0(wd_, "_estimatedgeno_genind.rds"))
  }
  
  # Save structure format file
  if (!file.exists(paste0(wd_, "_estimatedgeno.structure"))) {
    Export_Structure(polyrad_dat, file = paste0(wd_, "_estimatedgeno.structure"))
  }
  
  vcf_dat <- RADdata2VCF(polyrad_dat)
  genotype_dat <- GetProbableGenotypes(polyrad_dat)
  snowcrab <- genomic_converter(
    data =  genind_format_dat, 
    strata = read_delim(paste0("SNP_data/", spp_, "/popmap_", spp_, ".txt"), col_names = c("INDIVIDUALS", "STRATA"), delim = "\t"),
    output = c("vcf")
  )
}


convert_all <- function(){
  convert_spp("CD")
  convert_spp("DS")
  convert_spp("EC")
  convert_spp("LS")
  convert_spp("PA")
  convert_spp("TO")
}
