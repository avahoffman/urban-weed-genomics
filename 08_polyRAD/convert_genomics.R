#' Convert polyRAD output to common genetics formats
#' Using the check_coverage function visually plots HindHe stat versus coverage.
#' This allowed us to remove any samples with very low coverage that might
#' interfere with downstream applications like Structure.

library(polyRAD)
library(readr)
library(stringr)


check_coverage <- function(spp_) {
  wd_ <- paste0("SNP_data/", spp_, "/", spp_)
  polyrad_dat <- read_rds(paste0(wd_, "_estimatedgeno_RADdata.rds"))
  
  hh <- HindHe(polyrad_dat)
  
  HindHe <- rowMeans(hh, na.rm = TRUE)
  Depth <- rowSums(polyrad_dat$locDepth)
  
  dd <- data.frame(HindHe, Depth)
  
  ggplot(dd, aes(x = Depth, y = HindHe)) +
    geom_point() +
    ggtitle("Read depth and Hind/He across individuals") +
    scale_x_log10()
}


convert_spp <- function(spp_){
  # Read in polyRAD output (RADdata object)
  wd_ <- paste0("SNP_data/", spp_, "/", spp_)
  polyrad_dat <- read_rds(paste0(wd_,"_estimatedgeno_RADdata.rds"))
  
  # ! Remove the following low coverage or outlying samples.
  lowcov <- c(
    "CD.BA.PSP.M.1",
    "CD.BA.DHI.U.2",
    "CD.BA.DHI.U.3",
    "CD.BA.RG-1.M.5",
    "CD.BA.RG-1.M.4",
    "DS.BO.WL1.M.4",
    "DS.BO.I1.U.3",
    "EC.BO.R4.U.1",
    "LS.BO.HC2.M.5",
    "LS.BO.LC4.M.3",
    "LS.BO.R2.U.4",
    "LS.BO.R2.U.1",
    "PA.BA.LH-3.M.4",
    "PA.BA.AA.U.3",
    "PA.BA.AA.U.4",
    "PA.PX.RPP.U.2",
    "PA.BO.HC2.M.4",
    "PA.PX.RPP.U.1",
    "TO.BA.TRC.U.1",
    "TO.BA.TRC.U.3",
    "TO.BO.R4.U.1",
    "TO.BA.TRC.U.2",
    "TO.BO.R4.U.2",
    "TO.BO.R2.U.2"
  )

  # Drop these samples
  keepsamples <- names(polyrad_dat$taxaPloidy)[!(names(polyrad_dat$taxaPloidy) %in% lowcov) ]
  polyrad_dat <- SubsetByTaxon(polyrad_dat, keepsamples)
  
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
    
    # Extract population information
    struct_$pop <- stringr::str_extract(struct_$V1, "(?<=...)[:graph:]{2,4}(?=\\..)")
    struct_ <- struct_ %>% dplyr::relocate(pop, .after = V1)
    
    # Recode (Structure wants integers)
    struct_w_pop <- struct_ %>% 
      dplyr::mutate(pop = dplyr::case_when(
        pop == "BA" ~ 1,
        pop == "BO" ~ 2,
        pop == "LA" ~ 3,
        pop == "MN" ~ 4,
        pop == "PX" ~ 5
      ))
    
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
