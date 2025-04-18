#' This script reads in outputs from the polyRAD prep and performs population
#' levelstatistics for population genetics.

library(polysat) # calcPopDiff
library(polyRAD) # CRAN v2.0.0 # GetWeightedMeanGenotypes nLoci GetLoci
library(readr) # read_rds
library(magrittr) # %>%


do_continental_stats <- function(spp_) {
  # Read in polyrad object and genind object
  gen_ <-
    readRDS(paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno_genind.rds"))
  wd_ <- paste0("SNP_data/", spp_, "/", spp_)
  polyrad_dat <- read_rds(paste0(wd_, "_estimatedgeno_RADdata.rds"))
  
  # Extract population info
  popmap <-
    stringr::str_extract(names(polyrad_dat$taxaPloidy),
                         "(?<=...)[:graph:]{2,4}(?=\\..)")
  
  # aggregate by population
  matrix_all <-
    GetWeightedMeanGenotypes(polyrad_dat, omit1allelePerLocus = FALSE)
  freq_by_pop <- do.call(rbind, by(matrix_all, popmap, colMeans))
  
  # Make sure everything sums to 1 for frequencies -- see https://lvclark.r-universe.dev/articles/polyRADtutorials/population_genetics.html
  for (i in seq_len(nLoci(polyrad_dat))) {
    thesealleles <- which(polyrad_dat$alleles2loc == i)
    # normalize to sum to one
    temp <- freq_by_pop[, thesealleles]
    temp <- sweep(temp, 1, rowSums(temp), "/")
    freq_by_pop[, thesealleles] <- temp
    # change allele names
    thisloc <- GetLoci(polyrad_dat)[i]
    thisloc <-
      gsub("\\.", "_", thisloc) # remove periods from locus names
    colnames(freq_by_pop)[thesealleles] <-
      paste(thisloc, seq_along(thesealleles), sep = ".")
  }
  
  genomes <-
    data.frame(pop = popmap) %>% dplyr::count(pop) %>% dplyr::mutate(n = n * gen_$ploidy[1])
  freq_df <- data.frame(Genomes = genomes$n, freq_by_pop)
  
  # Make sure things look ok
  freq_df[, 1:5]
  
  # Drop BO for EC -- only one sample from Boston
  freq_df <- freq_df[-2,]
  
  # Calc stats (continental)
  statistic_ <- c("JostD", "Gst", "Fst")
  my_JostD <-
    calcPopDiff(freq_df, metric = "Jost's D", global = TRUE)
  my_JostD
  my_Gst <- calcPopDiff(freq_df, metric = "Gst", global = TRUE)
  my_Gst
  my_Fst <- calcPopDiff(freq_df, metric = "Fst", global = TRUE)
  my_Fst
  values_ <- c(my_JostD, my_Gst, my_Fst)
  
  out_ <- data.frame(statistic = statistic_,
                     value = values_)
  
  write.csv(out_,
            paste0("output/population_stats/", spp_, "_continental_stats.csv"))
}


do_all_continental_stats <- function(){
  do_continental_stats(spp_ = "CD")
  do_continental_stats(spp_ = "DS")
  do_continental_stats(spp_ = "EC")
  do_continental_stats(spp_ = "LS")
  do_continental_stats(spp_ = "PA")
  do_continental_stats(spp_ = "TO")  
}




