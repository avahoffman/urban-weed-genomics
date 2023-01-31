library(adegenet)
library(hierfstat)
library(poppr)

read_dat <- function() {
  CD <-
    read.genepop("SNP_data/CD/populations_20_pct/populations.snps.gen")
  DS <-
    read.genepop("SNP_data/DS/populations_20_pct/populations.snps.gen")
  EC <-
    read.genepop("SNP_data/EC/populations_20_pct/populations.snps.gen")
  LS <-
    read.genepop("SNP_data/LS/populations_20_pct/populations.snps.gen")
  PA <-
    read.genepop("SNP_data/PA/populations_20_pct/populations.snps.gen")
  TO <-
    read.genepop("SNP_data/TO/populations_20_pct/populations.snps.gen")
  
  genind_ <- list(
    "CD" = CD,
    "DS" = DS,
    "EC" = EC,
    "LS" = LS,
    "PA" = PA,
    "TO" = TO
  )
  
  return(genind_)
}


calc_fst <- function() {
  read_dat()
  
  fst_dat <- tibble::tibble()
  # Loop through species
  for (i in 1:6) {
    fst_dat <- fst_dat %>% rbind(tibble::tibble(
    "spp" = names(genind_)[i],
    "Fst" = round(wc(genind_[[i]])$FST, 3)
    ))
  }
  
  readr::write_csv(fst_dat, "output/population_stats-fst.csv")
}



calc_rbarD <- function() {
  # Nice explanation here!
  # https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html
  genind_ <- read_dat()
  
  # Create an empty dataframe so values can be appended
  rbarD_dat <- tibble::tibble()
  # Loop through species
  for (i in 1:6) {
    genind_dat <- genind_[[i]]
    
    # Clean up population names 
    pops_clean <-
      genind_dat$pop %>%
      stringr::str_extract("(?<=..)[:alpha:]{2}")
    genind_dat$pop <-
      factor(pops_clean,
             levels = unique(pops_clean))
    
    # Loop thru cities in each species
    # Takes a long time..
    for (pop_ in levels(genind_dat$pop)) {
      genind_dat_city <- popsub(genind_dat, pop_)
      set.seed(1234)
      spp_ia <- poppr::ia(genind_dat_city, sample = 999)
      spp_ia_ <- cbind(
        spp = names(genind_)[i],
        city = pop_,
        n = length(genind_dat_city$pop),
        tibble::tibble(!!!spp_ia)
      )
      # Append data
      rbarD_dat <- rbarD_dat %>% rbind(spp_ia_)
    }
  }
  
  readr::write_csv(rbarD_dat, "output/population_stats-rbarD.csv")
}
