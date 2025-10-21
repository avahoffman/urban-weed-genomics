# library(adegenet)
# library(hierfstat)
library(poppr)

read_dat <- function() {
  CD <-
    readRDS("SNP_data/CD/CD_estimatedgeno_genind.rds")
  DS <-
    readRDS("SNP_data/DS/DS_estimatedgeno_genind.rds")
  EC <-
    readRDS("SNP_data/EC/EC_estimatedgeno_genind.rds")
  LS <-
    readRDS("SNP_data/LS/LS_estimatedgeno_genind.rds")
  PA <-
    readRDS("SNP_data/PA/PA_estimatedgeno_genind.rds")
  TO <-
    readRDS("SNP_data/TO/TO_estimatedgeno_genind.rds")
  
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


calc_rbarD <- function() {
  # Nice explanation here!
  # https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html
  genind_ <- read_dat()
  
  # Create an empty dataframe so values can be appended
  rbarD_dat <- tibble::tibble()
  # Loop through species
  for (i in 1:6) {
    genind_dat <- genind_[[i]]
    
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
  
  readr::write_csv(rbarD_dat, "output/population_stats/rbarD.csv")
}
