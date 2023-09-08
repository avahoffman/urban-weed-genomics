# library(adegenet)
# library(hierfstat)
# library(poppr)

read_dat <- function() {
  CD <- 
    read.csv("output/population_stats/genodive_output_rho_CD.csv", row.names = 1)
  DS <-
    read.csv("output/population_stats/genodive_output_rho_DS.csv", row.names = 1)
  EC <-
    read.csv("output/population_stats/genodive_output_rho_EC.csv", row.names = 1)
  LS <-
    read.csv("output/population_stats/genodive_output_rho_LS.csv", row.names = 1)
  PA <-
    read.csv("output/population_stats/genodive_output_rho_PA.csv", row.names = 1)
  TO <-
    read.csv("output/population_stats/genodive_output_rho_TO.csv", row.names = 1)
  
  rhos_ <- list(
    "CD" = as.matrix(CD),
    "DS" = as.matrix(DS),
    "EC" = as.matrix(EC),
    "LS" = as.matrix(LS),
    "PA" = as.matrix(PA),
    "TO" = as.matrix(TO)
  )
  
  return(rhos_)
}

