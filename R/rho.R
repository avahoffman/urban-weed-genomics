require(tidyverse)

clean_rho_table <- function(delim = "") {
  CD <-
    read.csv(paste0("output/population_stats/genodive_output_rho_CD",delim,".csv")) %>% column_to_rownames("X")
  DS <-
    read.csv(paste0("output/population_stats/genodive_output_rho_DS",delim,".csv")) %>% column_to_rownames("X")
  EC <-
    read.csv(paste0("output/population_stats/genodive_output_rho_EC",delim,".csv")) %>% column_to_rownames("X")
  LS <-
    read.csv(paste0("output/population_stats/genodive_output_rho_LS",delim,".csv")) %>% column_to_rownames("X")
  PA <-
    read.csv(paste0("output/population_stats/genodive_output_rho_PA",delim,".csv")) %>% column_to_rownames("X")
  TO <-
    read.csv(paste0("output/population_stats/genodive_output_rho_TO",delim,".csv")) %>% column_to_rownames("X")
  
  CD[upper.tri(CD, diag = T)] <-
    DS[upper.tri(DS, diag = T)] <-
    EC[upper.tri(EC, diag = T)] <-
    LS[upper.tri(LS, diag = T)] <-
    PA[upper.tri(PA, diag = T)] <- TO[upper.tri(TO, diag = T)] <- NA
  
  rho_stats <- rbind(
    reshape2::melt(CD %>% rownames_to_column() %>% mutate(spp = "CD")),
    reshape2::melt(DS %>% rownames_to_column() %>% mutate(spp = "DS")),
    reshape2::melt(EC %>% rownames_to_column() %>% mutate(spp = "EC")),
    reshape2::melt(LS %>% rownames_to_column() %>% mutate(spp = "LS")),
    reshape2::melt(PA %>% rownames_to_column() %>% mutate(spp = "PA")) %>% filter(rowname != "MN", variable != "MN"),
    reshape2::melt(TO %>% rownames_to_column() %>% mutate(spp = "TO"))
  ) %>% na.omit(value)
  
  rho_stats <- rho_stats %>% 
    relocate(spp, .before = rowname) %>% 
    rename(Species = spp, City1 = rowname, City2 = variable, rho = value) %>% 
    arrange(Species, -rho)
  
  return(rho_stats)
}


compile_rho_table <- function(){
  # Read in the stats only
  rho_vals <- clean_rho_table(delim = "")
  # Get p-values
  rho_p <- clean_rho_table(delim = "_pval")
  
  readr::write_csv(rho_stats, "output/population_stats/rho_all.csv")
}
