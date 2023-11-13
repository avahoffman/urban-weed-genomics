require(tidyverse)

compile_rho_table <- function() {
  CD <-
    read.csv("output/population_stats/genodive_output_rho_CD.csv") %>% column_to_rownames("X")
  DS <-
    read.csv("output/population_stats/genodive_output_rho_DS.csv") %>% column_to_rownames("X")
  EC <-
    read.csv("output/population_stats/genodive_output_rho_EC.csv") %>% column_to_rownames("X")
  LS <-
    read.csv("output/population_stats/genodive_output_rho_LS.csv") %>% column_to_rownames("X")
  PA <-
    read.csv("output/population_stats/genodive_output_rho_PA.csv") %>% column_to_rownames("X")
  TO <-
    read.csv("output/population_stats/genodive_output_rho_TO.csv") %>% column_to_rownames("X")
  
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
  
  rho_stats[,c(1,3,2,4)]
  
  readr::write_csv(rho_stats, "output/population_stats/rho_all.csv")
}
