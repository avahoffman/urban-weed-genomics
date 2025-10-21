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
  rho_p <- clean_rho_table(delim = "_pval") %>% rename(`p-value` = rho)
  
  rho_stats <- full_join(rho_vals, rho_p)
  
  rho_stats$`adjusted p-value` <- round(c(
    p.adjust(rho_stats$`p-value`[1:3], method = "BH"), # CD
    p.adjust(rho_stats$`p-value`[4:9], method = "BH"), # DS
    p.adjust(rho_stats$`p-value`[10:12], method = "BH"), # EC
    p.adjust(rho_stats$`p-value`[13:22], method = "BH"), # LS
    p.adjust(rho_stats$`p-value`[23:28], method = "BH"), # PA
    p.adjust(rho_stats$`p-value`[29:38], method = "BH") # TO
  ), 4)
  
  readr::write_csv(rho_stats, "output/population_stats/rho_all.csv")
}


read_rhos_dat <- function() {
  CD <-
    read.csv("output/population_stats/genodive_output_rho_CD.csv")
  DS <-
    read.csv("output/population_stats/genodive_output_rho_DS.csv")
  EC <-
    read.csv("output/population_stats/genodive_output_rho_EC.csv")
  LS <-
    read.csv("output/population_stats/genodive_output_rho_LS.csv")
  PA <-
    read.csv("output/population_stats/genodive_output_rho_PA.csv")
  TO <-
    read.csv("output/population_stats/genodive_output_rho_TO.csv")
  
  rhos_ <- list(
    "CD" = CD,
    "DS" = DS,
    "EC" = EC,
    "LS" = LS,
    "PA" = PA,
    "TO" = TO
  )
  
  return(rhos_)
}


# Function currently not being used:
plot_rhos <- function() {
  rhos <- read_rhos_dat()
  
  rho_i <- rhos$EC
  rownames(rho_i) <- rho_i$X 
  rho_i <- rho_i[,-1]
  rho_i[upper.tri(rho_i)] <- NA
  rho_i[rho_i == 0] <- NA
  rho_i$v1 <- rownames(rho_i)
  rho_i <- melt(rho_i)
  
  ggplot(rho_i, aes(v1, variable)) +
    theme_minimal() +
    geom_tile(aes(fill = value), color='white') +
    scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', limit = c(0,0.5)) +
    theme(axis.text.x=element_text(angle=90),
          axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'))
}