# How many samples by city and species?
library(tidyverse)

get_sample_size <- function(){
  # Data output from PolyRAD. This is the last stage at which samples were dropped.
  
  dat <- bind_rows(
    read_tsv(paste0("SNP_data/CD/CD_estimatedgeno.structure"))[,1],
    read_tsv(paste0("SNP_data/DS/DS_estimatedgeno.structure"))[,1],
    read_tsv(paste0("SNP_data/EC/EC_estimatedgeno.structure"))[,1],
    read_tsv(paste0("SNP_data/LS/LS_estimatedgeno.structure"))[,1],
    read_tsv(paste0("SNP_data/PA/PA_estimatedgeno.structure"))[,1],
    read_tsv(paste0("SNP_data/TO/TO_estimatedgeno.structure"))[,1]
    )
  
  dat <- dat %>% distinct()
  
  dat <-
    dat %>% 
    mutate(spp = stringr::str_sub(`...1`, 1, 2),
           city = stringr::str_sub(`...1`, 4, 5))
  
  dat %>% count(spp)
  
  dat %>% count(spp, city) %>% print(n = 30)

  }
