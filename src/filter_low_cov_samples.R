# This script filters out low coverage samples to be used in the next steps

library(readr)
library(dplyr)

retained_reads_threshold <- 1000000 # need more than 1000000 reads per sample
prop_sample_threshold <- 0.01 # need more than 1% of the sublibrary (e.g., 1% of 1-1, 1-2 etc.)

read_csv("output/process_radtags-sample_output.csv") %>% 
  filter(retained_reads > retained_reads_threshold &  
         prop_sample_per_library > prop_sample_threshold) %>% 
  select(Filename) %>% arrange(Filename) %>% 
  write_csv("output/process_radtags-kept_samples.csv")

read_csv("output/process_radtags-sample_output.csv") %>% 
  filter(retained_reads <= retained_reads_threshold |  
         prop_sample_per_library <= prop_sample_threshold) %>% 
  write_csv("output/process_radtags-discarded_samples.csv")