# This script filters out low coverage samples to be used in the next steps

library(readr) # read_csv write_csv
library(dplyr) # %>% filter select arrange mutate case_when

retained_reads_threshold <- 1000000 # need more than 1000000 reads per sample
prop_sample_threshold <- 0.01 # need more than 1% of the sublibrary (e.g., 1% of 1-1, 1-2 etc.)

read_csv("output/process_radtags-sample_output.csv") %>%
  filter(
    retained_reads > retained_reads_threshold &
      prop_sample_per_library > prop_sample_threshold
  ) %>%
  select(Filename) %>% arrange(Filename) %>%
  mutate(
    Population = case_when(
      grepl("\\.BA\\.", Filename) ~ "Baltimore",
      grepl("\\.BO\\.", Filename) ~ "Boston",
      grepl("\\.LA\\.", Filename) ~ "Los Angeles",
      grepl("\\.MN\\.", Filename) ~ "Minneapolis",
      grepl("\\.PX\\.", Filename) ~ "Phoenix",
      TRUE ~ "None"
    )
  ) %>%
  write_csv("output/process_radtags-kept_samples.csv")

read_csv("output/process_radtags-sample_output.csv") %>%
  filter(
    retained_reads <= retained_reads_threshold |
      prop_sample_per_library <= prop_sample_threshold
  ) %>%
  select(!(`...1`)) %>%
  write_csv("output/process_radtags-discarded_samples.csv")
