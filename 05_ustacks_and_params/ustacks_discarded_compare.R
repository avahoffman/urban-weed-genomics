library(tidyverse)

compare <- read_csv("output/ustacks-discarded_samples_notes.csv")

compare <- 
  compare %>% 
  group_by(`ustacks discarded`) %>% 
  select(!c(sample, population)) %>% 
  summarize(across(everything(), mean))

write_csv(compare, "output/ustacks-discared_samples_summary.csv")
