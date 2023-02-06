library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)


make_filterplot <- function() {
  # Read in output from the radtag filtering part of sequence preparation
  df <- readr::read_csv("output/process_radtags-sample_output.csv") %>% 
    dplyr::select(Filename, low_quality, retained_reads, prop_sample_per_library) %>% 
    pivot_longer(
      cols = c(low_quality, retained_reads, prop_sample_per_library),
      names_to = "feature",
      values_to = "value"
    )
  
  # Reorder factor and add plotworthy names
  df$feature <-
    factor(
      df$feature,
      levels = c("low_quality", "retained_reads", "prop_sample_per_library"),
      labels = c(
        "Number of low quality reads (removed)",
        "Retained reads",
        "Proportion of sample makeup in sublibrary"
      )
    )
  
  # Create plot
  gg <- ggplot(data = df) +
    geom_density(aes(x = value), fill = "#2171b5", alpha = 0.6) +
    facet_wrap( ~ feature, ncol = 1, scales = "free") +
    theme_minimal() +
    xlab(NULL)
  gg
  
  # Write plot
  ggsave("figures/radtag_filter_stats.jpg")
  
  return(gg)
  
}


make_manual_discard_plot <- function(){
  # Read in output from the radtags filtering part of sequence preparation
  df <- 
    readr::read_csv("output/process_radtags-sample_output.csv") %>% 
    mutate(species = case_when (
      startsWith(Filename, "DS") ~ "DS",
      startsWith(Filename, "CD") ~ "CD",
      startsWith(Filename, "PA") ~ "PA",
      startsWith(Filename, "EC") ~ "EC",
      startsWith(Filename, "LS") ~ "LS",
      startsWith(Filename, "TO") ~ "TO",
      startsWith(Filename, "TE") ~ "TE")) %>% 
    filter(species != "TE") %>% 
    mutate(discarded = factor(case_when(
      prop_sample_per_library < 0.01 ~ "Samples discarded: Less than 1% of sublibrary",
      retained_reads < 1000000 ~ "Samples discarded: Less than 1 million reads",
      TRUE ~ "Samples kept"
    ))) %>% 
    dplyr::select(Filename, discarded, species) %>% 
    group_by(discarded, species) %>% count()
  
  # Create plot
  gg <- ggplot(data = df) +
    geom_col(aes(x = species, y = n), fill = "#2171b5", alpha = 0.6) +
    facet_wrap( ~ discarded, ncol = 1) +
    theme_minimal() +
    xlab(NULL)
  gg
  
  # Write plot
  ggsave("figures/radtag_samples_discarded.jpg")
  
  return(gg)
  
}
