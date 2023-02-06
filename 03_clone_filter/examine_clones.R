library(ggplot2)
library(tidyr)
library(magrittr)


make_cloneplot <- function() {
  # Read in output from the clone filtering part of sequence preparation
  df <- readr::read_csv("output/clone_filter_out.csv") %>%
    dplyr::mutate(pct_clones = PCR_clones_removed/input * 100) %>% 
    pivot_longer(
      cols = c(input, output, pct_clones),
      names_to = "feature",
      values_to = "value"
    )
  
  # Reorder factor and add plotworthy names
  df$feature <-
    factor(
      df$feature,
      levels = c("input", "pct_clones", "output"),
      labels = c(
        "Raw reads input (per sublibrary)",
        "Percent PCR clones (removed)",
        "Reads retained (per sublibrary)"
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
  ggsave("figures/pcr_clones.jpg")
  
  return(gg)
  
}
