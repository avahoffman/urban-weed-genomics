library(ggplot2)
library(tidyr)
library(magrittr)


make_cloneplot <- function() {
  # Read in output from the clone filtering part of sequence preperation
  df <- readr::read_csv("output/clone_filter_out.csv") %>%
    pivot_longer(
      cols = c(input, output, PCR_clones_removed),
      names_to = "feature",
      values_to = "value"
    )
  
  # Reorder factor and add plotworthy names
  df$feature <-
    factor(
      df$feature,
      levels = c("input", "PCR_clones_removed", "output"),
      labels = c(
        "Reads Input (per sublibrary)",
        "PCR clones removed",
        "Reads Output (per sublibrary)"
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
