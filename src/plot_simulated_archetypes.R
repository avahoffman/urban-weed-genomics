# SIMULATED AA plots (not real data - useful for presenting this work to lay audience)

library(tidyverse)
library(ggplot2)
library(ggh4x)

make_simulated_archetypes_plot <- function(out_file,
                                           primary_dir = here::here(),
                                           admixed = TRUE,
                                           width = 12,
                                           height = 3.5) {
  primary_dir <- here::here()
  
  k_ = 3
  
  df <-
    data.frame(
      sample = factor(paste0("Sample ", seq(1:300)), levels = paste0("Sample ", seq(1:300))),
      breed = c(rep("Breed 1", 100), rep("Breed 2", 100), rep("Breed 3", 100)),
      k = k_
    )
  
  if (admixed == T) {
    df$X1 <- runif(n = 300, min = 0, max = 1)
    df$X2 <- 1 - df$X1
    df$X3 <- df$X2 / runif(n = 300, min = 1, max = 10)
    df$X2 <- 1 - (df$X1 + df$X3)
  }
  if (admixed == F) {
    df$X1 <-
      c(
        runif(n = 100, min = 0.9, max = 1),
        runif(n = 100, min = 0, max = 0.1),
        runif(n = 100, min = 0, max = 0.1)
      )
    df$X2 <- 1 - df$X1
    df$X3 <-
      df$X2 / c(
        runif(n = 100, min = 1, max = 10),
        runif(n = 100, min = 1, max = 1.05),
        runif(n = 100, min = 10, max = 1000)
      )
    df$X2 <- 1 - (df$X1 + df$X3)
  }
  
  long_d <-
    tidyr::pivot_longer(df, cols = !c(sample, breed, `k`))
  
  grp.labs <- paste("K =", k_)
  names(grp.labs) <- k_
  my_pal <- RColorBrewer::brewer.pal(n = 8, name = "Set1")
  
  ggplot(data = long_d, aes(x = sample, y = value, fill = name)) +
    geom_col(width = 1.1, color = NA) +
    facet_nested(
      k ~ breed,
      scales = "free_x",
      space = "free",
      switch = "both",
      labeller = labeller(k = grp.labs)
    ) +
    theme_classic() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = c(my_pal)) +
    theme(
      legend.position = "none",
      axis.text.x.top = element_text(
        angle = 90,
        hjust = 0,
        size = 4,
        vjust = 0.5
      ),
      panel.spacing = unit(0.2, "lines"),
      ggh4x.facet.nestline = element_line(linetype = 3),
      axis.ticks.length = unit(-1, "inch"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x.top = element_blank(),
      axis.title.y = element_blank()
    )
  
  ggsave(out_file,
         dpi = "print",
         width = width,
         height = height)
}
