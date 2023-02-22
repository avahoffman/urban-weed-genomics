# # Archetype plots
#
library(here)
library(readr)
library(ggplot2)
library(cowplot)


make_pca_plot <- function(spp_, species_name) {
  pca_dat <-
    read_csv(paste0("output/pca/", spp_, "_IteratePopStructPCA.csv"))
  
  colors_ <- viridis::viridis(n = 5,
                              option = "H",
                              begin = 0.2)
  my_pal <- setNames(colors_,
                     c(
                       "Baltimore",
                       "Boston",
                       "Los Angeles",
                       "Minneapolis",
                       "Phoenix"
                     ))
  
  shapes_ <- c(21, 24, 22, 23, 25)
  
  shape_pal <- setNames(shapes_,
                        c(
                          "Baltimore",
                          "Boston",
                          "Los Angeles",
                          "Minneapolis",
                          "Phoenix"
                        ))
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      mapping = aes(
        x = PC1,
        y = PC2,
        fill = V2,
        shape = V2
      ),
      size = 3
    ) +
    scale_fill_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal) +
    guides(fill = guide_legend(title = ""),
           shape = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = "none",
          legend.text = element_text(size = 10)) +
    ggtitle(paste0("", species_name))
  
  # Save plot
  ggsave(
    paste0("figures/pca/pca_", spp_, ".png"),
    dpi = "print",
    width = 4,
    height = 4
  )
  
  return(gg)
  
}



plot_pcas <- function() {
  # Make plots
  p1 <- make_pca_plot("CD", "Bermuda grass")
  p2 <- make_pca_plot("DS", "crabgrass")
  p3 <- make_pca_plot("EC", "horseweed")
  p4 <- make_pca_plot("LS", "prickly lettuce")
  p5 <- make_pca_plot("PA", "bluegrass")
  p6 <- make_pca_plot("TO", "dandelion")
  
  legend <- get_legend(# create some space to the left of the legend
    p6 + theme(legend.position = "bottom",
               legend.direction = "horizontal"))
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    legend,
    align = 'v',
    axis = "l",
    hjust = -1,
    ncol = 1,
    rel_heights = c(1, 1, 1, 1, 1, 1, 0.5),
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here::here())
  ggsave(
    paste0("figures/pca/pca_all.png"),
    dpi = "print",
    width = 5,
    height = 18
  )
  
}

plot_pcas()
