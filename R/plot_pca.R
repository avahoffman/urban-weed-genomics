# # pca plots
#
library(here)
library(readr)
library(ggplot2)
library(cowplot)
library(tidyverse)


gather_pca_dat <- function(spp_) {
  pca_dat <-
    read_csv(paste0("output/pca/", spp_, "_IteratePopStructPCA.csv"))
  pca_dat$spp <- spp_
  return(pca_dat)
}


make_pca_plot2 <- function(){
  pca_dat <- bind_rows(
    gather_pca_dat("CD"),
    gather_pca_dat("DS"),
    gather_pca_dat("EC"),
    gather_pca_dat("LS"),
    gather_pca_dat("PA"),
    gather_pca_dat("TO")
  )  %>% 
    mutate(spp = case_when(
      spp == "CD" ~ "Bermuda grass",
      spp == "DS" ~ "crabgrass",
      spp == "EC" ~ "horseweed",
      spp == "LS" ~ "prickly lettuce",
      spp == "PA" ~ "bluegrass",
      spp == "TO" ~ "dandelion"
    ))
  pca_dat$spp <- factor(
    pca_dat$spp,
    levels = c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")
  )
  pca_dat$V2<- factor(
    pca_dat$V2,
    levels = c("Minneapolis", "Boston", "Baltimore", "Los Angeles", "Phoenix")
  )
  
  colors_ <- viridis::viridis(5, end = 0.8, option = "A") 
  my_pal <- setNames(colors_,
                     c(
                       "Minneapolis",
                       "Boston",
                       "Baltimore",
                       "Los Angeles",
                       "Phoenix"
                     ))
  
  shapes_ <- c(16, 4, 17, 3, 18)
  
  shape_pal <- setNames(shapes_,
                        c(
                          "Minneapolis",
                          "Boston",
                          "Baltimore",
                          "Los Angeles",
                          "Phoenix"
                        ))
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      alpha = 0.8,
      mapping = aes(
        x = PC1,
        y = PC2,
        color = V2,
        shape = V2
      )
    ) +
    facet_wrap(~spp, scales = "free") +
    scale_color_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal) +
    guides(color = guide_legend(title = ""),
           shape = guide_legend(title = "")) +
    theme_bw() +
    theme(text = element_text(size = 8))
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/pca/pca_all.png"),
    dpi = "print",
    width = 169,
    height = 100,
    units = "mm"
  )
}


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
      size = 2
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
    align = 'vh',
    #hjust = -1,
    ncol = 3,
    rel_heights = c(1, 1, 1),
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot_legend <- plot_grid(
    mega_plot,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.1)
  )
  mega_plot_legend
  setwd(here::here())
  ggsave(
    paste0("figures/pca/pca_all.png"),
    dpi = "print",
    width = 7,
    height = 5
  )
  
}


make_pca_plot_urbanness <- function(spp_, species_name) {
  pca_dat <-
    read_csv(paste0("output/pca/", spp_, "_IteratePopStructPCA.csv")) %>% 
    rename(sample = V1, city = V2)
  
  # Get the site info
  site_info <-
    read_csv(paste0(here(),
                    "/data/site_data_urban_cov.csv")) %>%
    dplyr::select(site_abbv, city, management_type, nlcd_urban_pct) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))
  
  pca_dat <- pca_dat %>% tidyr::separate(
    sample,
    into = c("spp", "ct", "site_abbv", "management_type", "id"),
    sep = "\\.",
    remove = F
  ) %>%
    left_join(site_info,
              by = c("city", "site_abbv", "management_type")) %>% 
    mutate(nlcd_urban_pct = case_when(
      nlcd_urban_pct > 50 ~ 100,
      nlcd_urban_pct <= 50 ~ 0
    ))
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      mapping = aes(
        x = PC1,
        y = PC2,
        fill = nlcd_urban_pct,
      ),
      shape = 21,
      size = 2
    ) +
    scale_fill_viridis_c(option = "A", limits = c(0, 100)) +
    guides(fill = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = "none",
          legend.text = element_text(size = 10)) +
    ggtitle(paste0("", species_name))
  
  # Save plot
  ggsave(
    paste0("figures/pca/pca_urbanness_", spp_, ".png"),
    dpi = "print",
    width = 4,
    height = 4
  )
  
  return(gg)
  
}


plot_pcas_urbanness <- function() {
  # Make plots
  p1 <- make_pca_plot_urbanness("CD", "Bermuda grass")
  p2 <- make_pca_plot_urbanness("DS", "crabgrass")
  p3 <- make_pca_plot_urbanness("EC", "horseweed")
  p4 <- make_pca_plot_urbanness("LS", "prickly lettuce")
  p5 <- make_pca_plot_urbanness("PA", "bluegrass")
  p6 <- make_pca_plot_urbanness("TO", "dandelion")
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    align = 'vh',
    #hjust = -1,
    ncol = 3,
    rel_heights = c(1, 1, 1),
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here::here())
  ggsave(
    paste0("figures/pca/pca_urbanness_all.png"),
    dpi = "print",
    width = 7,
    height = 5
  )
  
  return(mega_plot)
}