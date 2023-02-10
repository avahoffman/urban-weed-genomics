# Archetype plots

library(here)
library(tidyr)
library(readr)
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggh4x)
library(cowplot)

make_archetype_plot <-
  function(spp,
           species_name,
           k_list,
           width = 12,
           height = 3.5) {
    subdir <- paste0(here(), "/SNP_data/", spp)
    setwd(paste0(subdir, "/populations_20_pct"))
    
    # Get the site info
    site_info <-
      read_csv(paste0(here(),
                             "/data/site_data_urban_cov.csv")) %>%
      dplyr::select(site_abbv, city, management_type, nlcd_urban_pct) %>% 
      mutate(site_abbv = str_replace(site_abbv, " ", "_"))
    
    # Get samples and spread out the name to get more details
    # Join to site info to get urban cover
    samp <-
      read_tsv(paste0(subdir, "/popmap_", spp, ".txt"),
                      col_names = c("sample", "city")) %>%
      tidyr::separate(
        sample,
        into = c("spp", "ct", "site_abbv", "management_type", "id"),
        sep = "\\.",
        remove = F
      ) %>% 
      left_join(site_info,
                             by = c("city", "site_abbv", "management_type")) %>%
      dplyr::select(sample, site_abbv, city, nlcd_urban_pct, id, management_type)
    
    # Create empty df
    long_df <- data.frame()
    
    # Pull in archetypal analysis results by k
    for (k in k_list) {
      tmp_df <-
        read_delim(paste0("AA.", k, ".Q"),
                          delim = ' ',
                          col_names = FALSE) %>% mutate(k = k)
      long_df_k <-
        pivot_longer(
          cbind(samp, tmp_df),
          cols = !c(
            sample,
            city,
            `k`,
            site_abbv,
            id,
            management_type,
            nlcd_urban_pct
          )
        )
      long_df <- rbind(long_df, long_df_k)
    }
    
    # Create a new column with cover percent, make a factor, and sort
    long_df <-
      long_df %>%
      mutate(site_n_cov = paste0(round(nlcd_urban_pct, 0), "%")) %>% 
      mutate(site_n_cov = fct_reorder(site_n_cov, nlcd_urban_pct)) %>% 
      arrange(city, nlcd_urban_pct)
    
    # Create a set of labels for the x axis
    x_lbl_ <-
      long_df %>% 
      select(sample, city, site_n_cov, nlcd_urban_pct) %>% 
      unique()
    
    # Keep only the first label in the list to keep things tidy
    x_lbl <- x_lbl_ %>% 
      mutate(dupe = x_lbl_ %>% select(-sample) %>% duplicated()) %>% 
      mutate(site_n_cov = as.character(site_n_cov)) %>% 
        mutate(x = case_when(
        dupe == FALSE ~ site_n_cov,
        TRUE ~ "-"
      )) %>% 
      mutate(x = as_factor(x))
    
    # Create labels and set palette
    grp.labs <- paste("K =", k_list)
    names(grp.labs) <- k_list
    # my_pal <-
    #   rev(RColorBrewer::brewer.pal(n = max(long_df$k), name = "Set3"))
    colors_ <- viridis::viridis(n = 10, option = "H", begin = 0)
    my_pal <- setNames(colors_, 
                       c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"))
    
    gg <-
      ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
      geom_col(width = 1, color = NA) +
      facet_nested(
        k ~ city,
        scales = "free_x",
        space = "free",
        switch = "both",
        labeller = labeller(k = grp.labs)
      ) +
      theme_classic() +
      facetted_pos_scales(
        x = list(
          city == "Baltimore" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Baltimore"), ]$x),
          city == "Boston" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Boston"), ]$x),
          city == "Los Angeles" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Los Angeles"), ]$x),
          city == "Minneapolis" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Minneapolis"), ]$x),
          city == "Phoenix" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Phoenix"), ]$x)
        )
      ) +
      scale_fill_manual(values = c(my_pal)) +
      theme(
        legend.position = "none",
        axis.text.x.top = element_text(
          angle = 90,
          hjust = 0,
          size = 5,
          vjust = 0.5
        ),
        #axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        ggh4x.facet.nestline = element_line(linetype = 3),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y = element_blank()
      ) +
      ggtitle(paste0("  ", species_name))
    
    setwd(here())
    ggsave(
      paste0("figures/genetics/Archetypes_", species_name, ".png"),
      dpi = "print",
      width = width,
      height = height
    )
    return(gg)
  }


make_archetype_multi_plot <- function(width = 12,
                                      height = 12) {
  p1 <- make_archetype_plot(spp = "CD",
                            species_name = "Bermuda grass",
                            k_list = c(3, 6))
  p2 <- make_archetype_plot(spp = "DS",
                            species_name = "crabgrass",
                            k_list = c(4, 8))
  p3 <- make_archetype_plot(spp = "EC",
                            species_name = "horseweed",
                            k_list = c(3, 6))
  p4 <- make_archetype_plot(spp = "LS",
                            species_name = "prickly lettuce",
                            k_list = c(5, 10))
  p5 <- make_archetype_plot(spp = "PA",
                            species_name = "bluegrass",
                            k_list = c(5, 10))
  p6 <- make_archetype_plot(spp = "TO",
                            species_name = "dandelion",
                            k_list = c(5, 10))
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    align = 'v',
    axis = "l",
    #hjust = -1,
    ncol = 1,
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here())
  ggsave(
    paste0("figures/genetics/Archetypes_ALL.png"),
    dpi = "print",
    width = width,
    height = height
  )
}
