# AA plots

library(tidyverse)
library(ggplot2)
library(ggh4x)

make_ADMIXTURE_plot <-
  function(spp,
           k_list,
           out_file,
           primary_dir = here::here(),
           width = 12,
           height = 3.5) {
    primary_dir <- here::here()
    temp_dir <-
      paste0(primary_dir,
             "/data/SNP_data/by-species/",
             spp,
             "/populations_20_pct")
    setwd(temp_dir)
    
    path_popmap <-
      paste0(primary_dir,
             "/data/SNP_data/by-species/",
             spp,
             "/popmap_",
             spp,
             ".txt")
  
    long_df <- data.frame()
    samp <-
      readr::read_tsv(path_popmap, col_names = c("sample", "city"))
    
    samp <- 
      samp %>% 
      separate(sample, into = c("spp","ct","site_abbv","management_type","id"), sep = "\\.", remove = F)
    site_info <- 
      read_csv(paste0(primary_dir,
                      "/data/macrosystems_urban_site_data_nlcd_pct_urban_by_site.csv")) %>% 
      dplyr::select(site_abbv, city, management_type, nlcd_urban_pct)
    samp <- 
      samp %>% left_join(site_info, by = c("city", "site_abbv", "management_type")) %>% 
      dplyr::select(sample, site_abbv, city, nlcd_urban_pct)
    for (k in k_list) {
      tmp_df <-
        readr::read_delim(paste0("AA.", k, ".Q"),
                          delim = ' ',
                          col_names = FALSE) %>% mutate(k = k)
      long_df_k <-
        tidyr::pivot_longer(cbind(samp, tmp_df), cols = !c(sample, city, `k`, site_abbv, nlcd_urban_pct))
      long_df <- rbind(long_df, long_df_k)
    }
    long_df$name <-
      factor(long_df$name, levels = rev(levels(as.factor(long_df$name))))
    long_df$sample <- 
      fct_reorder(factor(long_df$sample), long_df$nlcd_urban_pct, mean)
    
    grp.labs <- paste("K =", k_list)
    names(grp.labs) <- k_list
    my_pal <- RColorBrewer::brewer.pal(n = 8, name = "Set1")
    
    ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
      geom_col(width = 1.1, color = NA) +
      facet_nested(
        k ~ city,
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
    
    setwd(primary_dir)
    ggsave(out_file, dpi = "print", width = width, height = height)
  }

# TODO: Combine this with previous function
make_ADMIXTURE_plot_by_City <-
  function(spp,
           city,
           k_list,
           out_file,
           primary_dir = here::here(),
           width = 12,
           height = 3.5) {
    primary_dir <- here::here()
    temp_dir <-
      paste0(
        primary_dir,
        "/data/SNP_data/by-city/",
        spp,
        "_",
        city,
        "/populations_20_pct"
      )
    setwd(temp_dir)
    
    path_popmap <-
      paste0(primary_dir,
             "/data/SNP_data/by-city/",
             spp,
             "_",
             city,
             "/popmap_",
             spp,
             "_",
             city,
             ".txt")
    
    long_df <- data.frame()
    samp <-
      readr::read_tsv(path_popmap, col_names = c("sample", "city"))
    samp <- 
      samp %>% 
      separate(sample, into = c("spp","ct","site_abbv","type","id"), sep = "\\.", remove = F)
    site_info <- 
      read_csv(paste0(primary_dir,
                      "/data/macrosystems_urban_site_data_nlcd_pct_urban_by_site.csv"))
    samp <- 
      samp %>% left_join(site_info, by = c("city", "site_abbv")) %>% 
      dplyr::select(sample, site_abbv, city, nlcd_urban_pct)
    for (k in k_list) {
      tmp_df <-
        readr::read_delim(paste0("AA.", k, ".Q"),
                          delim = ' ',
                          col_names = FALSE) %>% mutate(k = k)
      long_df_k <-
        tidyr::pivot_longer(cbind(samp, tmp_df), cols = !c(sample, city, `k`, site_abbv, nlcd_urban_pct))
      long_df <- rbind(long_df, long_df_k)
    }
    long_df$name <-
      factor(long_df$name, levels = rev(levels(as.factor(long_df$name))))
    long_df$site_abbv <- 
      fct_reorder(factor(long_df$site_abbv), long_df$nlcd_urban_pct, mean)
    
    readr::write_csv(long_df, paste0("Q_Data_",
                                     spp,
                                     "_",
                                     city,
                                     ".csv"))
    
    grp.labs <- paste("K =", k_list)
    names(grp.labs) <- k_list
    my_pal <- RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")
    
    ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
      geom_col(width = 1.1, color = NA) +
      facet_nested(
        k ~ city + site_abbv,
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
          vjust = 0
        ),
        panel.spacing = unit(0, "lines"),
        #ggh4x.facet.nestline = element_line(linetype = 10),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y = element_blank()
      )
    
    setwd(primary_dir)
    ggsave(out_file, dpi = "print", width = width, height = height)
  }
