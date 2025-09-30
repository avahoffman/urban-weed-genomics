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


gather_pca_dat_w_pctimp <- function(spp_) {
  pca_dat <-
    read_csv(paste0("output/pca/", spp_, "_IteratePopStructPCA.csv")) %>%
    rename(sample = V1, city = V2)
  
  # Get the site info
  site_info <-
    read_csv(paste0(here(), "/data/site_data_urban_cov.csv")) %>%
    dplyr::select(site_abbv, city, management_type, nlcd_urban_pct) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))
  
  pca_dat <- pca_dat %>% tidyr::separate(
    sample,
    into = c("spp", "ct", "site_abbv", "management_type", "id"),
    sep = "\\.",
    remove = F
  ) %>%
    left_join(site_info, by = c("city", "site_abbv", "management_type")) 
  pca_dat$spp <- spp_
  return(pca_dat)
}


city_lab_vector <- function() {
  return(c(
    "Minneapolis-\nSaint Paul",
    "Boston",
    "Baltimore",
    "Los Angeles",
    "Phoenix"
  ))
}


spp_labels <- function() {
  return(
    c(
      CD =  "_C. dactylon_ (Bermuda grass)",
      DS = "_D. sanguinalis_ (crabgrass)",
      EC = "_E. canadensis_ (horseweed)",
      LS = "_L. serriola_ (prickly lettuce)",
      PA = "_P. annua_ (bluegrass)",
      TO = "_T. officinale_ (dandelion)"
    )
  )
}


# color PCA by city only (no urban labeling) - includes all species
# can optionally set labels_ to TRUE to sub for sample labels (help with troubleshooting/outliers)
make_pca_city_only_all <- function(labels_ = F) {
  pca_dat <- bind_rows(
    gather_pca_dat("CD"),
    gather_pca_dat("DS"),
    gather_pca_dat("EC"),
    gather_pca_dat("LS"),
    gather_pca_dat("PA"),
    gather_pca_dat("TO")
  )  
  pca_dat <-
    pca_dat %>%
    mutate(V2 = case_when(V2 == "Minneapolis" ~ "Minneapolis-\nSaint Paul", TRUE ~ V2)) %>% 
    mutate(V2 = factor(V2, levels = city_lab_vector()))
  
  # Palettes
  colors_ <- viridis::turbo(n = 5)
  my_pal <- setNames(colors_, city_lab_vector())
  shapes_ <- c(16, 4, 17, 3, 18)
  shape_pal <- setNames(shapes_, city_lab_vector())
  
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
    )
  
  if (labels_) {
    gg <- gg + geom_label(data = pca_dat %>% filter(V2 == "Boston"), aes(x = PC1, y = PC2, label = V1))
  }
  
  gg <- gg +
    facet_wrap(~ spp, scales = "free", labeller = labeller(spp = spp_labels())) +
    scale_color_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal) +
    guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) +
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.text = ggtext::element_markdown())
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required,
  # all lettering and symbols will be clear and easy to read,
  if (labels_) {
    ggsave(
      paste0("figures/Fig4_pca/pca_all_wlabels.png"),
      dpi = "screen",
      width = 1690,
      height = 1000,
      units = "mm",
      limitsize = F
    )
  } else {
    ggsave(
      paste0("figures/Fig4_pca/pca_city_all.png"),
      dpi = "print",
      width = 169,
      height = 100,
      units = "mm"
    )
  }
}


# color PCA by %impervious only (no city labeling) - includes all species
make_pca_pctimp_only_all <- function() {
  pca_dat <- bind_rows(
    gather_pca_dat_w_pctimp("CD"),
    gather_pca_dat_w_pctimp("CD"),
    gather_pca_dat_w_pctimp("DS"),
    gather_pca_dat_w_pctimp("EC"),
    gather_pca_dat_w_pctimp("LS"),
    gather_pca_dat_w_pctimp("PA"),
    gather_pca_dat_w_pctimp("TO")
  )  
  pca_dat <-
    pca_dat %>%
    mutate(city = case_when(city == "Minneapolis" ~ "Minneapolis-\nSaint Paul", TRUE ~ city)) %>% 
    mutate(city = factor(city, levels = city_lab_vector()))
  
  # Palettes
  colors_ <- viridis::turbo(n = 5)
  my_pal <- setNames(colors_, city_lab_vector())
  shapes_ <- c(16, 4, 17, 3, 18)
  shape_pal <- setNames(shapes_, city_lab_vector())
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      mapping = aes(x = PC1, y = PC2, fill = nlcd_urban_pct, color = nlcd_urban_pct),
      shape = 21,
      size = 1.5
    ) +
    scale_fill_viridis_c(option = "A", limits = c(0, 100)) +
    scale_color_viridis_c(option = "A", limits = c(0, 100)) +
    guides(fill = guide_legend(title = "% Impervious"), color = guide_legend(title = "% Impervious")) +
    theme_bw() +
    theme(legend.text = element_text(size = 10)) +
    facet_wrap(~ spp, scales = "free", labeller = labeller(spp = spp_labels())) +
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.text = ggtext::element_markdown())
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required,
  # all lettering and symbols will be clear and easy to read,
  ggsave(
      paste0("figures/Fig4_pca/pca_pctimp_all.png"),
      dpi = "print",
      width = 169,
      height = 100,
      units = "mm"
    )
}


make_pca_city_and_pctimp_all <- function() {
  pca_dat <- bind_rows(
    gather_pca_dat_w_pctimp("CD"),
    gather_pca_dat_w_pctimp("CD"),
    gather_pca_dat_w_pctimp("DS"),
    gather_pca_dat_w_pctimp("EC"),
    gather_pca_dat_w_pctimp("LS"),
    gather_pca_dat_w_pctimp("PA"),
    gather_pca_dat_w_pctimp("TO")
  )  
  pca_dat <-
    pca_dat %>%
    mutate(city = case_when(city == "Minneapolis" ~ "Minneapolis-\nSaint Paul", TRUE ~ city)) %>% 
    mutate(city = factor(city, levels = city_lab_vector())) %>% 
    group_by(spp) %>% 
    mutate(median_pctimp = median(nlcd_urban_pct)) %>% 
    mutate(is_urban = nlcd_urban_pct >= median_pctimp)
  
  # Palettes
  colors_ <- viridis::turbo(n = 5)
  my_pal <- setNames(colors_, city_lab_vector())
  shapes_ <- c(16, 1)
  shape_pal <- setNames(shapes_, c(T, F))
  
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      alpha = 0.8,
      mapping = aes(
        x = PC1,
        y = PC2,
        color = city,
        shape = is_urban
      )
    ) +
    facet_wrap(~ spp, scales = "free", labeller = labeller(spp = spp_labels())) +
    scale_color_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal, labels = c("Low", "High")) +
    guides(color = guide_legend(title = "City"), shape = guide_legend(title = "% Impervious")) +
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.text = ggtext::element_markdown())
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required,
  # all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/Fig4_pca/pca_city_and_pctimp_all.png"),
    dpi = "print",
    width = 169,
    height = 100,
    units = "mm"
  )
}


make_urban_and_PC1_PC2 <- function() {
  pca_dat <- bind_rows(
    gather_pca_dat("CD"),
    gather_pca_dat("DS"),
    gather_pca_dat("EC"),
    gather_pca_dat("LS"),
    gather_pca_dat("PA"),
    gather_pca_dat("TO")
  )  %>%
    mutate(
      spp = case_when(
        spp == "CD" ~ "Bermuda grass",
        spp == "DS" ~ "crabgrass",
        spp == "EC" ~ "horseweed",
        spp == "LS" ~ "prickly lettuce",
        spp == "PA" ~ "bluegrass",
        spp == "TO" ~ "dandelion"
      )
    )
  
  pca_dat$spp <- factor(
    pca_dat$spp,
    levels = c(
      "Bermuda grass",
      "crabgrass",
      "horseweed",
      "prickly lettuce",
      "bluegrass",
      "dandelion"
    )
  )
  
  # Get the site info
  site_info <-
    read_csv(paste0(here(), "/data/site_data_urban_cov.csv")) %>%
    dplyr::select(site_abbv, city, management_type, nlcd_urban_pct) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))
  
  pca_dat <- pca_dat %>%
    rename(sample = V1, city = V2) %>%
    tidyr::separate(
      sample,
      into = c("spp_abbv", "ct", "site_abbv", "management_type", "id"),
      sep = "\\.",
      remove = F
    ) %>%
    left_join(site_info, by = c("city", "site_abbv", "management_type")) %>%
    mutate(is_urban = nlcd_urban_pct > 40) %>%
    drop_na(is_urban)
  
  pca_dat <-
    pca_dat %>%
    mutate(city = case_when(city == "Minneapolis" ~ "Minneapolis-\nSaint Paul", TRUE ~ city))
  
  pca_dat$city <- factor(
    pca_dat$city,
    levels = c(
      "Minneapolis-\nSaint Paul",
      "Boston",
      "Baltimore",
      "Los Angeles",
      "Phoenix"
    )
  )
  
  colors_ <- viridis::turbo(n = 5)
  my_pal <- setNames(
    colors_,
    c(
      "Minneapolis-\nSaint Paul",
      "Boston",
      "Baltimore",
      "Los Angeles",
      "Phoenix"
    )
  )
  
  shapes_ <- c(16, 1)
  
  shape_pal <- setNames(shapes_, c(T, F))
  
  gg <- ggplot() +
    geom_point(
      data = pca_dat,
      alpha = 0.8,
      mapping = aes(
        x = PC1,
        y = PC2,
        color = city,
        shape = is_urban
      )
    ) +
    facet_wrap(~ spp, scales = "free") +
    scale_color_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal, labels = c("Rural", "Urban")) +
    guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) +
    theme_bw() +
    theme(text = element_text(size = 8))
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required,
  # all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/pca/pc1_pc2_with_urban.png"),
    dpi = "print",
    width = 169,
    height = 100,
    units = "mm"
  )
}