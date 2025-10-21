library(algatr)    # [github::TheWangLab/algatr] v1.0.0 # %>% / mmrr_run / mmrr_table
library(adegenet)  # CRAN v2.1.10 # as.matrix
library(cowplot)   # CRAN v1.1.3 # plot_grid
library(tidyverse) # CRAN v2.0.0 


spp_labels <- function() {
  return(
    c(
      CD =  "_C. dactylon_<br>(Bermuda grass)",
      DS = "_D. sanguinalis_<br>(crabgrass)",
      EC = "_E. canadensis_<br>(horseweed)",
      LS = "_L. serriola_<br>(prickly lettuce)",
      PA = "_P. annua_<br>(bluegrass)",
      TO = "_T. officinale_<br>(dandelion)"
    )
  )
}


# First goal is to get the genetic distance matrix -- by site within city
# Need to correct the structure file to reflect this.
correct_structure_files <- function(spp_){
  infile_stem <- paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno")
  struct_ <-
    read_delim(paste0(infile_stem, ".structure")) %>% 
    rename(SampleName = `...1`, pop = `...2`)
  struct_ <-
    struct_ %>% 
    mutate(pop = as.numeric(as.factor(str_replace(struct_$SampleName, "..$", "")))) %>% 
    arrange(pop)
  write_delim(struct_, paste0(infile_stem, "_sitesaspops.structure"), delim = "\t")
  
  file <- readLines(paste0(infile_stem, "_sitesaspops.structure"))
  file[1] <- stringr::str_remove_all(file[1], "SampleName|pop")
  writeLines(file, paste0(infile_stem, "_sitesaspops.structure")) # This file will be opened in GenoDive
  # Chord distance is used in genodive
}

get_gen_dist_matrix <- function(spp_){
  infile_stem <- paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno")
  names_ <-
    read_delim(paste0(infile_stem, ".structure")) %>% 
    rename(SampleName = `...1`, pop = `...2`) %>% 
    mutate(pop = str_replace(SampleName, "..$", "")) %>%
    arrange(pop) %>% 
    distinct(pop) %>% 
    pull(pop)
  
  mat_ <- read_delim(paste0("output/genodive/", spp_, "_genetic_distance_bysite_rho.txt"))
  mat_ <- as.matrix(mat_ %>% dplyr::select(-`...1`), labels = T)
  colnames(mat_) <- rownames(mat_) <- names_
  mat_ <- as.dist(mat_)
  return(mat_)
}

get_geo_dist_matrix <- function(spp_){
  infile_stem <- paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno")
  names_ <-
    read_delim(paste0(infile_stem, ".structure")) %>% 
    rename(SampleName = `...1`, pop = `...2`) %>% 
    mutate(pop = str_replace(SampleName, "..$", "")) %>%
    arrange(pop) %>% 
    distinct(pop) %>% 
    dplyr::select(pop)
  
  # Get lat long from the site info
  # Make sure to include management type as otherwise site names aren't unique
  # Get the site info
  site_info <-
    read_csv("data/site_data_urban_cov.csv") %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".") %>%
    mutate(pop = paste0(spp_,".",site_abbv))

  d_geo <-
    names_ %>%
    left_join(site_info) %>%
    dplyr::select(pop, lat, long) %>%
    column_to_rownames(var = "pop")
  
  # Generate dissimilarity matrix
  Dgeo <- geodist::geodist(d_geo, measure = "geodesic")
  Dgeo <- Dgeo / 1000 # make kilometers
  Dgeo <- as.matrix(Dgeo, labels = T)
  colnames(Dgeo) <- rownames(Dgeo) <- rownames(d_geo)
  Dgeo <- as.dist(Dgeo)
  return(Dgeo)
}

get_env_dist_matrix_list <- function(spp_) {
  infile_stem <- paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno")
  names_ <-
    read_delim(paste0(infile_stem, ".structure")) %>% 
    rename(SampleName = `...1`, pop = `...2`) %>% 
    mutate(pop = str_replace(SampleName, "..$", "")) %>%
    arrange(pop) %>% 
    distinct(pop) %>% 
    dplyr::select(pop)
  
  # Get the site environmental info 
  # Make sure to include management type as otherwise site names aren't unique
  # Get the site info
  site_info <-
    read_csv("data/site_data_DUC_environvars.csv") 
  site_info <- 
    site_info %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  nlcd_urban_pct,
                  distance_to_city_center,
                  soiltemp_2.5cm_Apr_12pm,
                  soiltemp_2.5cm_Jul_12pm
                  ) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".") %>%
    mutate(pop = paste0(spp_,".",site_abbv)) %>% 
    rename(soiltemp_Apr = soiltemp_2.5cm_Apr_12pm,
           soiltemp_Jul = soiltemp_2.5cm_Jul_12pm)
  
  env_vars <- colnames(site_info)[2:5]
  
  env_list <- list()
  for(var in env_vars){
    d_env <-
      names_ %>%
      left_join(site_info) %>%
      distinct() %>%
      dplyr::select(pop, as.name(var)) %>%
      column_to_rownames(var = "pop")
    
    Denv <- dist(d_env)
    env_list[[var]] <- Denv
  }

  return(env_list)
}


do_MMRR <- function(spp_){
  Y <- get_gen_dist_matrix(spp_)
  X <- get_env_dist_matrix_list(spp_)
  X[["geodist"]] <- get_geo_dist_matrix(spp_)
  set.seed(444)
  results_full <- mmrr_run(Y, X, nperm = 9999, stdz = TRUE, model = "full")
  out_ <- mmrr_table(results_full, digits = 2, summary_stats = TRUE)
  out_ <- out_[["_data"]]
  out_$spp <- spp_
  return(out_)
}


make_mmrr_plot <- function(){
  
  if(!file.exists("output/MMRR/MMRR_total.csv")) {
    mmrr_ <- rbind(
      do_MMRR("CD"),
      do_MMRR("DS"),
      do_MMRR("EC"),
      do_MMRR("LS"),
      do_MMRR("PA"),
      do_MMRR("TO")
    )
    readr::write_csv(mmrr_, "output/MMRR/MMRR_total.csv")
  } else {
    mmrr_ <- readr::read_csv("output/MMRR/MMRR_total.csv")
  }
  
  if(!file.exists("output/MMRR/MMRR_total_overall_models.csv")) {
    p_vals_overall_models <-
      mmrr_ %>% filter(var %in% c("R-Squared:", "F-Statistic:", "F p-value:")) %>% 
      dplyr::select(var, estimate, spp) %>% 
      pivot_wider(names_from = var, values_from = estimate)
    readr::write_csv(p_vals_overall_models, "output/MMRR/MMRR_total_overall_models.csv")
  } else {
    p_vals_overall_models <- readr::read_csv("output/MMRR/MMRR_total_overall_models.csv")
  }
    
  mmrr_plot_ <-
    mmrr_ %>% filter(var != "F-Statistic:" &
                       var != "F p-value:" & var != "Intercept") %>%
    mutate(
      estimate = case_when(p > 0.05 ~ NA, TRUE ~ estimate),
      facet_ = case_when(
        var == "geodist" ~ "Distance",
        var %in% c(
          "distance_to_city_center",
          "nlcd_urban_pct",
          "soiltemp_Apr",
          "soiltemp_Jul"
        ) ~ "Environment"
      ),
      var = case_when(var == "R-Squared:" ~ "R-Squared", TRUE ~ var),
      estimate_rounded = round(estimate, digits = 2),
    ) %>%
    mutate(Species = spp_labels()[as.character(spp)])
  
  heat_plot <- function(in_dat, params = T) {
    plot_ <-
      ggplot() +
      geom_tile(
        data = in_dat ,
        aes(x = Species, y = var, fill = estimate)
      ) +
      geom_text(
        data = in_dat ,
        size = 2,
        aes(x = Species, y = var, label = estimate_rounded)
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 8)) +
      labs(x = NULL, y = NULL)
   
    if (params) {
      plot_ <- plot_ + scale_fill_gradient2(
        limit = c(-2, 1.5),
        high = "#641a80",
        low = "#fd6d21",
        mid = "white",
        midpoint = 0,
        na.value = "white"
      )
    }
    if (!params) {
      plot_ <- plot_ + scale_fill_gradient(
        limit = c(0, 1),
        high = "#20A387FF",
        low = "white"
      )
    }
    
    plot_ <-
      plot_ +
      ggnewscale::new_scale_fill() +
      geom_tile(data = in_dat %>% filter(spp != "EC"), aes(x = Species, y = var, fill = estimate)) +
      geom_text(
        data = in_dat %>% filter(spp != "EC"),
        size = 2,
        aes(x = Species, y = var, label = estimate_rounded)
      ) +
      scale_fill_gradient2(
        limit = c(-2, 1.5),
        high = "black",
        low = "black",
        mid = "white",
        midpoint = 0,
        na.value = "white"
      )
    
    return(plot_)
  }
  
  env_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(facet_ == "Environment"))  +
    labs(y = "Model variable") +
    scale_y_discrete(labels = c("Distance to city center", "% Impervious", "April soil temp.", "July soil temp."))  +
    theme(legend.title = element_blank(), axis.text.x = element_blank()) +
    guides(fill = FALSE)
  dist_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(facet_ == "Distance")) +
    scale_y_discrete(labels = c("Geographic distance")) +
    theme(legend.title = element_blank(), axis.text.x = element_blank()) +
    theme(legend.position = "none")
  r2_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(var == "R-Squared"), params = F) +
    guides(fill = FALSE) +
    theme(legend.title = element_blank(), 
          axis.text.x = ggtext::element_markdown())
    
  
  gg <- cowplot::plot_grid(
    env_plot_,
    dist_plot_,
    r2_plot_,
    ncol = 1,
    rel_heights = c(1,0.35, 0.85),
    align = "v"
  )
  
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, all lettering and symbols will be clear and easy to read,
  ggsave(paste0("figures/Fig7_MMRR/MMRR_total.png"),
         dpi = "print",
         height = 80,
         width = 112,
         units = "mm") 
}

do_MMRR_by_city <- function(spp_, city_){
  print(spp_); print(city_)
  
  # Get full matrices
  Y <- get_gen_dist_matrix(spp_)
  X <- get_env_dist_matrix_list(spp_)
  X[["geodist"]] <- get_geo_dist_matrix(spp_)
  
  # Filter matrices
  Y_city <- as.matrix(Y, labels = T)
  labels_to_keep <- labels(Y_city)[[1]][str_detect(labels(Y_city)[[1]], paste0(spp_,".",city_))]
  Y_city <- Y_city[rownames(Y_city) %in% labels_to_keep, colnames(Y_city) %in% labels_to_keep]
  Y_city <- as.dist(Y_city)
  
  X_city <- list()
  for(item in names(X)){
    X_temp_ <- X[[item]]
    X_temp_ <- as.matrix(X_temp_, labels = T)
    X_temp_ <- X_temp_[rownames(X_temp_) %in% labels_to_keep, colnames(X_temp_) %in% labels_to_keep]
    X_temp_ <- as.dist(X_temp_)
    X_city[[item]] <- X_temp_
  }
  set.seed(444)
  results_full <- mmrr_run(Y_city, X_city, nperm = 999, stdz = TRUE, model = "full")
  out_ <- mmrr_table(results_full, digits = 2, summary_stats = TRUE)
  out_ <- out_[["_data"]]
  out_$spp <- spp_
  out_$city <- city_
  return(out_)
}

make_mmrr_plot_bycity_ <- function(){
  mmrr_city_ <- rbind(
    do_MMRR_by_city("CD", "BA"),
    do_MMRR_by_city("CD", "LA"),
    do_MMRR_by_city("CD", "PX"),
    do_MMRR_by_city("DS", "BA"),
    do_MMRR_by_city("DS", "BO"),
    do_MMRR_by_city("DS", "MN"),
    do_MMRR_by_city("DS", "PX"),
    do_MMRR_by_city("EC", "BA"),
    do_MMRR_by_city("EC", "LA"),
    do_MMRR_by_city("EC", "PX"),
    do_MMRR_by_city("LS", "BA"),
    do_MMRR_by_city("LS", "BO"),
    do_MMRR_by_city("LS", "LA"),
    do_MMRR_by_city("LS", "MN"),
    do_MMRR_by_city("LS", "PX"),
    do_MMRR_by_city("PA", "BA"),
    do_MMRR_by_city("PA", "BO"),
    do_MMRR_by_city("PA", "LA"),
    do_MMRR_by_city("PA", "PX"),
    do_MMRR_by_city("TO", "BA"),
    do_MMRR_by_city("TO", "BO"),
    do_MMRR_by_city("TO", "LA"),
    do_MMRR_by_city("TO", "MN"),
    do_MMRR_by_city("TO", "PX")
  )
  
  readr::write_csv(mmrr_city_, "output/MMRR/MMRR_splitbycity.csv")
  
  p_vals_overall_models <-
    mmrr_city_ %>% filter(var %in% c("R-Squared:", "F-Statistic:", "F p-value:")) %>% 
    dplyr::select(var, estimate, spp, city) %>% 
    pivot_wider(names_from = var, values_from = estimate)
  
  p_vect <- pull(p_vals_overall_models, `F p-value:`)
  
  p_vals_overall_models$`p adjusted` <- round(c(
    p.adjust(p_vect[1:3], method = "BH"), #CD
    p.adjust(p_vect[4:7], method = "BH"), #DS
    p.adjust(p_vect[8:10], method = "BH"), #EC
    p.adjust(p_vect[11:15], method = "BH"), #LS
    p.adjust(p_vect[16:19], method = "BH"), #PA
    p.adjust(p_vect[20:24], method = "BH")  #TO
  ), 4)
  
  readr::write_csv(p_vals_overall_models, "output/MMRR/MMRR_splitbycity_overall_models.csv")
  
  mmrr_plot_city_ <- 
    mmrr_city_ %>% filter(var != "F-Statistic:" & var != "F p-value:" & var != "Intercept") %>% 
    mutate(
      estimate = case_when(p > 0.05 ~ NA, TRUE ~ estimate),
      facet_ = case_when(var == "geodist" ~ "Distance", var %in% c("distance_to_city_center", "nlcd_urban_pct", "soiltemp_Apr", "soiltemp_Jul") ~ "Environment"),
      var = case_when(var == "R-Squared:" ~ "R-Squared", TRUE ~ var),
      estimate_rounded = round(estimate, digits = 2),
      spp_city = paste0(spp, ".", city)
    )
  
  env_plot_ <- 
    ggplot(mmrr_plot_city_ %>% filter(facet_ == "Environment"), aes(x = spp_city, y = var)) +
    geom_tile(aes(fill = estimate)) +
    geom_text(aes(label = estimate_rounded)) +
    scale_fill_gradient2(limit = c(-2,1.5), low = "red", high =  "blue", mid = "white", midpoint = 0, na.value = "white") +
    scale_y_discrete(labels = c("Distance to city center", "% Urban cover", "April soil temp.", "July soil temp.")) +
    theme_classic() +
    labs(y = "Model variable", x = NULL) +
    theme(legend.title = element_blank(), axis.text.x = element_blank())
  
  dist_plot_ <- 
    ggplot(mmrr_plot_city_ %>% filter(facet_ == "Distance"), aes(x = spp_city, y = var)) +
    geom_tile(aes(fill = estimate)) +
    geom_text(aes(label = estimate_rounded)) +
    scale_fill_gradient2(limit = c(-2,1.5), low = "red", high =  "blue", mid = "white", midpoint = 0, na.value = "white") +
    scale_y_discrete(labels = c("Geographic distance")) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_blank()) +
    labs(x = NULL, y = NULL)
  
  r2_plot_ <- 
    ggplot(mmrr_plot_city_ %>% filter(var == "R-Squared"), aes(x = spp_city, y = var)) +
    geom_tile(aes(fill = estimate)) +
    geom_text(aes(label = estimate_rounded)) +
    scale_fill_gradient2(limit = c(-2,1.5), low = "red", high =  "blue", mid = "white", midpoint = 0, na.value = "white") +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Species + city combination", y = NULL) 
  
  gg <- cowplot::plot_grid(
    env_plot_,
    dist_plot_,
    r2_plot_,
    ncol = 1,
    rel_heights = c(1,0.35, 0.5),
    align = "v"
  )
  
  gg
  
  ggsave(paste0("figures/Fig7_MMRR/MMRR_splitbycity.png"),
         dpi = "print",
         height = 4,
         width = 12,
         units = "in") 
}


