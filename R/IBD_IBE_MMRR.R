library(algatr)
library(adegenet)
library(tidyverse)
library(cowplot)

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
  
  mat_ <- read_delim(paste0("output/genodive/", spp_, "_genetic_distance_bysite.txt"))
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
  mmrr_ <- rbind(
    do_MMRR("CD"),
    do_MMRR("DS"),
    do_MMRR("EC"),
    do_MMRR("LS"),
    do_MMRR("PA"),
    do_MMRR("TO")
  )
  
  readr::write_csv(mmrr_, "output/MMRR/MMRR_total.csv")
    
  mmrr_plot_ <- 
    mmrr_ %>% filter(var != "F-Statistic:" & var != "F p-value:" & var != "Intercept") %>% 
    mutate(
      estimate = case_when(p > 0.05 ~ NA, TRUE ~ estimate),
      facet_ = case_when(var == "geodist" ~ "Distance", var %in% c("distance_to_city_center", "nlcd_urban_pct", "soiltemp_Apr", "soiltemp_Jul") ~ "Environment"),
      var = case_when(var == "R-Squared:" ~ "R-Squared", TRUE ~ var),
      estimate_rounded = round(estimate, digits = 2),
      )  %>% 
    mutate(spp = case_when(
      spp == "CD" ~ "Bermuda grass",
      spp == "DS" ~ "crabgrass",
      spp == "EC" ~ "horseweed",
      spp == "LS" ~ "prickly lettuce",
      spp == "PA" ~ "bluegrass",
      spp == "TO" ~ "dandelion"
    ))
  mmrr_plot_$spp <- factor(
    mmrr_plot_$spp,
    levels = c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")
  )
  
  heat_plot <- function(in_dat) {
    plot_ <- 
      ggplot() +
      geom_tile(data = in_dat %>% filter(spp != "bluegrass", spp != "dandelion"),
                aes(x = spp, y = var, fill = estimate)) +
      geom_text(data = in_dat %>% filter(spp != "bluegrass", spp != "dandelion"),
                size = 2,
                aes(x = spp, y = var, label = estimate_rounded)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 8)
      ) +
      labs(x = NULL, y = NULL) +
      scale_fill_gradient2(limit = c(-2,1.5), high = "#641a80", low = "#fd6d21", mid = "white", midpoint = 0, na.value = "white") +
      ggnewscale::new_scale_fill() +
      geom_tile(data = in_dat %>% filter(spp %in% c("bluegrass", "dandelion")), aes(x = spp, y = var, fill = estimate)) +
      geom_text(data = in_dat %>% filter(spp %in% c("bluegrass", "dandelion")),
                size = 2,
                aes(x = spp, y = var, label = estimate_rounded)) +
      scale_fill_gradient2(limit = c(-2,1.5), high = "black", low = "black", mid = "white", midpoint = 0, na.value = "white") 
    
    return(plot_)
  }
  
  env_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(facet_ == "Environment"))  +
    labs(y = "Model variable") +
    scale_y_discrete(labels = c("Distance to city center", "% Urban cover", "April soil temp.", "July soil temp."))  +
    theme(legend.title = element_blank(), axis.text.x = element_blank()) +
    guides(fill = FALSE)
  dist_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(facet_ == "Distance")) +
    scale_y_discrete(labels = c("Geographic distance")) +
    theme(legend.title = element_blank(), axis.text.x = element_blank()) +
    theme(legend.position = "none")
  r2_plot_ <- 
    heat_plot(mmrr_plot_ %>% filter(var == "R-Squared")) +
    theme(legend.position = "none")
  
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
  ggsave(paste0("figures/MMRR/MMRR_total.png"),
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
  results_full <- mmrr_run(Y_city, X_city, nperm = 9999, stdz = TRUE, model = "full")
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
    select(var, estimate, spp, city) %>% 
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
  
  ggsave(paste0("figures/MMRR/MMRR_splitbycity.png"),
         dpi = "print",
         height = 4,
         width = 12,
         units = "in") 
}

get_gen_dist_matrix_Fst <- function(spp_){
  infile_stem <- paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno")
  names_ <-
    read_delim(paste0(infile_stem, ".structure")) %>% 
    rename(SampleName = `...1`, pop = `...2`) %>% 
    mutate(pop = str_replace(SampleName, "..$", "")) %>%
    arrange(pop) %>% 
    distinct(pop) %>% 
    pull(pop)
  
  mat_ <- read_delim(paste0("output/genodive/", spp_, "_genetic_distance_bysite_Fst.txt"))
  mat_ <- as.matrix(mat_ %>% dplyr::select(-`...1`))
  colnames(mat_) <- names_
  rownames(mat_) <- names_
  mat_ <- as.dist(mat_)
  
  # Set all negative Fst to zero
  mat_[mat_ < 0] <- 0
  
  return(mat_)
}

make_scatter_Fst_dist <- function(spp_){
  Dgen <- get_gen_dist_matrix_Fst(spp_)
  Dgeo <- get_geo_dist_matrix(spp_)
  
  gen_df <- 
    reshape2::melt(as.matrix(Dgen)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Dgen = value) %>% 
    mutate(city1 = substr(Var1, start = 4, stop = 5),
           city2 = substr(Var2, start = 4, stop = 5),
           type = case_when(
             city1 == city2 ~ "Within-city",
             TRUE ~ "Among-city"
           ))
  
  geo_df <- 
    reshape2::melt(as.matrix(Dgeo)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Dgeo = value) %>% 
    mutate(log_geodist = log10(Dgeo))
  
  gen_geo_df <- 
    gen_df %>% 
    left_join(geo_df, by = join_by(Var1, Var2))
  
  gen_geo_df$spp <- spp_

  return(gen_geo_df)
}

get_mantel_geo_Fst <- function(spp_){
  Dgen <- get_gen_dist_matrix_Fst(spp_)
  Dgeo <- get_geo_dist_matrix(spp_)
  
  set.seed(444)
  mantel_ <- mantel.randtest(Dgeo, Dgen, nrepet = 9999)
  
  mantel_df <- tibble(
    spp = spp_,
    obs = mantel_$obs,
    p = mantel_$pvalue
  )
  
  return(mantel_df)
}

make_scatter_Fst_dist_plot <- function() {
  spp_dat_comb <- 
    rbind(
      make_scatter_Fst_dist("CD"),
      make_scatter_Fst_dist("DS"),
      make_scatter_Fst_dist("EC"),
      make_scatter_Fst_dist("LS"),
      make_scatter_Fst_dist("PA"),
      make_scatter_Fst_dist("TO")
    )
  
  spp_mantels <-
    rbind(
      get_mantel_geo_Fst("CD"),
      get_mantel_geo_Fst("DS"),
      get_mantel_geo_Fst("EC"),
      get_mantel_geo_Fst("LS"),
      get_mantel_geo_Fst("PA"),
      get_mantel_geo_Fst("TO")
    )
  spp_mantels$spp <- factor(c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion"))
  spp_mantels$type <- "Within-city"
  
  colors_ <- viridis::viridis(2, begin = 0.3, end = 0.8, option = "A") #c("grey10", "grey50")
  my_pal <- setNames(colors_, c("Among-city", "Within-city"))
  shapes_ <- c(4, 1)
  shape_pal <- setNames(shapes_, c("Among-city", "Within-city"))
  
  spp_dat <- spp_dat_comb %>% 
    left_join(spp_mantels) %>% 
    mutate(spp = case_when(
      spp == "CD" ~ "Bermuda grass",
      spp == "DS" ~ "crabgrass",
      spp == "EC" ~ "horseweed",
      spp == "LS" ~ "prickly lettuce",
      spp == "PA" ~ "bluegrass",
      spp == "TO" ~ "dandelion"
    ))
  spp_dat$spp <- factor(
    spp_dat$spp,
    levels = c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")
  )
  
  gg <-
    ggplot(data = spp_dat, aes(
      x = log_geodist,
      y = Dgen,
      color = type,
      shape = type
    )) +
    facet_wrap(~spp, ) +
    geom_point() +
    geom_text(
      data    = spp_mantels,
      mapping = aes(x = -0.75, y = 0.55, label = paste0("r = ",round(obs, 3))),
      hjust   = 0,
      color = "black",
      size = 2
    ) +
    geom_text(
      data    = spp_mantels,
      mapping = aes(x = -0.75, y = 0.50, label = paste0("p = ",round(p, 4))),
      hjust   = 0,
      color = "black",
      size = 2
    ) +
    scale_color_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal) +
    guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) +
    theme_bw() +
    labs(x = "Log Euclidean distance (km)", y = expression("F"["ST"])) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      text = element_text(size = 8)
    )
    
  gg
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/MMRR/Fst_geo_all.png"),
    dpi = "print",
    width = 169,
    height = 100,
    units = "mm"
  )
}
