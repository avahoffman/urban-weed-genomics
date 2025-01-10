library(magrittr) # CRAN v2.0.3 # %>%
library(adegenet) # CRAN v2.1.10 # as.matrix

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


