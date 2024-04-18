#' Test for isolation by distance
#' Mantel test for correlation between two distance matrices
#' https://www.uwyo.edu/dbmcd/molmark/lect06/lect6.html has a nice summary of
#' distance metrics in genomics.
#'
#' Also found the adegenet docs helpful! 
#' https://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf

library(adegenet)
library(tidyverse)


do_ibe <- function(spp_, env_var_to_use) {
  # Load in genind format data
  gen_ <-
    readRDS(paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno_genind.rds"))
  
  # Extract site names (Will use these as pops)
  cities_ <- pop(gen_)
  pop(gen_) <-
    paste0(pop(gen_), ".", str_replace(rownames(gen_$tab), "......", "") %>% str_replace("..$", ""))
  pop(gen_)
  site_city <-
    data.frame(city_abbv = cities_, site_abbv = pop(gen_))
  
  # Get the site environmental info
  # Make sure to include management type as otherwise site names aren't unique
  # Get the site info
  site_info <-
    read_csv("data/site_data_DUC_environvars.csv") 
  site_info <- site_info %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  as.name(env_var_to_use)) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".")
  
  # Clean and order geographic data
  d_env <-
    site_city %>%
    left_join(site_info) %>%
    distinct() %>%
    dplyr::select(site_abbv, as.name(env_var_to_use)) %>%
    column_to_rownames(var = "site_abbv")
  
  # Convert to genpop
  gen_genpop <- genind2genpop(gen_)
  
  # Generate dissimilarity matrices
  Dgen <- dist.genpop(gen_genpop, method = 2)
  Denv <- dist(d_env)
  
  # Run monte-carlo test
  message("Running mantel test...")
  ibe <- mantel.randtest(Dgen, Denv, nrepet = 9999)
  return(list(ibe, Denv, Dgen))
}

make_matrix_df <- function(spp_, Denv, Dgen) {
  Denv_df <- reshape2::melt(as.matrix(Denv)) %>% rename(Denv = value)
  Dgen_df <-
    reshape2::melt(as.matrix(Dgen)) %>% rename(Dgen = value)
  df <- full_join(Denv_df, Dgen_df)
  df$spp <- spp_
  
  return(df)
}


extract_ibe_stats_and_plots <- function(city_level = FALSE, city_choice = "BA", env_var_to_use = "nlcd_urban_pct") {
  CD <- do_ibe(spp_ = "CD", env_var_to_use = env_var_to_use)
  CD_ <- CD[[1]]
  CD_ <- c(CD_$obs, CD_$alter, CD_$rep, CD_$expvar, CD_$pvalue)
  
  DS <- do_ibe(spp_ = "DS", env_var_to_use = env_var_to_use)
  DS_ <- DS[[1]]
  DS_ <- c(DS_$obs, DS_$alter, DS_$rep, DS_$expvar, DS_$pvalue)
  
  EC <- do_ibe(spp_ = "EC", env_var_to_use = env_var_to_use)
  EC_ <- EC[[1]]
  EC_ <- c(EC_$obs, EC_$alter, EC_$rep, EC_$expvar, EC_$pvalue)
  
  LS <- do_ibe(spp_ = "LS", env_var_to_use = env_var_to_use)
  LS_ <- LS[[1]]
  LS_ <- c(LS_$obs, LS_$alter, LS_$rep, LS_$expvar, LS_$pvalue)
  
  PA <- do_ibe(spp_ = "PA", env_var_to_use = env_var_to_use)
  PA_ <- PA[[1]]
  PA_ <- c(PA_$obs, PA_$alter, PA_$rep, PA_$expvar, PA_$pvalue)
  
  TO <- do_ibe(spp_ = "TO", env_var_to_use = env_var_to_use)
  TO_ <- TO[[1]]
  TO_ <- c(TO_$obs, TO_$alter, TO_$rep, TO_$expvar, TO_$pvalue)
  
  out <- data.frame(rbind(CD_, DS_, EC_, LS_, PA_, TO_))
  out$spp <- c("Bermuda grass (CD)", "crabgrass (DS)", "horseweed (EC)", "prickly lettuce (LS)", "bluegrass (PA)", "dandelion (TO)")
  out$env <- env_var_to_use
  colnames(out) <-
    c(
      "Observation",
      "Hypothesis",
      "Reps",
      "Std.Obs",
      "Expectation",
      "Variance",
      "p-value",
      "Species",
      "Environmental variable"
    )
  out <- relocate(out, "Species", .before = "Observation")
  # readr::write_csv(out, paste0("output/isolation-by-envt-", env_var_to_use,"-mantel-test.csv"))
  
  df <- rbind(
    make_matrix_df("CD", CD[[2]], CD[[3]]),
    make_matrix_df("DS", DS[[2]], DS[[3]]),
    make_matrix_df("EC", EC[[2]], EC[[3]]),
    make_matrix_df("LS", LS[[2]], LS[[3]]),
    make_matrix_df("PA", PA[[2]], PA[[3]]),
    make_matrix_df("TO", TO[[2]], TO[[3]])
  )
  
  df <- df %>% 
    mutate(spp = case_when(
      spp == "CD" ~ "Bermuda grass",
      spp == "DS" ~ "crabgrass",
      spp == "EC" ~ "horseweed",
      spp == "LS" ~ "prickly lettuce",
      spp == "PA" ~ "bluegrass",
      spp == "TO" ~ "dandelion"
    ))
  df$spp <- factor(
    df$spp,
    levels = c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")
  )
  
  site_info <-
    read_csv("data/site_data_DUC_environvars.csv") %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  as.name(env_var_to_use)) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".", remove = F)
  
  if(city_level){
    only_one_city_sites <- site_info %>% filter(city_abbv == city_choice) %>% pull(site_abbv)
    df <- df %>% 
      filter(Var1 %in% only_one_city_sites, Var2 %in% only_one_city_sites) %>% 
      filter(Denv > 0)
  } else {
    city_choice <- "all"
  }
  
  gg <- ggplot(data = df, aes(Denv, Dgen)) +
    geom_density_2d_filled() +
    facet_wrap(~spp) +
    geom_point(shape = 1, alpha = 0.1) +
    scale_fill_manual(values = c("white",  rev(viridisLite::magma(n = 13)))) +
    theme_classic() +
    theme(legend.position = "none") +
    #geom_smooth(formula = y~x, se = FALSE, linetype = 2) +
    labs(
      x = paste0("Environmental distance (", env_var_to_use, ")"),
      y = "Genetic distance"
    )
  
  gg
  
  ggsave(paste0("figures/IBE/isolation_by_envt_citylevel", city_level, "_", city_choice, "_", env_var_to_use,".png"),
         dpi = "print",
         height = 5,
         width = 8,
         units = "in") 
}

