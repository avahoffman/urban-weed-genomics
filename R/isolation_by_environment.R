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
    dplyr::select(city_abbv, site_abbv, as.name(env_var_to_use)) %>%
    column_to_rownames(var = "site_abbv")
  
  # Convert to genpop
  gen_genpop <- genind2genpop(gen_)
  
  # Generate dissimilarity matrices
  Dgen <- dist.genpop(gen_genpop, method = 2)
  Denv <- dist(d_env %>% select(-city_abbv))
  
  # Run monte-carlo mantel test
  message("Running mantel test...")
  ibe <- mantel.randtest(Dgen, Denv, nrepet = 9999)
  
  # subsetting so that the mantel test can be run by city - the matrix values themselves don't change.
  ibe_results <- list()
  for(city_i in 1:length(levels(cities_))) {
    total_envs <- d_env %>% filter(city_abbv == levels(cities_)[city_i]) %>% select(-city_abbv)
    if(nrow(total_envs) > 1){
      # Subset the genind object by city and make dist matrix
      gen_city <- poppr::popsub(gen_, sublist = popNames(gen_)[stringr::str_detect(popNames(gen_), paste0("^",levels(cities_)[city_i]) )])
      # Convert to genpop
      gen_genpop_city <- genind2genpop(gen_city)
      # Generate dissimilarity matrices
      Dgen_city <- dist.genpop(gen_genpop_city, method = 2)
      Denv_city <- dist(total_envs)
      
      ibe_city <- mantel.randtest(Dgen_city, Denv_city, nrepet = 9999)
      ibe_results[[levels(cities_)[city_i]]] <- ibe_city
    } else {
      ibe_results[[levels(cities_)[city_i]]] <- NULL
    }
  }
  
  # Combine all ibe
  ibe_results_all <- c(list(ibe), ibe_results)
  
  return(list(ibe_results_all, Denv, Dgen))
}

make_matrix_df <- function(spp_, Denv, Dgen) {
  Denv_df <- 
    reshape2::melt(as.matrix(Denv)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Denv = value)
  
  Dgen_df <- 
    reshape2::melt(as.matrix(Dgen)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Dgen = value)
  
  df <- full_join(Denv_df, Dgen_df)
  df$spp <- spp_
  
  return(df)
}


extract_ibe_stats_and_plots <- function(env_var_to_use = "nlcd_urban_pct") {
  CD <- do_ibe(spp_ = "CD", env_var_to_use = env_var_to_use)
  CD_ <- CD[[1]][[1]]
  CD_ <- c(CD_$obs, CD_$alter, CD_$rep, CD_$expvar, CD_$pvalue)
  
  DS <- do_ibe(spp_ = "DS", env_var_to_use = env_var_to_use)
  DS_ <- DS[[1]][[1]]
  DS_ <- c(DS_$obs, DS_$alter, DS_$rep, DS_$expvar, DS_$pvalue)
  
  EC <- do_ibe(spp_ = "EC", env_var_to_use = env_var_to_use)
  EC_ <- EC[[1]][[1]]
  EC_ <- c(EC_$obs, EC_$alter, EC_$rep, EC_$expvar, EC_$pvalue)
  
  LS <- do_ibe(spp_ = "LS", env_var_to_use = env_var_to_use)
  LS_ <- LS[[1]][[1]]
  LS_ <- c(LS_$obs, LS_$alter, LS_$rep, LS_$expvar, LS_$pvalue)
  
  PA <- do_ibe(spp_ = "PA", env_var_to_use = env_var_to_use)
  PA_ <- PA[[1]][[1]]
  PA_ <- c(PA_$obs, PA_$alter, PA_$rep, PA_$expvar, PA_$pvalue)
  
  TO <- do_ibe(spp_ = "TO", env_var_to_use = env_var_to_use)
  TO_ <- TO[[1]][[1]]
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
  readr::write_csv(out, paste0("output/IBE/isolation-by-envt-", env_var_to_use,"-mantel-test.csv"))
  
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
  
  mantel_stat_lab <- data.frame(
    spp = factor(c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")),
    r = round(as.numeric(out$Observation), 4)
  )
  
  gg <- ggplot(data = df, aes(Denv, Dgen)) +
    geom_density_2d_filled() +
    geom_text(x = max(df$Denv) * 0.75, y = max(df$Dgen) * 0.2, aes(label = paste0("r = ", r)), data = mantel_stat_lab) +
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
  
  ggsave(paste0("figures/IBE/isolation_by_envt_", env_var_to_use,".png"),
         dpi = "print",
         height = 5,
         width = 8,
         units = "in") 
  
  plot_data <- list(df, mantel_stat_lab)
  saveRDS(plot_data, paste0("output/IBE/isolation_by_envt_", env_var_to_use,"_plot_data.rds"))
}


ibe_mega_plot <- function(){
  d1 <- readRDS("output/IBE/isolation_by_envt_nlcd_urban_pct_plot_data.rds")
  d1[[1]]$env_var <- "% Urban cover"
  d1[[2]]$env_var <- "% Urban cover"
  d1[[2]]$x_coord <- max(d1[[1]]$Denv) * 0.75
  d2 <- readRDS("output/IBE/isolation_by_envt_distance_to_city_center_km_plot_data.rds")
  d2[[1]]$env_var <- "Dist. to city center (km)"
  d2[[2]]$env_var <- "Dist. to city center (km)"
  d2[[2]]$x_coord <- max(d2[[1]]$Denv) * 0.75
  d3 <- readRDS("output/IBE/isolation_by_envt_soiltemp_2.5cm_Apr_12pm_plot_data.rds")
  d3[[1]]$env_var <- "April soil temp. (C)"
  d3[[2]]$env_var <- "April soil temp. (C)"
  d3[[2]]$x_coord <- max(d3[[1]]$Denv) * 0.75
  d4 <- readRDS("output/IBE/isolation_by_envt_soiltemp_2.5cm_Jul_12pm_plot_data.rds")
  d4[[1]]$env_var <- "July soil temp. (C)"
  d4[[2]]$env_var <- "July soil temp. (C)"
  d4[[2]]$x_coord <- max(d4[[1]]$Denv) * 0.75
  
  df <- rbind(d1[[1]], d2[[1]], d3[[1]], d4[[1]])
  df$env_var <- factor(
    df$env_var,
    levels = c(
      "% Urban cover", "Dist. to city center (km)", "April soil temp. (C)", "July soil temp. (C)"
    )
  )
  
  df2 <- rbind(d1[[2]], d2[[2]], d3[[2]], d4[[2]])
  df2$y_coord <- max(d1[[1]]$Dgen) * 0.2
  df2$env_var <- factor(
    df2$env_var,
    levels = c(
      "% Urban cover", "Dist. to city center (km)", "April soil temp. (C)", "July soil temp. (C)"
    )
  )
  
  df_smooth <- df %>% filter(
    (spp == "Bermuda grass" & env_var == "Dist. to city center (km)") |
      env_var == "April soil temp. (C)" |
      env_var == "July soil temp. (C)")
  
  ggplot(data = df, aes(Denv, Dgen)) +
    geom_density_2d_filled() +
    geom_point(shape = 1, alpha = 0.1) +
    geom_text(
      aes(
        x = x_coord,
        y = y_coord,
        label = paste0("r = ", r)), 
      data = df2
    ) +
    facet_grid(spp~env_var, scales = "free_x") +
    scale_fill_manual(values = c("white",  rev(viridisLite::magma(n = 13)))) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_smooth(formula = y~x, se = FALSE, data = df_smooth, color = "red") +
    labs(
      x = paste0("Environmental distance"),
      y = "Genetic distance"
    )
  
  ggsave(paste0("figures/IBE/isolation_by_envt_all.png"),
         dpi = "print",
         height = 10,
         width = 8,
         units = "in") 
}