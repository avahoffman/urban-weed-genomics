# Test for isolation by distance
# Mantel test for correlation between two distance matrices
# https://www.uwyo.edu/dbmcd/molmark/lect06/lect6.html has a nice summary of
# distance metrics in genomics.
library(adegenet)
library(tidyverse)


do_ibd <- function(spp_) {
  # Load in genind format data
  gen_ <-
    readRDS(paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno_genind.rds"))
  
  # Extract site names (Will use these as pops)
  cities_ <- pop(gen_)
  pop(gen_) <-
    str_replace(rownames(gen_$tab), "......", "") %>% str_replace("..$", "")
  site_city <-
    data.frame(city_abbv = cities_, site_abbv = pop(gen_))
  
  # Get lat long from the site info
  # Make sure to include management type as otherwise site names aren't unique
  # Get the site info
  site_info <-
    read_csv("data/site_data_urban_cov.csv") %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  nlcd_urban_pct,
                  lat,
                  long) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", site_abbv, management_type, sep = ".")
  
  # Clean and order geographic data
  d_geo <-
    site_city %>%
    left_join(site_info) %>%
    distinct() %>%
    dplyr::select(site_abbv, lat, long) %>%
    column_to_rownames(var = "site_abbv")
  
  # Convert to genpop
  gen_genpop <- genind2genpop(gen_)
  
  # Generate dissimilarity matrices
  Dgen <- dist.genpop(gen_genpop, method = 2)
  Dgeo <- dist(d_geo)
  
  # Run monte-carlo test
  message("Running mantel test...")
  ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 9999)
  return(list(ibd, Dgeo, Dgen))
}

make_matrix_df <- function(spp_, Dgeo, Dgen) {
  Dgeo_df <- reshape2::melt(as.matrix(Dgeo)) %>% rename(Dgeo = value)
  Dgen_df <-
    reshape2::melt(as.matrix(Dgen)) %>% rename(Dgen = value)
  df <- full_join(Dgeo_df, Dgen_df)
  df$spp <- spp_
  
  return(df)
}


extract_ibd_stats_and_plots <- function() {
  CD <- do_ibd(spp_ = "CD")
  CD_ <- CD[[1]]
  CD_ <- c(CD_$obs, CD_$alter, CD_$rep, CD_$expvar, CD_$pvalue)
  
  DS <- do_ibd(spp_ = "DS")
  DS_ <- DS[[1]]
  DS_ <- c(DS_$obs, DS_$alter, DS_$rep, DS_$expvar, DS_$pvalue)
  
  EC <- do_ibd(spp_ = "EC")
  EC_ <- EC[[1]]
  EC_ <- c(EC_$obs, EC_$alter, EC_$rep, EC_$expvar, EC_$pvalue)
  
  LS <- do_ibd(spp_ = "LS")
  LS_ <- LS[[1]]
  LS_ <- c(LS_$obs, LS_$alter, LS_$rep, LS_$expvar, LS_$pvalue)
  
  PA <- do_ibd(spp_ = "PA")
  PA_ <- PA[[1]]
  PA_ <- c(PA_$obs, PA_$alter, PA_$rep, PA_$expvar, PA_$pvalue)
  
  TO <- do_ibd(spp_ = "TO")
  TO_ <- TO[[1]]
  TO_ <- c(TO_$obs, TO_$alter, TO_$rep, TO_$expvar, TO_$pvalue)
  
  out <- data.frame(rbind(CD_, DS_, EC_, LS_, PA_, TO_))
  out$spp <- c("Bermuda grass (CD)", "crabgrass (DS)", "horseweed (EC)", "prickly lettuce (LS)", "bluegrass (PA)", "dandelion (TO)")
  colnames(out) <-
    c(
      "Observation",
      "Hypothesis",
      "Reps",
      "Std.Obs",
      "Expectation",
      "Variance",
      "p-value",
      "Species"
    )
  out <- relocate(out, "Species", .before = "Observation")
  readr::write_csv(out, "output/isolation-by-distance-mantel-test.csv")
  
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
  
  ggplot(data = df, aes(Dgeo, Dgen)) +
    geom_density_2d_filled() +
    facet_wrap(~spp) +
    geom_point(shape = 1) +
    scale_fill_manual(values = c("white",  rev(viridisLite::magma(n = 13)))) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_smooth(formula = y~x, se = FALSE, linetype = 2) +
    labs(
      x = "Geographic distance",
      y = "Genetic distance"
    )
 
  ggsave(paste0("figures/isolation_by_distance.png"),
         dpi = "print",
         height = 5,
         width = 8,
         units = "in") 
}

