#' Test for isolation by distance
#' Mantel test for correlation between two distance matrices
#' https://www.uwyo.edu/dbmcd/molmark/lect06/lect6.html has a nice summary of
#' distance metrics in genomics.
#'
#' Also found the adegenet docs helpful! 
#' https://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf

library(adegenet)
library(tidyverse)


do_ibd <- function(spp_) {
  # Load in genind format data
  gen_ <-
    readRDS(paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno_genind.rds"))
  
  # Extract site names (Will use these as pops)
  cities_ <- pop(gen_)
  pop(gen_) <-
    paste0(pop(gen_), ".", str_replace(rownames(gen_$tab), "......", "") %>% str_replace("..$", ""))
  site_city <-
    data.frame(city_abbv = cities_, site_abbv = pop(gen_))
  
  # Get lat long from the site info
  # Make sure to include management type as otherwise site names aren't unique
  # Get the site info
  site_info <-
    read_csv("data/site_data_urban_cov.csv") %>%
    dplyr::select(city,
                  site_abbv,
                  city_abbv,
                  management_type,
                  nlcd_urban_pct,
                  lat,
                  long) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".")
  
  # Clean and order geographic data
  d_geo <-
    site_city %>%
    left_join(site_info) %>%
    distinct() %>%
    dplyr::select(city_abbv, site_abbv, lat, long) %>%
    column_to_rownames(var = "site_abbv")
  
  # Convert to genpop
  gen_genpop <- genind2genpop(gen_)
  
  # Generate dissimilarity matrices
  Dgen <- dist.genpop(gen_genpop, method = 2)
  Dgeo <- dist(d_geo %>% select(-city_abbv))
  
  # Run monte-carlo mantel test
  message("Running mantel test...")
  set.seed(444)
  ibd <- mantel.randtest(Dgen, Dgeo, nrepet = 9999)
  
  # subsetting so that the mantel test can be run by city - the matrix values themselves don't change.
  ibd_results <- list()
  for(city_i in 1:length(levels(cities_))) {
    total_geos <- d_geo %>% filter(city_abbv == levels(cities_)[city_i]) %>% select(-city_abbv)
    if(nrow(total_geos) > 1){
      # Subset the genind object by city and make dist matrix
      gen_city <- poppr::popsub(gen_, sublist = popNames(gen_)[stringr::str_detect(popNames(gen_), paste0("^",levels(cities_)[city_i]) )])
      # Convert to genpop
      gen_genpop_city <- genind2genpop(gen_city)
      # Generate dissimilarity matrices
      Dgen_city <- dist.genpop(gen_genpop_city, method = 2)
      Dgeo_city <- dist(total_geos)
      
      set.seed(444)
      ibd_city <- mantel.randtest(Dgen_city, Dgeo_city, nrepet = 9999)
      ibd_results[[levels(cities_)[city_i]]] <- ibd_city
    } else {
      ibd_results[[levels(cities_)[city_i]]] <- NULL
    }
  }
  
  # Combine all ibd
  ibd_results_all <- c(list(ibd), ibd_results)
  
  return(list(ibd_results_all, Dgeo, Dgen))
}

make_matrix_df <- function(spp_, Dgeo, Dgen) {
  Dgeo_df <- 
    reshape2::melt(as.matrix(Dgeo)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Dgeo = value)
  
  Dgen_df <- 
    reshape2::melt(as.matrix(Dgen)) %>% 
    filter(as.numeric(Var1) > as.numeric(Var2)) %>% 
    rename(Dgen = value)
  
  df <- full_join(Dgeo_df, Dgen_df)
  df$spp <- spp_
  
  return(df)
}


extract_ibd_stats_and_plots <- function() {
  
  # DO ALL CITIES FIRST
  
  CD <- do_ibd(spp_ = "CD")
  CD_ <- CD[[1]][[1]]
  CD_ <- c(CD_$obs, CD_$alter, CD_$rep, CD_$expvar, CD_$pvalue)
  
  DS <- do_ibd(spp_ = "DS")
  DS_ <- DS[[1]][[1]]
  DS_ <- c(DS_$obs, DS_$alter, DS_$rep, DS_$expvar, DS_$pvalue)
  
  EC <- do_ibd(spp_ = "EC")
  EC_ <- EC[[1]][[1]]
  EC_ <- c(EC_$obs, EC_$alter, EC_$rep, EC_$expvar, EC_$pvalue)
  
  LS <- do_ibd(spp_ = "LS")
  LS_ <- LS[[1]][[1]]
  LS_ <- c(LS_$obs, LS_$alter, LS_$rep, LS_$expvar, LS_$pvalue)
  
  PA <- do_ibd(spp_ = "PA")
  PA_ <- PA[[1]][[1]]
  PA_ <- c(PA_$obs, PA_$alter, PA_$rep, PA_$expvar, PA_$pvalue)
  
  TO <- do_ibd(spp_ = "TO")
  TO_ <- TO[[1]][[1]]
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
  readr::write_csv(out, "output/IBD/isolation-by-distance-mantel-test.csv")
  
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
    read_csv("data/site_data_urban_cov.csv") %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  nlcd_urban_pct,
                  lat,
                  long) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite("site_abbv", city_abbv, site_abbv, management_type, sep = ".", remove = F)
  
  mantel_stat_lab <- data.frame(
    spp = factor(c("Bermuda grass", "crabgrass", "horseweed", "prickly lettuce", "bluegrass", "dandelion")),
    r = round(c(0.447643046, 0.329999156, 0.402833889, 0.193960714, 0.282156217, 0.310145748), 3)
  )

  gg <- ggplot(data = df, aes(Dgeo, Dgen)) +
    geom_density_2d_filled() +
    geom_text(x = max(df$Dgeo) * 0.75, y = max(df$Dgen) * 0.2, aes(label = paste0("r = ", r)), data = mantel_stat_lab) +
    facet_wrap(~spp) +
    geom_point(shape = 1, alpha = 0.1) +
    scale_fill_manual(values = c("white",  rev(viridisLite::magma(n = 13)))) +
    theme_classic() +
    theme(legend.position = "none") +
    geom_smooth(formula = y~x, se = FALSE, color = "red") +
    labs(
      x = "Geographic distance",
      y = "Genetic distance"
    )

  gg
   
  ggsave(paste0("figures/IBD/isolation_by_distance.png"),
         dpi = "print",
         height = 5,
         width = 8,
         units = "in") 
  
  # Now do within city
  # TODO DRY this up!
  CD_BA <- c(CD[[1]]$BA$obs, CD[[1]]$BA$alter, CD[[1]]$BA$rep, CD[[1]]$BA$expvar, CD[[1]]$BA$pvalue)
  CD_LA <- c(CD[[1]]$LA$obs, CD[[1]]$LA$alter, CD[[1]]$LA$rep, CD[[1]]$LA$expvar, CD[[1]]$LA$pvalue)
  CD_PX <- c(CD[[1]]$PX$obs, CD[[1]]$PX$alter, CD[[1]]$PX$rep, CD[[1]]$PX$expvar, CD[[1]]$PX$pvalue)
  
  DS_BA <- c(DS[[1]]$BA$obs, DS[[1]]$BA$alter, DS[[1]]$BA$rep, DS[[1]]$BA$expvar, DS[[1]]$BA$pvalue)
  DS_BO <- c(DS[[1]]$BO$obs, DS[[1]]$BO$alter, DS[[1]]$BO$rep, DS[[1]]$BO$expvar, DS[[1]]$BO$pvalue)
  DS_MN <- c(DS[[1]]$MN$obs, DS[[1]]$MN$alter, DS[[1]]$MN$rep, DS[[1]]$MN$expvar, DS[[1]]$MN$pvalue)
  DS_PX <- c(DS[[1]]$PX$obs, DS[[1]]$PX$alter, DS[[1]]$PX$rep, DS[[1]]$PX$expvar, DS[[1]]$PX$pvalue)
  
  EC_BA <- c(EC[[1]]$BA$obs, EC[[1]]$BA$alter, EC[[1]]$BA$rep, EC[[1]]$BA$expvar, EC[[1]]$BA$pvalue)
  EC_LA <- c(EC[[1]]$LA$obs, EC[[1]]$LA$alter, EC[[1]]$LA$rep, EC[[1]]$LA$expvar, EC[[1]]$LA$pvalue)
  EC_PX <- c(EC[[1]]$PX$obs, EC[[1]]$PX$alter, EC[[1]]$PX$rep, EC[[1]]$PX$expvar, EC[[1]]$PX$pvalue)
  
  LS_BA <- c(LS[[1]]$BA$obs, LS[[1]]$BA$alter, LS[[1]]$BA$rep, LS[[1]]$BA$expvar, LS[[1]]$BA$pvalue)
  LS_BO <- c(LS[[1]]$BO$obs, LS[[1]]$BO$alter, LS[[1]]$BO$rep, LS[[1]]$BO$expvar, LS[[1]]$BO$pvalue)
  LS_LA <- c(LS[[1]]$LA$obs, LS[[1]]$LA$alter, LS[[1]]$LA$rep, LS[[1]]$LA$expvar, LS[[1]]$LA$pvalue)
  LS_MN <- c(LS[[1]]$MN$obs, LS[[1]]$MN$alter, LS[[1]]$MN$rep, LS[[1]]$MN$expvar, LS[[1]]$MN$pvalue)
  LS_PX <- c(LS[[1]]$PX$obs, LS[[1]]$PX$alter, LS[[1]]$PX$rep, LS[[1]]$PX$expvar, LS[[1]]$PX$pvalue)
  
  PA_BA <- c(PA[[1]]$BA$obs, PA[[1]]$BA$alter, PA[[1]]$BA$rep, PA[[1]]$BA$expvar, PA[[1]]$BA$pvalue)
  PA_BO <- c(PA[[1]]$BO$obs, PA[[1]]$BO$alter, PA[[1]]$BO$rep, PA[[1]]$BO$expvar, PA[[1]]$BO$pvalue)
  PA_LA <- c(PA[[1]]$LA$obs, PA[[1]]$LA$alter, PA[[1]]$LA$rep, PA[[1]]$LA$expvar, PA[[1]]$LA$pvalue)
  PA_PX <- c(PA[[1]]$PX$obs, PA[[1]]$PX$alter, PA[[1]]$PX$rep, PA[[1]]$PX$expvar, PA[[1]]$PX$pvalue)
  
  TO_BA <- c(TO[[1]]$BA$obs, TO[[1]]$BA$alter, TO[[1]]$BA$rep, TO[[1]]$BA$expvar, TO[[1]]$BA$pvalue)
  TO_BO <- c(TO[[1]]$BO$obs, TO[[1]]$BO$alter, TO[[1]]$BO$rep, TO[[1]]$BO$expvar, TO[[1]]$BO$pvalue)
  TO_LA <- c(TO[[1]]$LA$obs, TO[[1]]$LA$alter, TO[[1]]$LA$rep, TO[[1]]$LA$expvar, TO[[1]]$LA$pvalue)
  TO_MN <- c(TO[[1]]$MN$obs, TO[[1]]$MN$alter, TO[[1]]$MN$rep, TO[[1]]$MN$expvar, TO[[1]]$MN$pvalue)
  TO_PX <- c(TO[[1]]$PX$obs, TO[[1]]$PX$alter, TO[[1]]$PX$rep, TO[[1]]$PX$expvar, TO[[1]]$PX$pvalue)
  
  out <- data.frame(rbind(CD_BA, CD_LA, CD_PX, DS_BA, DS_BO, DS_MN, DS_PX, EC_BA, EC_LA, EC_PX, 
                          LS_BA, LS_BO, LS_LA, LS_MN, LS_PX, PA_BA, PA_BO, PA_LA, PA_PX, TO_BA, TO_BO, TO_LA, TO_MN, TO_PX))
  
  # Do p-value adjustment for however many cities, not across species since they are independent.
  
  out$padj <- c(
    p.adjust(out$V7[1:3], method = "BH"), # CD
    p.adjust(out$V7[4:7], method = "BH"), # DS
    p.adjust(out$V7[8:10], method = "BH"), # EC
    p.adjust(out$V7[11:15], method = "BH"), # LS
    p.adjust(out$V7[16:19], method = "BH"), # PA
    p.adjust(out$V7[20:24], method = "BH") # TO
  )
  out$spp <- c(
    rep("Bermuda grass (CD)", 3), 
    rep("crabgrass (DS)", 4), 
    rep("horseweed (EC)", 3), 
    rep("prickly lettuce (LS)", 5), 
    rep("bluegrass (PA)", 4), 
    rep("dandelion (TO)", 5)
    )
  out$city <- c(
    "BA", "LA", "PX", "BA", "BO", "MN", "PX", "BA", "LA", "PX", "BA", "BO", "LA", "MN", "PX", "BA", "BO", "LA", "PX", "BA", "BO", "LA", "MN", "PX"
  )
  colnames(out) <-
    c(
      "Observation",
      "Hypothesis",
      "Reps",
      "Std.Obs",
      "Expectation",
      "Variance",
      "p-value",
      "Adjusted p-value",
      "Species",
      "City"
    )
  out <- relocate(out, "Species", .before = "Observation")
  readr::write_csv(out, "output/IBD/isolation-by-distance-within-city-mantel-test.csv")
  
}

