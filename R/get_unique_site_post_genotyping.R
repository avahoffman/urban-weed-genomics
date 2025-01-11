get_cities_and_sites <- function(spp_, k){
  file <-
    paste0("output/structure/structure_out_", spp_, k, "_naive_final_f")
  df <- parseStructure(file)
  
  # Import population information / fix sample names
  popmap <-
    read.table(paste0("SNP_data/", spp_, "/popmap_", spp_, "_polyrad.txt"),
               sep = "\t")
  popmap_ordered <- dplyr::arrange(popmap, V1) # Sort by sample name
  df$sample <- popmap_ordered$V1
  df <- df[, -1]
  
  # Join to cover info and other relevant sample details
  samp <- get_sample_info(spp_)
  df <- df %>% left_join(samp)
  
  # Reorder by urban % cover
  df <-
    df %>%
    mutate(sample = fct_reorder(sample, nlcd_urban_pct)) %>%
    arrange(sample)
  
  # Pivot so that Ks (however many there are) and associated assignments are in 
  # long form
  long_df <-  
    pivot_longer(df, cols = !c(sample, pct_miss, city, site_abbv, nlcd_urban_pct, id, management_type))
  
  # ----- Filter out some data and do some cleanup -----
  
  # Remove an outlying sample
  long_df <- long_df %>% filter(sample != "EC.BO.R4.U.1")
  
  # Only 2 PA from Minneapolis. Include when tallying sites.
  long_df <- long_df %>% mutate(spp = spp_)
  
  return(long_df)
}

calc_unique_sites_ <- function(){
  # All unique sites across all species w/in city
  all_dat <-
    rbind(
      get_cities_and_sites("CD", k = 3),
      get_cities_and_sites("DS", k = 3),
      get_cities_and_sites("EC", k = 2),
      get_cities_and_sites("LS", k = 3),
      get_cities_and_sites("PA", k = 4),
      get_cities_and_sites("TO", k = 3)
    )
    
  all_dat %>%
    dplyr::select(city, site_abbv) %>%
    group_by(city) %>%
    distinct() %>%
    count()
  
  ddd <- all_dat %>%
    dplyr::select(spp, city, site_abbv) %>%
    group_by(spp,city) %>%
    distinct() %>%
    count()
}