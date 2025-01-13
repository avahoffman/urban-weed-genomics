library(tidyverse)
library(here)
library(ggh4x)
library(cowplot)

parseStructure <- function(file) {
  # initial parse
  input_1 <- suppressWarnings(readLines(file))
  
  # Get ancestry part
  input <-
    input_1[(grep("Inferred ancestry of individuals", input_1) + 1):length(input_1)]
  input <-
    input[2:(grep("Estimated Allele Frequencies", input) - 3)]
  
  # Get k (pops)
  k <-
    as.numeric(gsub(pattern = ".*([0-9]+) .*",
                    replacement = "\\1",
                    input_1[grep("populations assumed", input_1)]))
  
  ancestry_df <- data.frame()
  for (i in 1:length(input)) {
    split_line <- str_split(input[i], pattern = "[ ]+")[[1]]
    split_line <- split_line[split_line != "" & split_line != ":"]
    ancestry_df <- rbind(ancestry_df, split_line)
  }
  
  colnames(ancestry_df)[5:(5 + k - 1)] <- paste0("K", seq(1:k))
  colnames(ancestry_df)[1:4] <-
    c("id", "sample", "pct_miss", "city")
  for (colname_ in colnames(ancestry_df)[5:(5 + k - 1)]) {
    ancestry_df[, colname_] <- as.numeric(ancestry_df[, colname_])
  }
  dat <- ancestry_df
  
  # Reinstitute meaningful city names
  if (any(dat$city == "1"))
    dat[dat$city == "1", ]$city <- "Baltimore"
  if (any(dat$city == "2"))
    dat[dat$city == "2", ]$city <- "Boston"
  if (any(dat$city == "3"))
    dat[dat$city == "3", ]$city <- "Los Angeles"
  if (any(dat$city == "4"))
    dat[dat$city == "4", ]$city <- "Minneapolis"
  if (any(dat$city == "5"))
    dat[dat$city == "5", ]$city <- "Phoenix"
  
  return(dat)
}


get_sample_info <- function(spp_) {
  # Get the site info
  site_info <-
    read_csv(paste0(here(),
                    "/data/site_data_urban_cov.csv")) %>%
    dplyr::select(site_abbv, city, management_type, nlcd_urban_pct) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))
  
  # Get samples and spread out the name to get more details
  # Join to site info to get urban cover
  samp <-
    read_tsv(paste0("08_polyRAD/popmap_", spp_, "_polyrad.txt"),
             col_names = c("sample", "city")) %>%
    tidyr::separate(
      sample,
      into = c("spp", "ct", "site_abbv", "management_type", "id"),
      sep = "\\.",
      remove = F
    ) %>%
    left_join(site_info,
              by = c("city", "site_abbv", "management_type")) %>%
    dplyr::select(sample, site_abbv, city, nlcd_urban_pct, id, management_type)
}


make_structure_plot <- function(spp_,
                                k = "",
                                species_name,
                                structure_plot = T,
                                width = 12,
                                height = 3.5) {
  # Read in structure results
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
  
  # Only 2 PA from Minneapolis. Looks bad on the plot. Replace with asterisk!
  long_df <- long_df %>%
    mutate(city = case_when(
      sample %in% c("PA.MN.L01-TO_PA.M.4", "PA.MN.L01-TO_PA.M.5") ~ "*",
      TRUE ~ city
    ))  %>%
    mutate(city = case_when(
      city == "Minneapolis" ~ "Minneapolis-Saint Paul",
      TRUE ~ city
    ))
  
  # ----- Reorder cities -----
  
  if (spp_ == "DS")
    long_df$city <- factor(long_df$city, levels = c("Minneapolis-Saint Paul", "Boston", "Baltimore", "Phoenix"))
  if (spp_ %in% c("LS", "TO"))
    long_df$city <- factor(long_df$city, levels = c("Minneapolis-Saint Paul", "Boston", "Baltimore", "Los Angeles", "Phoenix"))
  if (spp_ == "PA")
    long_df$city <- factor(long_df$city, levels = c("*", "Boston", "Baltimore", "Los Angeles", "Phoenix"))
  
  # Make plot
  if (structure_plot) {
    gg <-
      ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
      geom_col(width = 1, color = NA) +
      facet_nested(~ city,
                   scales = "free_x",
                   space = "free") +
      theme_classic() +
      labs(y = species_name) +
      scale_fill_manual(values = c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF")) + # viridis::turbo(n = 5)
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        ggh4x.facet.nestline = element_line(linetype = 3),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(hjust = 0.1),
        plot.margin = unit(c(0,5,0,5), "pt"),
        text = element_text(size = 8)
      )
    
  } else {
    gg <-
      ggplot(data = long_df, aes(x = sample, y = value, fill = nlcd_urban_pct)) +
      geom_col(width = 1, color = NA) +
      facet_nested(~ city,
                   scales = "free_x",
                   space = "free",
                   switch = "both") +
      theme_void() +
      scale_fill_viridis_c(option = "A") +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        #ggh4x.facet.nestline = element_line(linetype = 3),
        strip.text.x = element_blank(),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(15,5,0,5), "pt")
      )
  }
  
  setwd(here())
  ggsave(
    paste0("figures/genetics/structure_", species_name, ".png"),
    dpi = "print",
    width = width,
    height = height
  )
  return(gg)
  
}

make_structure_multi_plot <- function() {
  # Optimal K=3
  p1 <- plot_grid(
    make_structure_plot(spp = "CD", k = 3, structure_plot = F, species_name = "Bermuda grass"),
    make_structure_plot(spp = "CD", k = 3, species_name = "Bermuda grass"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Optimal K=3
  p2 <- plot_grid(
    make_structure_plot(spp = "DS", k = 3, structure_plot = F, species_name = "crabgrass"),
    make_structure_plot(spp = "DS", k = 3, species_name = "crabgrass"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Optimal K=2
  p3 <- plot_grid(
    make_structure_plot(spp = "EC", k = 2, structure_plot = F, species_name = "horseweed"),
    make_structure_plot(spp = "EC", k = 2, species_name = "horseweed"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Optimal K=3
  p4 <- plot_grid(
    make_structure_plot(spp = "LS", k = 3, structure_plot = F, species_name = "prickly lettuce"),
    make_structure_plot(spp = "LS", k = 3, species_name = "prickly lettuce"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Optimal K=4
  p5 <- plot_grid(
    make_structure_plot(spp = "PA", k = 4, structure_plot = F, species_name = "bluegrass"),
    make_structure_plot(spp = "PA", k = 4, species_name = "bluegrass"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Optimal K=3
  p6 <- plot_grid(
    make_structure_plot(spp = "TO", k = 3, structure_plot = F, species_name = "dandelion"),
    make_structure_plot(spp = "TO", k = 3, species_name = "dandelion"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  
  # Add K plots
  # p1k <- Pr_plot(spp_ = "CD") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  # p2k <- Pr_plot(spp_ = "DS") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  # p3k <- Pr_plot(spp_ = "EC") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  # p4k <- Pr_plot(spp_ = "LS") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  # p5k <- Pr_plot(spp_ = "PA") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  # p6k <- Pr_plot(spp_ = "TO") + theme(plot.margin = unit(c(20,5,0,5), "pt"))
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    #align = 'h',
    #axis = "b",
    hjust = 0,
    ncol = 1
    #rel_widths = c(1, 0.2),
    # labels = c(
    #   "(a)",
    #   "(b)",
    #   "(c)",
    #   "(d)",
    #   "(e)",
    #   "(f)"
    # ),
    # label_size = 10,
    # label_y = 0.8
  )
  mega_plot
  setwd(here())
  
  # Prepare figures such that, after reduction to fit across one column, two-thirds page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/genetics/structure_ALL.png"),
    dpi = "print",
    width = 169,
    height = 190,
    units = "mm"
  )
}
