

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
  colnames(ancestry_df)[1:4] <- c("id", "sample", "pct_miss", "city")
  for (colname_ in colnames(ancestry_df)[5:(5 + k - 1)]) {
    ancestry_df[, colname_] <- as.numeric(ancestry_df[, colname_])
  }
  dat <- ancestry_df
  
  # Reinstitute meaningful city names
  if(any(dat$city == "1")) dat[dat$city == "1",]$city <- "Baltimore"
  if(any(dat$city == "2")) dat[dat$city == "2",]$city <- "Boston"
  if(any(dat$city == "3")) dat[dat$city == "3",]$city <- "Los Angeles"
  if(any(dat$city == "4")) dat[dat$city == "4",]$city <- "Minneapolis"
  if(any(dat$city == "5")) dat[dat$city == "5",]$city <- "Phoenix"
  
  return(dat)
}


get_sample_info <- function(spp_){
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
                                species_name,
                                width = 12,
                                height = 3.5){
  file <- paste0("output/structure/structure_out_", spp_, "_naive_f")
  df <- parseStructure(file)
  
  # Import population information / fix sample names
  popmap <- read.table(paste0("SNP_data/", spp_, "/popmap_", spp_, "_polyrad.txt"), sep = "\t")
  popmap_ordered <- dplyr::arrange(popmap, V1)
  df$sample <- popmap_ordered$V1
  df <- df[,-1]
  
  long_df <- pivot_longer(df, cols = !c(sample, pct_miss, city))
  
  # Get cover info and other relevant sample details
  samp <- get_sample_info(spp_)
  long_df <- long_df %>% left_join(samp)
  
  # Create a new column with cover percent, make a factor, and sort
  long_df <-
    long_df %>%
    mutate(site_n_cov = paste0(round(nlcd_urban_pct, 0), "%")) %>% 
    mutate(site_n_cov = fct_reorder(site_n_cov, nlcd_urban_pct)) %>% 
    arrange(city, nlcd_urban_pct)
  
  # Create a set of labels for the x axis
  x_lbl_ <-
    long_df %>% 
    select(sample, city, site_n_cov, nlcd_urban_pct) %>% 
    unique()
  
  # Keep only the first label in the list to keep things tidy
  x_lbl <- x_lbl_ %>% 
    mutate(dupe = x_lbl_ %>% select(-sample) %>% duplicated()) %>% 
    mutate(site_n_cov = as.character(site_n_cov)) %>% 
    mutate(x = case_when(
      dupe == FALSE ~ site_n_cov,
      TRUE ~ "-"
    )) %>% 
    mutate(x = as_factor(x))
  
  # Create labels and set palette
  colors_ <- viridis::viridis(n = 5, option = "H", begin = 0.2)
  my_pal <- setNames(colors_, 
                     c("K1","K2","K3","K4","K5"))
  
  gg <- ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
    geom_col(width = 1, color = NA) +
    facet_nested(
      ~ city,
      scales = "free_x",
      space = "free",
      switch = "both"
    ) +
    theme_classic() +
    facetted_pos_scales(
      x = list(
        city == "Baltimore" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Baltimore"), ]$x),
        city == "Boston" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Boston"), ]$x),
        city == "Los Angeles" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Los Angeles"), ]$x),
        city == "Minneapolis" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Minneapolis"), ]$x),
        city == "Phoenix" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Phoenix"), ]$x)
      )
    ) +
    scale_fill_manual(values = c(my_pal)) +
    theme(
      legend.position = "none",
      axis.text.x.top = element_text(
        angle = 90,
        hjust = 0,
        size = 5,
        vjust = 0.5
      ),
      #axis.text.x = element_blank(),
      panel.spacing = unit(0, "lines"),
      ggh4x.facet.nestline = element_line(linetype = 3),
      axis.ticks.length = unit(-1, "inch"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x.top = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(paste0("      ", species_name))
  
  setwd(here())
  ggsave(
    paste0("figures/genetics/structure_", species_name, ".png"),
    dpi = "print",
    width = width,
    height = height
  )
  return(gg)
  
}


make_structure_multi_plot <- function(width = 12,
                                      height = 12) {
  p1 <- make_structure_plot(spp = "CD",
                            species_name = "Bermuda grass")
  p2 <- make_structure_plot(spp = "EC",
                            species_name = "crabgrass")
  p3 <- make_structure_plot(spp = "EC",
                            species_name = "horseweed")
  p4 <- make_structure_plot(spp = "LS",
                            species_name = "prickly lettuce")
  p5 <- make_structure_plot(spp = "PA",
                            species_name = "bluegrass")
  p6 <- make_structure_plot(spp = "TO",
                            species_name = "dandelion")
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    align = 'v',
    axis = "l",
    #hjust = -1,
    ncol = 1,
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here())
  ggsave(
    paste0("figures/genetics/structure_ALL.png"),
    dpi = "print",
    width = width,
    height = height
  )
}

