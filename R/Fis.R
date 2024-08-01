# Note that Fis is calculated in genodive. This is only for modifying the
# structure output appropriately to do so at the site (not city) level.

library(stringr) # CRAN v1.5.1
library(naniar)  # CRAN v1.1.0
library(dplyr)   # CRAN v1.1.4

source("R/plot_structure.R") # get_sample_info

#' Combine Fis by site and urban percent cover datasets.
#' This function uses the output of genodive for Fis.
#'
#' @param spp_ two letter string representing the species (e.g., "CD")
#'
#' @return tibble
#'
#' @examples gather_bysite_fis("CD")
gather_bysite_fis <- function(spp_) {
  # Get by site Fis stats from genodive output
  dat_ <- readLines(paste0(
    "output/population_stats/genodive_output_Fis_",
    spp_,
    "_bysite.gdv"
  ))
  
  # Get the Lines that we should keep
  file_indexes <- NULL
  for (i in 1:length(dat_)) {
    file_indexes <- c(file_indexes, nchar(dat_[i]) > 200)
  }
  dat_ <- dat_[file_indexes] # drop any lines that are short (aka don't have Fis info)
  
  # Split column and clean up the genodive results
  the_column_names <- str_split(dat_[1], "\t")[[1]] # gather column names by splitting apart by tab
  dat_ <- tibble(dat_) %>%
    separate(1, the_column_names, "\t") %>% # split apart by tab
    select(Population, `Multi-locus`) %>%
    filter(Population != "Population" &
             Population != "Overall") %>% # Remove extra headings
    mutate(`Multi-locus` = as.numeric(na_if(`Multi-locus`, "---")))
  dat_ <- dat_ %>%
    mutate(measure = c(rep("Fis", nrow(dat_) / 2), rep("p-value", nrow(dat_) / 2))) %>% # Make sure we distinguish between Fis and pvalues
    filter(measure != "p-value")
  
  # Get the structure sample names (to ensure correct order)
  struct_dat <- read_table(paste0("output/AMOVA/", spp_, "_estimatedgeno_bysite.structure"))
  colnames(struct_dat) <- c("sample", colnames(struct_dat)) # Fix lack of 1st column name
  # Clean pop names so they match the genodive output
  struct_dat <-
    struct_dat %>%
    rename(Population = 2) %>%
    mutate(Population = case_when(
      nchar(Population) == 1 ~ paste0("Pop00", Population),
      nchar(Population) == 2 ~ paste0("Pop0", Population)
    )) %>% select(sample, Population)
  
  # Get urban pct
  samp <- get_sample_info(spp_)
  
  # Merge everything
  dat_merged <-
    struct_dat %>%
    left_join(samp) %>%
    select(!c(sample, id)) %>%
    distinct() %>%
    right_join(dat_)
  
  return(dat_merged)
}


#' Make a scatterplot of % Urban Cover versus Fis Statistic
#'
#' @param spp_ two letter string representing the species (e.g., "CD")
#' @param fis_dat output of gather_bysite_fis()
#' @param species_name full text of species name for the plot title
#'
#' @return ggplot object
#'
#' @examples make_fis_plot("CD", CD_bysite_Fis, "Bermuda grass")
make_fis_plot <- function(spp_, fis_dat, species_name) {
  colors_ <- viridis::viridis(n = 5,
                              option = "H",
                              begin = 0.2)
  my_pal <- setNames(colors_,
                     c(
                       "Baltimore",
                       "Boston",
                       "Los Angeles",
                       "Minneapolis",
                       "Phoenix"
                     ))
  
  shapes_ <- c(21, 24, 22, 23, 25)
  
  shape_pal <- setNames(shapes_,
                        c(
                          "Baltimore",
                          "Boston",
                          "Los Angeles",
                          "Minneapolis",
                          "Phoenix"
                        ))
  gg <-
    ggplot() +
    geom_point(
      data = fis_dat,
      aes(
        x = nlcd_urban_pct,
        y = `Multi-locus`,
        fill = city,
        shape = city
      ),
      size = 2
    ) +
    labs(x = "% Urban Cover", y = "Multi-locus Fis Statistic") +
    scale_fill_manual(values = my_pal) +
    scale_shape_manual(values = shape_pal) +
    guides(fill = guide_legend(title = ""), shape = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = "none", legend.text = element_text(size = 10)) +
    ggtitle(paste0("", species_name))
  
  ggsave(gg,
         filename = paste0("figures/Fis/Fis_bysite_byurban_", spp_, ".png"))
  
  return(gg)
}


#' Wrapper to make by-site Fis vs. Urban Cover plots
#'
#' @return saves file & prints plot of all species combined.
#'
#' @examples plot_all_bysite_Fis()
plot_all_bysite_Fis <- function() {
  CD_bysite_Fis <- gather_bysite_fis("CD")
  DS_bysite_Fis <- gather_bysite_fis("DS")
  EC_bysite_Fis <- gather_bysite_fis("EC")
  LS_bysite_Fis <- gather_bysite_fis("LS")
  PA_bysite_Fis <- gather_bysite_fis("PA")
  TO_bysite_Fis <- gather_bysite_fis("TO")
  
  p1 <- make_fis_plot("CD", CD_bysite_Fis, "Bermuda grass")
  p2 <- make_fis_plot("DS", DS_bysite_Fis, "crabgrass")
  p3 <- make_fis_plot("EC", EC_bysite_Fis, "horseweed")
  p4 <- make_fis_plot("LS", LS_bysite_Fis, "prickly lettuce")
  p5 <- make_fis_plot("PA", PA_bysite_Fis, "bluegrass")
  p6 <- make_fis_plot("TO", TO_bysite_Fis, "dandelion")
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    align = 'vh',
    ncol = 3,
    rel_heights = c(1, 1, 1),
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here::here())
  ggsave(
    paste0("figures/Fis/Fis_bysite_all.png"),
    dpi = "print",
    width = 7,
    height = 5
  )
}
