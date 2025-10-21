# Prepare data for AMOVA using GenoDive and plot results
library(ggplot2)
library(tidyverse)

clean_in_prep_for_AMOVA <- function(spp_) {
  dat_ <- readr::read_delim(paste0("output/AMOVA/", spp_, "_estimatedgeno.structure"),
                            delim = "\t")
  dat_ <-
    dat_ %>%
    tidyr::separate(
      col = 1,
      into = c("spp", "city", "site", "mgmt", "rep"),
      "\\.",
      remove = F
    ) %>%
    tidyr::unite("site", c(city, site, mgmt), remove = F) %>%
    dplyr::mutate("site" = as.numeric(factor(site)))
  
  # Optional, Look at combos to gut check
  dat_ %>% select(city, site) %>% distinct() %>% arrange(city, site) %>% print(n = 100)
  
  dat_ <- dat_ %>%
    select(-c(spp, city, mgmt, rep, `...2`))
  
  write.table(
    dat_,
    sep = "\t",
    row.names = F,
    quote = F,
    file = paste0("output/AMOVA/", spp_, "_estimatedgeno_bysite.structure")
  )
}


spp_labels <- function() {
  return(
    c(
      CD =  "_C. dactylon_<br>(Bermuda grass)",
      DS = "_D. sanguinalis_<br>(crabgrass)",
      EC = "_E. canadensis_<br>(horseweed)",
      LS = "_L. serriola_<br>(prickly lettuce)",
      PA = "_P. annua_<br>(bluegrass)",
      TO = "_T. officinale_<br>(dandelion)"
    )
  )
}


# Plot AMOVA
make_amova_plot <- function() {
  amova <- read_csv("output/AMOVA/AMOVA_results.csv")
  
  amova <-
    amova %>%
    dplyr::filter(Species != "PA - no MN") %>%
    mutate(Species = case_when(
      Species == "PA - old" ~ "PA", 
      Species == "PA - w/ MN" ~ "PA",
      TRUE ~ Species
    )) %>%
    mutate(`%Var` = `%Var` * 100) %>%
    mutate(
      `Source of Variation` = case_when(
        `Source of Variation` == "Among Population" ~ "Within-city\nvariation",
        `Source of Variation` == "Within Population" ~ "Within-site\nvariation",
        `Source of Variation` == "Among City" ~ "Among-city\nvariation",
        TRUE ~ `Source of Variation`
      )
    ) %>%
    dplyr::select(Species, `Source of Variation`, `%Var`) %>% 
    mutate(Species = spp_labels()[as.character(Species)])

  gg <- 
    ggplot() +
    geom_col(data = amova, aes(x = `%Var`, y = Species, fill = `Source of Variation`)) +
    labs(x = "% of variation", y = NULL) +
    scale_fill_manual(
      values = c("#30123BFF", "#28BBECFF", "#A2FC3CFF"),
      breaks = c("Within-site\nvariation", "Within-city\nvariation", "Among-city\nvariation")
    ) +
    scale_y_discrete(limits = rev(enframe(spp_labels())$value)) + # Reverse due to the way ggplot be stackin'
    guides(fill = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.y = ggtext::element_markdown())
  
  gg
  # Prepare figures such that, after reduction to fit across one column, two-thirds 
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, 
  # all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/Fig6_amova/AMOVA.png"),
    dpi = "print",
    width = 120,
    height = 120,
    units = "mm"
  )
}