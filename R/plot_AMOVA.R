# Plot AMOVA

make_amova_plot <- function() {
  amova <- read_csv("output/AMOVA/AMOVA_results.csv")
  
  amova <-
    amova %>%
    filter(Species != "PA - no MN") %>%
    mutate(Species = case_when(Species == "PA - old" ~ "PA", TRUE ~ Species)) %>%
    mutate(
      Species = case_when(
        Species == "CD" ~ "Bermuda grass",
        Species == "DS" ~ "crabgrass",
        Species == "EC" ~ "horseweed",
        Species == "LS" ~ "prickly lettuce",
        Species == "PA" ~ "bluegrass",
        Species == "TO" ~ "dandelion"
      )
    ) %>%
    mutate(`%Var` = `%Var` * 100) %>%
    mutate(
      `Source of Variation` = case_when(
        `Source of Variation` == "Among Population" ~ "Among Site\n(Within City)",
        `Source of Variation` == "Within Population" ~ "Within Site",
        TRUE ~ `Source of Variation`
      )
    ) %>%
    dplyr::select(Species, `Source of Variation`, `%Var`)
  
  # Reverse due to the way ggplot be stackin'
  amova$Species <-
    factor(amova$Species, levels = rev(
      c(
        "Bermuda grass",
        "crabgrass",
        "horseweed",
        "prickly lettuce",
        "bluegrass",
        "dandelion"
      )
    ))
  
  gg <- 
    ggplot() +
    geom_col(data = amova, aes(x = `%Var`, y = Species, fill = `Source of Variation`)) +
    labs(x = "% of variation", y = NULL) +
    scale_fill_manual(
      values = c("#30123BFF", "#28BBECFF", "#A2FC3CFF"),
      breaks = c("Within Site", "Among Site\n(Within City)", "Among City")
    ) +
    guides(fill = guide_legend(title = "")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  gg
  # Prepare figures such that, after reduction to fit across one column, two-thirds 
  # page width, or two columns (80 mm, 112 mm, or 169 mm, respectively) as required, 
  # all lettering and symbols will be clear and easy to read,
  ggsave(
    paste0("figures/AMOVA.png"),
    dpi = "print",
    width = 120,
    height = 120,
    units = "mm"
  )
}