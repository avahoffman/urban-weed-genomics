# Normals can be accessed here: https://www.ncei.noaa.gov/access/us-climate-normals
library(tidyverse)
library(cowplot)
library(ggplot2)

get_climate_normals <- function(city_ = "BA", temp_ = T) {
  cn <- readr::read_csv(file = "data/normals-monthly-1991-2020-2025-07-31T13.csv")
  cn <- cn %>% mutate(
    city = case_when(
      STATION == "USW00093721" ~ "BA",
      STATION == "USW00014739" ~ "BO",
      STATION == "USW00093134" ~ "LA",
      STATION == "USW00014922" ~ "MN",
      STATION == "USW00023183" ~ "PX"
    )
  ) %>% mutate(
    # Convert to Celsius
    `MLY-TAVG-NORMAL` = (`MLY-TAVG-NORMAL` - 32) * 5 / 9,
    `MLY-TMAX-NORMAL` = (`MLY-TMAX-NORMAL` - 32) * 5 / 9,
    `MLY-TMIN-NORMAL` = (`MLY-TMIN-NORMAL` - 32) * 5 / 9
  ) %>% mutate(# Convert to mm
    `MLY-PRCP-NORMAL` = `MLY-PRCP-NORMAL` * 25.4,
    `MLY-SNOW-NORMAL` = `MLY-SNOW-NORMAL` * 25.4,) %>%
    pivot_longer(
      c(
        `MLY-PRCP-NORMAL`,
        `MLY-SNOW-NORMAL`,
        `MLY-TAVG-NORMAL`,
        `MLY-TMAX-NORMAL`,
        `MLY-TMIN-NORMAL`
      )
    ) %>%
    dplyr::filter(city == city_)
  
  if (temp_) {
    gg <-
      cn %>%
      dplyr::filter(name %in% c("MLY-TAVG-NORMAL", "MLY-TMAX-NORMAL", "MLY-TMIN-NORMAL")) %>%
      ggplot(aes(x = DATE, y = value, color = name)) +
      geom_line() +
      scale_color_manual(values = c("#e36829", "#9d0309", "#fcbf1c")) +
      ylim(c(-15, 43)) +
      theme_classic()
    gg
  } else {
    gg <-
      cn %>%
      dplyr::filter(name %in% c("MLY-PRCP-NORMAL", "MLY-SNOW-NORMAL")) %>%
      ggplot(aes(x = DATE, y = value, color = name)) +
      geom_line() +
      scale_color_manual(values = c("#476e91", "#5cdaf8")) +
      ylim(c(-1, 400)) +
      theme_classic()
    gg
  }
  
}


make_normals_all_cities <- function() {
  g1 <- get_climate_normals(city_ = "BA", temp_ = T)
  g2 <- get_climate_normals(city_ = "BO", temp_ = T)
  g3 <- get_climate_normals(city_ = "LA", temp_ = T)
  g4 <- get_climate_normals(city_ = "MN", temp_ = T)
  g5 <- get_climate_normals(city_ = "PX", temp_ = T)
  
  g6 <- get_climate_normals(city_ = "BA", temp_ = F)
  g7 <- get_climate_normals(city_ = "BO", temp_ = F)
  g8 <- get_climate_normals(city_ = "LA", temp_ = F)
  g9 <- get_climate_normals(city_ = "MN", temp_ = F)
  g10 <- get_climate_normals(city_ = "PX", temp_ = F)
  
  legend1 <- get_legend(g3)
  legend2 <- get_legend(g8)
  
  mega_plot <-
    plot_grid(
      g1 + theme(legend.position = "none"),
      g2 + theme(legend.position = "none"),
      g3 + theme(legend.position = "none"),
      g4 + theme(legend.position = "none"),
      g5 + theme(legend.position = "none"),
      legend1,
      g6 + theme(legend.position = "none"),
      g7 + theme(legend.position = "none"),
      g8 + theme(legend.position = "none"),
      g9 + theme(legend.position = "none"),
      g10 + theme(legend.position = "none"),
      legend2,
      align = 'v',
      axis = "l",
      #hjust = -1,
      nrow = 2
    )
  mega_plot
  
  ggsave(
    paste0("figures/Fig2_sampling_map/climate_norms_temp_ALL.png"),
    dpi = "print",
    height = 2,
    width = 11,
    units = "in"
  )
}
