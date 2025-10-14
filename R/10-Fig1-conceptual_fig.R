library(magrittr)
library(ggplot2)
library(cowplot)

panel1_dat <- function() {
  set.seed(420)
  dat <- data.frame(
    ID = seq(1:30),
    PC1 = sample(1:100, 30, replace = FALSE),
    PC2 = sample(1:100, 30, replace = FALSE),
    city = sample(1:100, 30, replace = FALSE),
    urban = sample(1:100, 30, replace = FALSE)
  )
  dat <-
    dat %>% arrange(city) %>%
    mutate(city = rep(c(
      "City 1", "City 2", "City 3", "City 4", "City 5"
    ), each = 6))  %>%
    mutate(urban = rep(c("More impervious", "Less impervious"), times = 15))
  
  dat <-
    dat %>% arrange(ID)
  
  return(dat)
}

panel2_dat <- function() {
  set.seed(423)
  dat <- data.frame(
    ID = seq(1:30),
    PC1 = c(
      sample(1:30, 6, replace = FALSE),
      sample(10:33, 6, replace = FALSE),
      sample(43:58, 6, replace = FALSE),
      sample(70:100, 6, replace = FALSE),
      sample(55:75, 6, replace = FALSE)
    ),
    PC2 = c(
      sample(1:30, 6, replace = FALSE),
      sample(40:70, 6, replace = FALSE),
      sample(75:100, 6, replace = FALSE),
      sample(53:70, 6, replace = FALSE),
      sample(10:30, 6, replace = FALSE)
    ),
    city = rep(c(
      "City 1", "City 2", "City 3", "City 4", "City 5"
    ), each = 6),
    urban = sample(1:100, 30, replace = FALSE)
  )
  dat <-
    dat %>% arrange(city, urban) %>%
    mutate(urban = rep(c("More impervious", "Less impervious"), times = 15))
  
  dat <-
    dat %>% arrange(ID)
  
  return(dat)
}

panel3_dat <- function() {
  set.seed(423)
  dat <- data.frame(
    ID = seq(1:30),
    PC1 = c(
      sample(1:30, 15, replace = FALSE),
      sample(60:100, 15, replace = FALSE)
    ),
    PC2 = c(
      sample(1:70, 15, replace = FALSE),
      sample(30:100, 15, replace = FALSE)
    ),
    city = sample(1:100, 30, replace = FALSE),
    urban = rep(c("More impervious", "Less impervious"), each = 15)
  )
  dat <-
    dat %>% arrange(urban, city) %>%
    mutate(city = rep(c(
      "City 1", "City 2", "City 3", "City 4", "City 5"
    ), times = 6))
  
  dat <-
    dat %>% arrange(ID)
  
  return(dat)
}


panel4_dat <- function() {
  set.seed(426)
  dat <- data.frame(
    ID = seq(1:30),
    PC1 = c(
      sample(1:35, 3, replace = FALSE),
      sample(25:50, 3, replace = FALSE),
      sample(5:25, 3, replace = FALSE),
      sample(25:50, 3, replace = FALSE),
      sample(65:85, 3, replace = FALSE),
      sample(58:72, 3, replace = FALSE),
      sample(90:100, 3, replace = FALSE),
      sample(60:80, 3, replace = FALSE),
      sample(60:85, 3, replace = FALSE),
      sample(50:70, 3, replace = FALSE)
    ),
    PC2 = c(
      sample(1:30, 3, replace = FALSE),
      sample(30:50, 3, replace = FALSE),
      sample(70:100, 3, replace = FALSE),
      sample(55:75, 3, replace = FALSE),
      sample(1:30, 3, replace = FALSE),
      sample(30:50, 3, replace = FALSE),
      sample(50:70, 3, replace = FALSE),
      sample(50:70, 3, replace = FALSE),
      sample(90:100, 3, replace = FALSE),
      sample(70:90, 3, replace = FALSE)
    ),
    city = rep(c(
      "City 1", "City 2", "City 3", "City 4", "City 5"
    ), each = 6),
    urban = factor(rep(
      c("More impervious", "Less impervious"),
      each = 3,
      times = 5
    ), levels = c("More impervious", "Less impervious"))
  ) %>%
    
    return(dat)
}


make_conceptual_plot <- function(data_) {
  gg <- data_ %>%
    ggplot(aes(
      x = PC1,
      y = PC2,
      color = city,
      shape = urban
    )) +
    geom_point(size = 6, stroke = 1.5) +
    scale_shape_manual(values = c(16, 1)) +
    scale_color_viridis_d() +
    ylim(c(-5, 105)) +
    xlim(c(-5, 105)) +
    theme_classic() +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 12)
    )
  gg
}


#####

make_fig1 <- function() {
  gg1 <- make_conceptual_plot(panel1_dat())
  gg2 <- make_conceptual_plot(panel2_dat())
  gg3 <- make_conceptual_plot(panel3_dat())
  gg4 <- make_conceptual_plot(panel4_dat())
  
  plot_grid(
    plot_grid(
      gg1 + theme(legend.position = "none") + ggtitle("\nNot distinct") + xlab(NULL) + theme(plot.margin = unit(c(-3, 0, 0, 5), "mm")),
      gg2 + theme(legend.position = "none") + ggtitle("\nCities distinct") + ylab(NULL) + xlab(NULL)  + theme(plot.margin = unit(c(-3, 0, 0, 10), "mm")),
      gg3 + theme(legend.position = "none") + ggtitle("Urban environment\ndistinct") + theme(plot.margin = unit(c(2, 0, 0, 5), "mm")),
      gg4 + theme(legend.position = "none") + ggtitle("City and urban\nenvironment distinct") + ylab(NULL) + theme(plot.margin = unit(c(2, 0, 0, 10), "mm")),
      
      #hjust = -1,
      nrow = 2,
      labels = c("(a)", "(b)", "(c)", "(d)"),
      rel_heights = c(95, 100),
      rel_widths = c(93, 100)
    ),
    get_legend(gg1),
    nrow = 1,
    rel_widths = c(3, 1)
  )
  
  ggsave(
    paste0("figures/Fig1_conceptual_fig/Fig1.png"),
    dpi = "print",
    height = 5,
    width = 7,
    units = "in"
  )
}
