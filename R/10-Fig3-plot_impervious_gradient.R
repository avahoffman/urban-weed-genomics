# This code creates a gradient of urbanness for Figure 1

library(ggplot2) # ggplot geom_tile aes scale_size_binned_area theme_void xlab theme element_text margin ggsave
library(viridis) # scale_fill_viridis scale_color_viridis

plot_gradient <-
  function() {
    
    urban_gradient <- 
      data.frame(
        x = seq(1:100),
        y = 1
      )

    gg <-
      ggplot() +
      geom_tile(
        data = urban_gradient,
        aes(
          x = x,
          y = y,
          fill = x,
          color = x
        )
      ) +
      scale_fill_viridis(option = "A") +
      scale_color_viridis(option = "A") + #prevents white borders on tile
      scale_size_binned_area(
        limits = c(0, 25),
        breaks = c(1, 5, 15)
      ) +
      theme_void() +
      xlab("% Impervious") +
      theme(
        axis.text.x = element_text(hjust = 0.5, size = 20, vjust = 0.8, margin = margin(10, 10, 10, 10)),
        axis.title.x.bottom = element_text(hjust = 0.5, size = 20, vjust = 0),
        legend.position = "none"
        ) 
    
    gg + theme(plot.margin = unit(c(0, 0, 0.2, 0), "cm"))
    
    ggsave(paste0("figures/Fig1_sampling_map/site_image_gradient.png"),
           height = 1.5,
           width = 20,
           dpi = "print")
    return(gg)
  }
