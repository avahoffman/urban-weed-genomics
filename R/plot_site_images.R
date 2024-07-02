# This code creates a gradient of urbanness with sample site images

library(ggplot2)
library(viridis)
library(cowplot) # theme_map()
library(magick)

plot_sites <-
  function() {
    # This function takes an urban pct cover raster layer and spatial data points and
    # makes a nice ggplot
    #
    # Args:
    # urban_data: output from trim_spatial() (an .rds file)
    #
    # Returns: plot to file
    #
    # Usage:
    # ggplot_urban_pct_cover_plot(
    #     readr::read_rds("spatial_data/trimmed_spatial_BA.rds")
    # )
    
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
      #ylim(0,10) +
      theme_void() +
      xlab("% Urban") +
      theme(
        axis.text.x = element_text(hjust = 0.5, size = 20, vjust = 0.8),
        axis.title.x.bottom = element_text(hjust = 0.5, size = 20, vjust = 0),
        legend.position = "none"
        ) 
    
    gg
    
    # gg +
    #   draw_image(
    #     "images/IMG_5536_BA.LA.U, 70 pct urban_.jpg",
    #     x = 60,
    #     y = 3,
    #     width = 20,
    #     height = 2
    #   ) +
    #   annotate("segment", 
    #            x = 70,
    #            xend = 69, 
    #            y = 1,
    #            yend = 3) +
    #   draw_image(
    #     "images/IMG_20180511_105449984_HDR_PX.ALA.M, 45 pct urban_.jpg",
    #     x = 40,
    #     y = 3,
    #     width = 20,
    #     height = 2
    #   ) +
    #   annotate("segment", 
    #            x = 45,
    #            xend = 43, 
    #            y = 1,
    #            yend = 3)
    
    ggsave(paste0("figures/site_image_gradient.png"),
           height = 1,
           width = 20,
           dpi = "print")
    return(gg)
  }
