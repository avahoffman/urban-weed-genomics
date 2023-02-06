# This code creates a map of samples overlaid on NLCD urban cover data

library(raster)
library(ggplot2)
library(viridis)
library(cowplot) # theme_map()

plot_urban_cover_and_sites <-
  function(spatial_data) {
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
    
    # De-list
    sites <- spatial_data[[1]]
    urban_data <- spatial_data[[2]]
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    # Crop the margins
    figure_margins <-
        raster::extent(
          xmin(sites) * 1.0003,
          xmax(sites) * 0.9997,
          ymin(sites) * 0.9997,
          ymax(sites) * 1.0003
        )
    cr <- crop(urban_data, figure_margins)
    site_points <- sites
    
    # !! Dodging the points a bit so that the sampling location is more on the point 
    # of the triangle point rather than the center, when plotted.
    range_factor <- ( max(site_points@data$lat) - min(site_points@data$lat) ) / 500
    site_points@data$lat <- site_points@data$lat+(range_factor*site_points@data$n)
    
    coords <- xyFromCell(cr, seq_len(ncell(cr)))
    urban_df <- stack(as.data.frame(getValues(cr)))
    names(urban_df) <- c('value', 'variable')
    urban_df <- cbind(coords, urban_df)
    
    # The "off map" area is coded as >100 %, so reset it to zero
    urban_df$value[(urban_df$value > 100)] <- 0
    
    gg <-
      ggplot() +
      geom_tile(
        data = urban_df,
        aes(
          x = x,
          y = y,
          fill = value,
          color = value
        )
      ) +
      coord_fixed(1.1) + #prevents weird warping of map due to tile size
      scale_fill_viridis(option = "A") +
      scale_color_viridis(option = "A") + #prevents white borders on tile
      geom_point(
        data = as.data.frame(site_points@data),
        aes(x = long, y = lat, size = n),
        color = "black",
        fill = "white",
        #size = 6,
        shape = 25
      ) +
      scale_size_binned_area(
        limits = c(0, 25),
        breaks = c(1, 5, 15)
      ) +
      theme_map() +
      ggtitle(label = as.character(as.factor(site_points@data$city[1]))) +
      guides(fill = guide_legend(title = "% Urban"),
             color = guide_legend(title = "% Urban")) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 20,
        vjust = 0
      )) 
    
    gg
    ggsave(paste0("figures/sampling_map/sites_", city, ".png"),
           dpi = "screen")
    return(gg)
  }


make_all_urban_site_plots <-
  function() {
    # This function serves as a wrapper to execute all of the cover plots
    
    g1 <- plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_BA.rds")
    ) + theme(plot.margin = unit(c(0, 0, -0.3, 0), "cm"))
    g2 <- plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_BO.rds")
    ) + theme(plot.margin = unit(c(0, 0, -0.3, 0), "cm"))
    g3 <- plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_LA.rds")
    ) + theme(plot.margin = unit(c(0, 0, -0.3, 0), "cm"))
    g4 <- plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_MN.rds")
    ) + theme(plot.margin = unit(c(-0.3, 0, 0, 0), "cm"))
    g5 <- plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_PX.rds")
    ) + theme(plot.margin = unit(c(-0.3, 0, 0, 0), "cm"))
    
    legend <- get_legend(
      # create some space to the left of the legend
      g1 + theme(legend.box.margin = unit(c(0, 0, 0, 2.5), "cm"),
                 legend.box = "horizontal")
    )
    
    mega_plot <- plot_grid(
      g1 + theme(legend.position="none"),
      g2 + theme(legend.position="none"),
      g3 + theme(legend.position="none"),
      g4 + theme(legend.position="none"),
      g5 + theme(legend.position="none"),
      legend,
      align = 'v',
      axis = "l",
      #hjust = -1,
      nrow = 2
    )
    mega_plot
    ggsave(paste0("figures/sampling_map/sites_ALL.png"),
           dpi = "screen")
    rm(g1,g2,g3,g4,g5,legend,mega_plot)
  }


