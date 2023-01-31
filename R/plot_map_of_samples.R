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
          xmin(sites) * 1.0001,
          xmax(sites) * 0.9999,
          ymin(sites) * 0.9999,
          ymax(sites) * 1.0001
        )
    cr <- crop(urban_data, figure_margins)
    site_points <- sites
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
      coord_fixed(1.4) + #prevents weird warping of map due to tile size
      scale_fill_viridis(option = "A") +
      scale_color_viridis(option = "A") + #prevents white borders on tile
      geom_point(
        data = as.data.frame(site_points@coords),
        aes(x = long, y = lat),
        color = "black",
        fill = "white",
        size = 6,
        shape = 25
      ) +
      theme_map() +
      ggtitle(label = as.character(as.factor(site_points@data$city[1]))) +
      guides(fill = guide_legend(title = "% Urban"),
             color = guide_legend(title = "% Urban")) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 30,
        vjust = 0
      ))
    
    gg
    ggsave(paste0("figures/sampling_map/sites_", city, ".png"),
           dpi = "print")
    return(gg)
  }


make_all_urban_site_plots <-
  function() {
    # This function serves as a wrapper to execute all of the cover plots
    
    plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_BA.rds")
    )
    plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_BO.rds")
    )
    plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_LA.rds")
    )
    plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_MN.rds")
    )
    plot_urban_cover_and_sites(
      readr::read_rds("spatial_data/trimmed_spatial_PX.rds")
    )
  }
