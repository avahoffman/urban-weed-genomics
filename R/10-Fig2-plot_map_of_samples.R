# This code creates a map of samples overlaid on NLCD urban cover data

library(raster)
library(ggplot2)
library(viridis)
library(cowplot) # theme_map()
source(file = "R/10-Fig2-wrangle-climate-normals.R")

plot_urban_cover_and_sites <-
  function(spatial_data, highres = F) {
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
    site_points <- sites
    urban_data <- spatial_data[[2]]
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    # Crop the margins
    figure_margins <-
      raster::extent(xmin(sites) * 1.0003,
                     xmax(sites) * 0.9997,
                     ymin(sites) * 0.9997,
                     ymax(sites) * 1.0003)
    cr <- crop(urban_data, figure_margins)
    
    # !! IMPORTANT -- this step aggregates pixels so it doesn't take a century
    # to plot. Aggregated data should be used for plotting only!
    if (!highres) {
      cr <- terra::aggregate(cr, fact = 6, fun = mean)
    }
    
    # !! Dodging the points a bit so that the sampling location is more on the point
    # of the triangle point rather than the center, when plotted.
    range_factor <- (max(site_points@data$lat) - min(site_points@data$lat)) / 500
    site_points@data$lat <- site_points@data$lat + (range_factor * site_points@data$n)
    
    coords <- xyFromCell(cr, seq_len(ncell(cr)))
    urban_df <- stack(as.data.frame(getValues(cr)))
    names(urban_df) <- c('value', 'variable')
    urban_df <- cbind(coords, urban_df)
    
    # The "off map" area is coded as >100 %, so reset it to zero
    urban_df$value[(urban_df$value > 100)] <- 0
    
    gg <-
      ggplot() +
      geom_tile(data = urban_df, aes(
        x = x,
        y = y,
        fill = value,
        color = value
      )) +
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
      scale_size_binned_area(limits = c(0, 25), breaks = c(5, 15)) +
      theme_classic() +
      labs(x = NULL, y = NULL) +
      ggtitle(label = as.character(as.factor(site_points@data$city[1]))) +
      guides(fill = guide_legend(title = "% Impervious"),
             color = guide_legend(title = "% Impervious")) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 15,
        vjust = 0
      ))
    
    gg
    # if(highres){
    #   ggsave(paste0("figures/Fig2_sampling_map/sites_", city, "_highres.png"),
    #          dpi = "print")
    # } else {
    #   ggsave(paste0("figures/Fig2_sampling_map/sites_", city, ".png"),
    #          dpi = "print")
    # }
    return(gg)
  }


make_all_urban_site_plots <-
  function() {
    # This function serves as a wrapper to execute all of the cover plots
    
    g1 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_BA.rds"))
    g2 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_BO.rds"))
    g3 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_LA.rds"))
    g4 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_MN.rds")) + theme(plot.margin = unit(c(0, 0, -0.3, 0), "cm")) +
      ggtitle(label = "Minneapolis-Saint Paul")
    g5 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_PX.rds")) + theme(plot.margin = unit(c(0, 0, -0.3, 0), "cm"))
    
    legend <- get_legend(# create some space to the left of the legend
      g3 + theme(
        legend.box.margin = unit(c(0, 0, 0, 2.5), "cm"),
        legend.box = "horizontal"
      ))
    
    mega_plot <- plot_grid(
      g1 + theme(legend.position = "none"),
      g2 + theme(legend.position = "none"),
      g3 + theme(legend.position = "none"),
      g4 + theme(legend.position = "none"),
      g5 + theme(legend.position = "none"),
      legend,
      align = 'v',
      axis = "l",
      #hjust = -1,
      ncol = 1
    )
    #mega_plot
    ggsave(
      paste0("figures/Fig2_sampling_map/sites_ALL.png"),
      dpi = "print",
      height = 13,
      width = 5,
      units = "in"
    )
    rm(g1, g2, g3, g4, g5, legend, mega_plot)
  }


# make_all_urban_site_plots_with_clim_normals <-
#   function() {
    g1 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_BA.rds")) + theme(plot.margin = unit(c(2, 2, 2, 15), "mm")) 
    g2 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_BO.rds")) + scale_y_continuous(breaks = c(42.0, 42.4)) + scale_x_continuous(breaks = c(-71.2,-70,-70.8)) + theme(plot.margin=unit(c(0,0,0,0),"mm"))
    g3 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_LA.rds")) + scale_y_continuous(breaks = c(33.8,34.0,34.2)) + scale_x_continuous(breaks = c(-118.4,-118.2)) 
    g4 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_MN.rds")) + ggtitle(label = "Minneapolis-Saint Paul")
    g5 <- plot_urban_cover_and_sites(readr::read_rds("spatial_data/trimmed_spatial_PX.rds")) + scale_y_continuous(breaks = c(33.4, 33.6)) + scale_x_continuous(breaks = c(-112, -111.8)) 
    
    legend1 <- get_legend(# create some space to the left of the legend
      g3 + theme(
        legend.box.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.box = "horizontal"
      ))
    
    gg1 <- get_climate_normals(city_ = "BA", temp_ = T) + xlab(NULL) + theme(plot.margin = unit(c(2, 2, 2, 15), "mm")) 
    gg2 <- get_climate_normals(city_ = "BO", temp_ = T) + ylab(NULL) + xlab(NULL) + theme(axis.text.y = element_blank())
    gg3 <- get_climate_normals(city_ = "LA", temp_ = T) + ylab(NULL) + xlab(NULL) + theme(axis.text.y = element_blank()) 
    gg4 <- get_climate_normals(city_ = "MN", temp_ = T) + ylab(NULL) + xlab(NULL) + theme(axis.text.y = element_blank()) 
    gg5 <- get_climate_normals(city_ = "PX", temp_ = T) + ylab(NULL) + xlab(NULL) + theme(axis.text.y = element_blank())
    
    gg6 <- get_climate_normals(city_ = "BA", temp_ = F) + xlab("") + theme(plot.margin = unit(c(2, 2, 2, 15), "mm")) 
    gg7 <- get_climate_normals(city_ = "BO", temp_ = F) + ylab(NULL) + xlab("") + theme(axis.text.y = element_blank())
    gg8 <- get_climate_normals(city_ = "LA", temp_ = F) + ylab(NULL) + theme(axis.text.y = element_blank())
    gg9 <- get_climate_normals(city_ = "MN", temp_ = F) + ylab(NULL) + xlab("") + theme(axis.text.y = element_blank())
    gg10 <- get_climate_normals(city_ = "PX", temp_ = F) + ylab(NULL) + xlab("") + theme(axis.text.y = element_blank()) 
    
    legend2 <- get_legend(gg3 + theme(legend.text = element_text(size = 10)))
    legend3 <- get_legend(gg8 + theme(legend.text = element_text(size = 10)))
    
    mega_plot <- plot_grid(
      g1 + theme(legend.position = "none"),
      g2 + theme(legend.position = "none"),
      g3 + theme(legend.position = "none"),
      g4 + theme(legend.position = "none"),
      g5 + theme(legend.position = "none"),
      legend1,
      gg1 + theme(legend.position = "none"),
      gg2 + theme(legend.position = "none"),
      gg3 + theme(legend.position = "none"),
      gg4 + theme(legend.position = "none"),
      gg5 + theme(legend.position = "none"),
      legend2,
      gg6 + theme(legend.position = "none"),
      gg7 + theme(legend.position = "none"),
      gg8 + theme(legend.position = "none"),
      gg9 + theme(legend.position = "none"),
      gg10 + theme(legend.position = "none"),
      legend3,
      # align = 'h',
      # axis = "b",
      #hjust = -1,
      nrow = 3,
      rel_heights = c(3, 2, 2),
      rel_widths = c(1.3, 1, 1, 1, 1, 0.7),
      labels = c("(a)", rep("", 5),
                 "(b)", rep("", 5),
                 "(c)", rep("", 5))
    )
    mega_plot
    ggsave(
      paste0("figures/Fig2_sampling_map/sites_climate_ALL.png"),
      dpi = "print",
      height = 6,
      width = 15,
      units = "in"
    )
  #   rm(g1, g2, g3, g4, g5, legend1, mega_plot)
  # }
