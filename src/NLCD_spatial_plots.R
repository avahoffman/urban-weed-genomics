# This script plots the site map with NLCD urbanized area data types
#
# The U.S. Geological Survey (USGS), in partnership with several federal agencies,
# has developed and released four National Land Cover Database (NLCD) products over
# the past two decades: NLCD 1992, 2001, 2006, and 2011. This one is for data from
# 2016 and describes urban imperviousness.
#
# https://www.mrlc.gov/data/type/urban-imperviousness
#
# NLCD imperviousness products represent urban impervious surfaces as a percentage of
# developed surface over every 30-meter pixel in the United States. NLCD 2016 updates all
# previously released versions of impervious products for CONUS (NLCD 2001, NLCD 2006,
# NLCD 2011) along with a new date of impervious surface for 2016. New for NLCD 2016 is
# an impervious surface descriptor layer. This descriptor layer identifies types of roads,
# core urban areas, and energy production sites for each impervious pixel to allow deeper
# analysis of developed features.
#
# https://www.mrlc.gov/data/nlcd-2016-developed-imperviousness-descriptor-conus
###########################################################################################
library(tidyverse)
library(raster)
library(ggplot2)
library(cowplot)
library(viridis)


NLCD_impervious_roads_color_key <-
  function() {
    # This function provides a hexidecimal color key for the urban roads data from NLCD.
    # I thought the colors were too bright so I made this one.
    #
    # Returns: data.frame
    #
    # Usage:
    # An argument to col_pal in ggplot_urban_roads_plot(), eg:
    # ggplot_urban_roads_plot(urban_data = r_BA, sites = s_BA,
    #    col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")
    
    return(data.frame(
      value = as.factor(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 127)),
      hex = c(
        '#000000',
        '#820000',
        '#490000',
        '#898905',
        '#464605',
        '#000082',
        '#000046',
        '#828282',
        '#464646',
        '#937176',
        '#894a89',
        '#004600',
        '#008200',
        '#000000'
      ),
      lab = c(
        'Natural area',
        'Primary road: urban',
        'Primary road: exurban',
        'Secondary road: urban',
        'Secondary road: exurban',
        'Tertiary road: urban',
        'Tertiary road: exurban',
        'Thinned road: urban',
        'Thinned road: exurban',
        'Impervious surface: urban',
        'Impervious surface: exurban',
        'Energy production: urban',
        'Energy production: exurban',
        'Natural area'
      )
    ))
  }


basic_urban_site_plot <-
  function(urban_data, sites) {
    # This function takes an urban environment raster layer and spatial data points and
    # makes a simple plot (can be used for debugging but doesn't look great)
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    #
    # Returns: plot to file
    #
    # Useage:
    # basic_urban_site_plot(urban_data = r_BA, sites = s_BA)
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    cr <- crop_raster_data(urban_data, sites)$raster_
    site_points <- crop_raster_data(urban_data, sites)$vector_
    
    pdf(paste("figures/", city, "_basic.pdf", sep = ""))
    plot(cr)
    plot(site_points,
         add = TRUE,
         col = 'green',
         lwd = 5)
    dev.off()
    dev.off()
    
  }


ggplot_urban_roads_plot <-
  function(urban_data, sites, col_pal, file_suffix) {
    # This function takes an urban roads raster layer and spatial data points and
    # makes a nice ggplot
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    # col_pal: output from NLCD_impervious_roads_color_key() -- can be tweaked if desired
    # file_suffix: suffix to add to any plot made -- to help flag what you plotted and also
    # allows you to change the file type (.jpg, .png, etc) created by ggsave()
    #
    # Returns: plot to file with file suffix
    #
    # Useage:
    # ggplot_urban_roads_plot(urban_data = r_BA, sites = s_BA,
    #   col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    cr <- crop_raster_data(urban_data, sites)$raster_
    site_points <- crop_raster_data(urban_data, sites)$vector_
    
    coords <- xyFromCell(cr, seq_len(ncell(cr)))
    urban_df <- stack(as.data.frame(getValues(cr)))
    names(urban_df) <- c('value', 'variable')
    urban_df <- cbind(coords, urban_df)
    urban_df$value <- as.factor(urban_df$value)
    
    urban_df <-
      dplyr::left_join(urban_df, col_pal)
    
    color_key <-
      setNames(as.character(col_pal$hex),
               col_pal$lab)
    
    gg <-
      ggplot() +
      geom_tile(data = urban_df, aes(
        x = x,
        y = y,
        fill = lab,
        color = lab
      )) +
      scale_fill_manual(values = color_key) +
      scale_color_manual(values = color_key) + #prevents white borders on tile
      geom_point(
        data = as.data.frame(site_points@coords),
        aes(x = long, y = lat),
        color = "green",
        size = 7,
        shape = '+'
      ) +
      theme_map() +
      ggtitle(label = as.character(as.factor(sites@data$city[1]))) +
      guides(fill = guide_legend(title = "Land type"),
             color = guide_legend(title = "Land type")) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 30,
        vjust = 0
      ))
    
    gg
    ggsave(paste("figures/", city, file_suffix, sep = ""),
           dpi = "print")
    
    return(gg)
  }


ggplot_urban_pct_cover_plot <-
  function(urban_data,
           sites,
           file_suffix,
           guide_title = NA) {
    # This function takes an urban pct cover raster layer and spatial data points and
    # makes a nice ggplot
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    # file_suffix: suffix to add to any plot made -- to help flag what you plotted and also
    # allows you to change the file type (.jpg, .png, etc) created by ggsave()
    # guide_title: Optional title for legend
    #
    # Returns: plot to file with file suffix
    #
    # Useage:
    # ggplot_urban_pct_cover_plot(urban_data = r_BA, sites = s_BA,
    # file_suffix = "urban_zone.jpg", guide_title = "% Urban")
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    # Need to have different figure buffers if the x values are negative
    # (both PX and LA are negative on the AEA projection)
    if (city == "BB") {
      figure_margins <-
        raster::extent(
          xmin(sites) * 0.999,
          xmax(sites) * 1.001,
          ymin(sites) * 0.999,
          ymax(sites) * 1.001
        )
    } else {
      figure_margins <-
        raster::extent(
          xmin(sites) * 1.0001,
          xmax(sites) * 0.9999,
          ymin(sites) * 0.9999,
          ymax(sites) * 1.0001
        )
    }
    
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
        color = "green",
        size = 7,
        shape = '+'
      ) +
      theme_map() +
      ggtitle(label = as.character(as.factor(sites@data$city[1]))) +
      guides(fill = guide_legend(title = guide_title),
             color = guide_legend(title = guide_title)) +
      theme(plot.title = element_text(
        hjust = 0.5,
        size = 30,
        vjust = 0
      ))
    
    gg
    ggsave(paste("figures/", city, file_suffix, sep = ""),
           dpi = "print")
    return(gg)
  }


make_all_urban_cover_plots <-
  function() {
    # This function serves as a wrapper to execute all of the cover plots
    
    BA <- readr::read_rds(paste0(nlcd_spatial_data, "_BA"))
    ggplot_urban_pct_cover_plot(
      urban_data = BA[[2]],
      sites = BA[[1]],
      file_suffix = "urban_zone.png",
      guide_title = "% Urban"
    )
    
    BO <- readr::read_rds(paste0(nlcd_spatial_data, "_BO"))
    ggplot_urban_pct_cover_plot(
      urban_data = BO[[2]],
      sites = BO[[1]],
      file_suffix = "urban_zone.png",
      guide_title = "% Urban"
    )
    
    MN <- readr::read_rds(paste0(nlcd_spatial_data, "_MN"))
    ggplot_urban_pct_cover_plot(
      urban_data = MN[[2]],
      sites = MN[[1]],
      file_suffix = "urban_zone.png",
      guide_title = "% Urban"
    )
    
    PX <- readr::read_rds(paste0(nlcd_spatial_data, "_PX"))
    ggplot_urban_pct_cover_plot(
      urban_data = PX[[2]],
      sites = PX[[1]],
      file_suffix = "urban_zone.png",
      guide_title = "% Urban"
    )
    
    LA <- readr::read_rds(paste0(nlcd_spatial_data, "_LA"))
    ggplot_urban_pct_cover_plot(
      urban_data = LA[[2]],
      sites = LA[[1]],
      file_suffix = "urban_zone.png",
      guide_title = "% Urban"
    )
    
  }
