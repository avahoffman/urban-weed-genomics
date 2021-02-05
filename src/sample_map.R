# This script plots the site map with NLCD urbanized area data types
# 
# The U.S. Geological Survey (USGS), in partnership with several federal agencies, 
# has developed and released four National Land Cover Database (NLCD) products over 
# the past two decades: NLCD 1992, 2001, 2006, and 2011. This one is for data from
# 2016 and describes urban imperviousness.
# 
# https://www.mrlc.gov/data/nlcd-2016-developed-imperviousness-descriptor-conus
###########################################################################################
library(raster)
library(ggplot2)
library(cowplot)


NLCD_impervious_roads_color_key <-
  function() {
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


trim_raster <- 
  function(city = "BA", data_source) {
  # This function trims the giant NLCD spatial data to make it more manageable for plotting
  # Args:
  # city: string denoting the desired city. possible values: c("BA","BO","LA","MN","PX")
  
  r <-
    raster::raster(data_source)
  # Determine projection and dimensions of the spatial data
  # crs(r)
  # raster::extent(r)
  
  BA_crop <-
    raster::extent(1610000,
                   1680000,
                   1950005,
                   2000005)
  
  BO_crop <-
    raster::extent(1910000,
                   2080000,
                   2350005,
                   2500005)
  
  MN_crop <-
    raster::extent(130000,
                   290000,
                   2400005,
                   2500005)
  
  PX_crop <-
    raster::extent(-1570000,
                   -1400000,
                   1220005,
                   1350005)
  
  LA_crop <-
    raster::extent(-2070000,
                   -1900000,
                   1340005,
                   1520005)
  
  if (city == "BA") {
    city_raster <- crop(r, BA_crop)
  } else if (city == "BO") {
    city_raster <- crop(r, BO_crop)
  } else if (city == "MN") {
    city_raster <- crop(r, MN_crop)
  } else if (city == "PX") {
    city_raster <- crop(r, PX_crop)
  } else {
    city_raster <- crop(r, LA_crop)
  }
  
  plot(city_raster)
  
  return(city_raster)
  
}


trim_site_data <- 
  function(city = "BA") {
  # This function
  #
  sitedata <-
    as.data.frame(read.csv(file = field_site_data,
                           header = T))
  
  sd <-
    sitedata[(sitedata$city_abbv == city),]
  xy <-
    sd[, c("long", "lat")]
  
  spdf <-
    SpatialPointsDataFrame(
      coords = xy,
      data = sd,
      proj4string = CRS(
        '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
      )
    )
  
  return(spdf)
}


crop_raster_data <-
  function(urban_data, sites) {
    # This function..
    # 
    
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    spdf_transformed <-
      spTransform(sites, crs(urban_data))
    
    # Need to have different figure buffers if the x values are negative
    # (both PX and LA are negative on the AEA projection)
    if (city != "PX" & city != "LA") {
      figure_margins <-
        raster::extent(
          xmin(spdf_transformed) * 0.999,
          xmax(spdf_transformed) * 1.001,
          ymin(spdf_transformed) * 0.999,
          ymax(spdf_transformed) * 1.001
        )
    } else {
      figure_margins <-
        raster::extent(
          xmin(spdf_transformed) * 1.001,
          xmax(spdf_transformed) * 0.999,
          ymin(spdf_transformed) * 0.999,
          ymax(spdf_transformed) * 1.001
        )
    }
  
    return(
      list(
      "raster_" = crop(urban_data, figure_margins),
      "vector_" = spdf_transformed
      )
      )
  }


basic_urban_site_plot <-
  function(urban_data, sites) {
    # This function
    #
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
  # This function
  #
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
}


ggplot_urban_pct_cover_plot <- 
  function(urban_data, sites, file_suffix, guide_title = NA) {
    # This function
    #
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    cr <- crop_raster_data(urban_data, sites)$raster_
    site_points <- crop_raster_data(urban_data, sites)$vector_
    
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
  }
