# This script plots the site map with NLCD urbanized area data types
###########################################################################################
library(raster)
library(ggplot2)
library(cowplot)


NLCD_impervious_color_key <-
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
  function(city = "BA") {
  # This function trims the giant NLCD spatial data to make it more manageable for plotting
  # Args:
  # city: string denoting the desired city. possible values: c("BA","BO","LA","MN","PX")
  
  r <-
    raster::raster(urbanization_cover_data)
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


basic_urban_site_plot <- 
  function(urban_data, sites) {
  # This function
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
  
  cr <- crop(urban_data, figure_margins)
  
  pdf(paste("figures/", city, "_basic.pdf", sep = ""))
  plot(cr)
  plot(spdf_transformed,
       add = TRUE,
       col = 'green',
       lwd = 5)
  dev.off()
  
}


ggplot_urban_site_plot <- 
  function(urban_data, sites) {
  # This function
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
  
  cr <- crop(urban_data, figure_margins)
  
  coords <- xyFromCell(cr, seq_len(ncell(cr)))
  urban_df <- stack(as.data.frame(getValues(cr)))
  names(urban_df) <- c('value', 'variable')
  urban_df <- cbind(coords, urban_df)
  urban_df$value <- as.factor(urban_df$value)
  
  urban_df <-
    dplyr::left_join(urban_df, NLCD_impervious_color_key())
  
  colz <-
    setNames(as.character(NLCD_impervious_color_key()$hex),
             NLCD_impervious_color_key()$lab)
  
  gg <-
    ggplot() +
    geom_tile(data = urban_df, aes(
      x = x,
      y = y,
      fill = lab,
      color = lab
    )) +
    scale_fill_manual(values = colz) +
    scale_color_manual(values = colz) +
    geom_point(
      data = as.data.frame(spdf_transformed@coords),
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
  ggsave(paste("figures/", city, "_ggplot.png", sep = ""),
  dpi = "print")
}


r_BA <- trim_raster(city = "BA")
s_BA <- trim_site_data(city = "BA")
#basic_urban_site_plot(r_BA, s_BA)
ggplot_urban_site_plot(r_BA, s_BA)

r_BO <- trim_raster(city = "BO")
s_BO <- trim_site_data(city = "BO")
#basic_urban_site_plot(r_BO, s_BO)
ggplot_urban_site_plot(r_BO, s_BO)

r_MN <- trim_raster(city = "MN")
s_MN <- trim_site_data(city = "MN")
#basic_urban_site_plot(r_MN, s_MN)
ggplot_urban_site_plot(r_MN, s_MN)

r_PX <- trim_raster(city = "PX")
s_PX <- trim_site_data(city = "PX")
#basic_urban_site_plot(r_PX, s_PX)
ggplot_urban_site_plot(r_PX, s_PX)

r_LA <- trim_raster(city = "LA")
s_LA <- trim_site_data(city = "LA")
#basic_urban_site_plot(r_LA, s_LA)
ggplot_urban_site_plot(r_LA, s_LA)