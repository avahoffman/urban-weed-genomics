# Utility functions for the project
#

##########
# General utils

right_string <-
  function(x, n) {
    # Pulls n number of characters from the end of a string
    #
    # Args: x: string in question (e.g., file path)
    # n: number of characters back to go
    #
    substring(x, nchar(x) - n + 1)
  }


##########
# Spatial functions

trim_spatial <-
  function(city = "None", data_source) {
    # This function trims the giant NLCD spatial data to make it more manageable for plotting.
    # Can also trim and transform shapefiles.
    # Args:
    # city: string denoting the desired city. possible values: c("BA","BO","LA","MN","PX") or
    # "None" if no filtering is desired.
    # data_source: file path containing an .img file
    #
    # Returns: RasterLayer or Shapefile
    #
    # Useage:
    # Output of trim_spatial is a trimmed spatial dataset. It should be used as input urban_data
    # in plotting functions, eg:
    # r_BA <- trim_spatial(city = "BA", data_source = urbanization_cover_data)
    # ggplot_urban_pct_cover_plot(urban_data = r_BA, sites = s_BA,
    #    file_suffix = "urban_zone.jpg", guide_title = "% Urban")
    
    # Read spatial data - check for raster or shapefile
    if (right_string(data_source, 4) == ".img") {
      r <-
        raster::raster(data_source)
    } else if (right_string(data_source, 4) == ".shp") {
      r <-
        raster::shapefile(data_source)
    }
    
    # Determine projection and dimensions of the spatial data
    # CRS(r)
    # raster::extent(r)
    #
    
    # The following narrow down the coordinates to the approximate area around each city
    # (to make data processing more manageable)
    # Return only the city that is specified by the args
    if (city == "BA") {
      BA_crop <-
        raster::extent(1610000,
                       1680000,
                       1950005,
                       2003005)
      out <- crop(r, BA_crop)
    } else if (city == "BO") {
      BO_crop <-
        raster::extent(1910000,
                       2080000,
                       2350005,
                       2500005)
      out <- crop(r, BO_crop)
    } else if (city == "MN") {
      MN_crop <-
        raster::extent(130000,
                       290000,
                       2400005,
                       2500005)
      out <- crop(r, MN_crop)
    } else if (city == "PX") {
      PX_crop <-
        raster::extent(-1570000,
                       -1400000,
                       1220005,
                       1350005)
      out <- crop(r, PX_crop)
    } else if (city == "LA") {
      LA_crop <-
        raster::extent(-2070000,
                       -1900000,
                       1340005,
                       1520005)
      out <- crop(r, LA_crop)
    } else if (city == "None") {
      out <- r
    }
    
    # Ensure correct projection
    out_projected <-
      raster::projectRaster(from = out, crs = data_projection_settings)
    
    return(out_projected)
  }


trim_site_data <-
  function(city = "None") {
    # This function imports and prepares sampling site data
    # Args:
    # city: string denoting the desired city. possible values: c("BA","BO","LA","MN","PX")
    # or "None" if no filtering is desired
    #
    # Returns: SpatialPointsDataFrame
    #
    # Usage:
    # s_BA <- trim_site_data(city = "BA")
    
    # Read in site data (sampling points on the map)
    sitedata <-
      as.data.frame(read.csv(file = field_site_data,
                             header = T))
    
    # Filter and keep only the city specified in the args
    if (city != "None") {
      sd <-
        sitedata[(sitedata$city_abbv == city),]
    } else {
      sd <- sitedata
    }
    
    # Extract Long (x) and lat (y) for the dataframe below
    xy <-
      sd[, c("long", "lat")]
    
    # Format points to be in spatial format (so points and spatial can be plotted
    # together!)
    spdf <-
      SpatialPointsDataFrame(
        coords = xy,
        data = sd,
        proj4string = CRS(data_projection_settings)
      )
    
    return(spdf)
  }
