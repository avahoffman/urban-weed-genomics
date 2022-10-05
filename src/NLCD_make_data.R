# This script extracts the data for easily making maps with NLCD urbanized area data types
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


crop_raster_data <-
  function(urban_data, sites) {
    # This function crops the raster data according to the sampling sites, plus a
    # little buffer for mapping purposes. Also transforms the sites so they can be
    # mapped together
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    #
    # Returns: list of ( raster_ = RasterLayer, vector_ = SpatialPointsDataFrame )
    #
    # Useage:
    # cr <- crop_raster_data(urban_data, sites)
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    # Transform sites to the same projection as
    # spdf_transformed <-
    #   spTransform(sites, crs(urban_data))
    
    spdf_transformed <- sites
    
    # Need to have different figure buffers if the x values are negative
    # (both PX and LA are negative on the AEA projection)
    if (city != "PX" | city != "LA") {
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
    
    return(list(
      "raster_" = crop(urban_data, figure_margins),
      "vector_" = spdf_transformed
    ))
  }


get_envt_data_from_site <-
  function(urban_data,
           sites) {
    # This function takes the NLCD percent cover data and finds the value for the
    # pixel that the sampling site is located within.
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    #
    # Returns: site info dataframe with an additional column containing nlcd_urban_pct
    # field, collected from the NLCD data
    #
    # Usage:
    # get_envt_data_from_site(urban_data = r_BA, sites = s_BA)
    
    # Automatically determine city from the input (sites)
    city <- as.character(as.factor(sites@data$city_abbv[1]))
    
    cr <- crop_raster_data(urban_data, sites)$raster_
    site_points <- crop_raster_data(urban_data, sites)$vector_
    
    coords <- xyFromCell(cr, seq_len(ncell(cr)))
    urban_df <- stack(as.data.frame(getValues(cr)))
    names(urban_df) <- c('value', 'variable')
    urban_df <- cbind(coords, urban_df)
    
    # The "off map" area (ocean, eg) is coded as >100 %, so reset it to zero
    urban_df$value[(urban_df$value > 100)] <- 0
    
    aea_site_points <- as.data.frame(site_points@coords)
    urban_cov <- data.frame()
    for (i in 1:nrow(aea_site_points)) {
      n <- urban_df %>%
        filter(abs(x - aea_site_points$long[i]) == min(abs(x - aea_site_points$long[i]))) %>%
        filter(abs(y - aea_site_points$lat[i]) == min(abs(y - aea_site_points$lat[i]))) %>%
        dplyr::select(value)
      urban_cov[i, 1] <- n
    }
    colnames(urban_cov) <- "nlcd_urban_pct"
    
    return(cbind(sites@data, urban_cov))
  }


make_all_urban_cover_data <-
  function() {
    # This function serves as a wrapper to execute all of the cover plots
    # Make sure you set the R_MAX_VSIZE=100Gb in ~/.Renviron
    
    s_BA <- trim_site_data(city = "BA")
    r_BA <-
      trim_spatial(city = "BA", data_source = urbanization_cover_data)
    readr::write_rds(list(s_BA, r_BA), paste0(nlcd_spatial_data, "_BA"))
    
    s_BO <- trim_site_data(city = "BO")
    r_BO <-
      trim_spatial(city = "BO", data_source = urbanization_cover_data)
    readr::write_rds(list(s_BO, r_BO), paste0(nlcd_spatial_data, "_BO"))
    
    s_MN <- trim_site_data(city = "MN")
    r_MN <-
      trim_spatial(city = "MN", data_source = urbanization_cover_data)
    readr::write_rds(list(s_MN, r_MN), paste0(nlcd_spatial_data, "_MN"))
    
    s_PX <- trim_site_data(city = "PX")
    r_PX <-
      trim_spatial(city = "PX", data_source = urbanization_cover_data)
    readr::write_rds(list(s_PX, r_PX), paste0(nlcd_spatial_data, "_PX"))
    
    s_LA <- trim_site_data(city = "LA")
    r_LA <-
      trim_spatial(city = "LA", data_source = urbanization_cover_data)
    readr::write_rds(list(s_LA, r_LA), paste0(nlcd_spatial_data, "_LA"))
    
    # Scale both datasets to the same projection, and extract pixel data
    BA_urbancov <-
      get_envt_data_from_site(urban_data = r_BA, sites = s_BA)
    BO_urbancov <-
      get_envt_data_from_site(urban_data = r_BO, sites = s_BO)
    MN_urbancov <-
      get_envt_data_from_site(urban_data = r_MN, sites = s_MN)
    PX_urbancov <-
      get_envt_data_from_site(urban_data = r_PX, sites = s_PX)
    LA_urbancov <-
      get_envt_data_from_site(urban_data = r_LA, sites = s_LA)
    
    # Concatenate the dataset
    d <-
      rbind(BA_urbancov,
            BO_urbancov,
            MN_urbancov,
            PX_urbancov,
            LA_urbancov)
    
    # Write combined data to file
    write.csv(d,
              file = nlcd_pct_urban_by_site_out,
              row.names = FALSE)
    
  }
