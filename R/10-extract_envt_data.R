

get_envt_data_from_site <-
  function(city) {
    # This function takes the NLCD percent cover data and finds the value for the
    # pixel that the sampling site is located within.
    #
    # Args:
    # urban_data: output from trim_spatial()
    # sites: output from trim_site_data(); used to narrow down spatial data for plotting
    #
    # Returns: site info dataframe with an additional column containing nlcd_urban_pct
    # field, collected from the NLCD data

    spatial_data <- readr::read_rds(paste0("spatial_data/trimmed_spatial_", city, ".rds"))
    urban_data <- spatial_data[[2]]
    sites <- spatial_data[[1]]
    
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
        dplyr::filter(abs(x - aea_site_points$long[i]) == min(abs(x - aea_site_points$long[i]))) %>%
        dplyr::filter(abs(y - aea_site_points$lat[i]) == min(abs(y - aea_site_points$lat[i]))) %>%
        dplyr::select(value)
      urban_cov[i, 1] <- n
    }
    colnames(urban_cov) <- "nlcd_urban_pct"
    
    return(cbind(sites@data, urban_cov))
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
    
    # Need figure buffers 
    figure_margins <-
      raster::extent(
        xmin(spdf_transformed) * 1.001,
        xmax(spdf_transformed) * 0.999,
        ymin(spdf_transformed) * 0.999,
        ymax(spdf_transformed) * 1.001
      )
    
    return(list(
      "raster_" = crop(urban_data, figure_margins),
      "vector_" = spdf_transformed
    ))
  }


write_urban_cover_data <- function(){
  urban_cov_data <- rbind(
    get_envt_data_from_site("BA"),
    get_envt_data_from_site("BO"),
    get_envt_data_from_site("LA"),
    get_envt_data_from_site("MN"),
    get_envt_data_from_site("PX")
  )
  readr::write_csv(urban_cov_data, "data/site_data_urban_cov.csv")
}
