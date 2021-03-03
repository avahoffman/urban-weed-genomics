# This script works with redlining data
#
# https://dsl.richmond.edu/panorama/redlining/#loc=5/39.1/-94.58
#
# https://dsl.richmond.edu/panorama/redlining/#loc=4/39.069/-101.206&text=downloads
#
# Note that there is also data on population makeup of different ethnic groups in each city
# (see website downloads)
#
###########################################################################################
library(raster)


ggplot_redlining_plot <-
  function(spatial_data,
           city,
           file_suffix) {
    # Plot redlining polygons with sampling sites
    # Args:
    # spatial_data: file path containing an .shp file
    # city: string denoting the desired city. possible values: c("BA","BO","LA","MN","PX")
    # file_suffix: suffix to add to any plot made -- to help flag what you plotted and also
    # allows you to change the file type (.jpg, .png, etc) created by ggsave()
    # 
    # Returns: map image to file
    
    spatial <-
      trim_spatial(city = city, data_source = spatial_data)
    
    sites_ <- trim_site_data(city = city)
    site_points <- spTransform(sites_, crs(spatial))
    
    # Compile the data frame so that polygons can be colored by grade
    gg_df <-
      fortify(spatial) %>%
      mutate(piece = as.numeric(piece)) %>%
      left_join(data.frame(
        "grade" = spatial@data$holc_grade,
        "id" = row.names(spatial@data)
      ),
      by = "id")
    
    gg <- 
      ggplot() +
      
      # Spatial data
      geom_polygon(data = gg_df,
                   aes(
                     x = long,
                     y = lat,
                     group = group,
                     fill = grade
                   )) +
      
      # Site data points
      geom_point(
        data = as.data.frame(site_points@coords),
        aes(x = long, y = lat),
        color = "black",
        size = 7,
        shape = '+'
      )
    
    gg
    ggsave(paste("figures/", city, file_suffix, sep = ""),
           dpi = "print")
  }


write_site_level_redlining_to_csv <- 
  function(spatial_data){
    # Collect redlining data from the sampling locations
    # Args: file path containing an .shp file
    # 
    # Returns: Writes a .csv containing site info dataframe with additional 
    # columns from redlining data
    
    # Load data
    spatial <-
      trim_spatial(city = "None", data_source = spatial_data)
    sites_ <- trim_site_data(city = "None")
    
    # Transform to fit spatial projection
    site_points <- spTransform(sites_, crs(spatial))
    
    # Find the polygon (if any) in which the site lies
    dat <- sp::over(site_points, spatial)
    
    # Write combined data to file
    write.csv(cbind(site_points@data, dat),
              file = redlining_by_site_out,
              row.names = FALSE)
  }


make_all_redlining_plots <- 
  function(){
    # Wrapper to run all redlining polygon plots
    
    ggplot_redlining_plot(spatial_data = redlining_data, city = "BA", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "BO", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "LA", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "MN", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "PX", file_suffix = "_redlining.jpg")
    
  }