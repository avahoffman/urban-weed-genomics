# This script works with redlining data
#
# https://dsl.richmond.edu/panorama/redlining/#loc=5/39.1/-94.58
#
# https://dsl.richmond.edu/panorama/redlining/#loc=4/39.069/-101.206&text=downloads
#
###########################################################################################
library(raster)


ggplot_redlining_plot <-
  function(spatial_data,
           city,
           file_suffix) {
    # Plot redlining polygons with sampling sites
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


make_all_redlining_plots <- 
  function(){
    # Wrapper
    # 
    ggplot_redlining_plot(spatial_data = redlining_data, city = "BA", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "BO", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "LA", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "MN", file_suffix = "_redlining.jpg")
    ggplot_redlining_plot(spatial_data = redlining_data, city = "PX", file_suffix = "_redlining.jpg")
    
  }