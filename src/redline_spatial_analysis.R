# This script works with redlining data
# 
# https://dsl.richmond.edu/panorama/redlining/#loc=5/39.1/-94.58
# 
# https://dsl.richmond.edu/panorama/redlining/#loc=4/39.069/-101.206&text=downloads
# 
###########################################################################################
library(raster)


s_BA <- trim_site_data(city = "BA")

city_raster <- trim_spatial(city = "BA", data_source = redlining_data)

site_points <- spTransform(s_BA, crs(city_raster))

gg_df <- 
  fortify(city_raster) %>%
  mutate(piece = as.numeric(piece)) %>%
  left_join(data.frame("grade" = city_raster@data$holc_grade,
           "id" = row.names(city_raster@data)),
           by = "id")

ggplot() +
  geom_polygon(data = gg_df,
               aes(
                 x = long,
                 y = lat,
                 group = group,
                 fill = grade
               )) +
  geom_point(
    data = as.data.frame(site_points@coords),
    aes(x = long, y = lat),
    color = "green",
    size = 7,
    shape = '+'
  ) 
  #scale_fill_manual(values = c("red","green","blue","cyan"))



outfff <- trim_raster(city = "BA", redlining_data)


aea_project_crs <- crs(r_BA)
save(aea_project_crs, file = data_projection_settings)
