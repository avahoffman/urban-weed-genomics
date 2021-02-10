###########################################################################################
# Execution script for project pipeline

###########################################################################################
# Set working directory for the repository (should be the git repo)
setwd(
  "/Users/avahoffman/Dropbox/Research/Urban_evo/Macrosystems/urban-weed-genomics/"
)

# Load source code
load("src/config.R")
load("src/sample_map.R")


s_BA <- trim_site_data(city = "BA")

r_BA <- trim_raster(city = "BA", data_source = urbanization_roads_data)
basic_urban_site_plot(urban_data = r_BA, sites = s_BA)
ggplot_urban_roads_plot(urban_data = r_BA, sites = s_BA, col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")

r_BA <- trim_raster(city = "BA", data_source = urbanization_cover_data)
basic_urban_site_plot(urban_data = r_BA, sites = s_BA)
ggplot_urban_pct_cover_plot(urban_data = r_BA, sites = s_BA, file_suffix = "urban_zone.jpg", guide_title = "% Urban")



s_BO <- trim_site_data(city = "BO")

r_BO <- trim_raster(city = "BO", data_source = urbanization_roads_data)
basic_urban_site_plot(urban_data = r_BO, sites = s_BO)
ggplot_urban_roads_plot(urban_data = r_BO, sites = s_BO, col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")

r_BO <- trim_raster(city = "BO", data_source = urbanization_cover_data)
basic_urban_site_plot(urban_data = r_BO, sites = s_BO)
ggplot_urban_pct_cover_plot(urban_data = r_BO, sites = s_BO, file_suffix = "urban_zone.jpg", guide_title = "% Urban")



s_MN <- trim_site_data(city = "MN")

r_MN <- trim_raster(city = "MN", data_source = urbanization_roads_data)
basic_urban_site_plot(urban_data = r_MN, sites = s_MN)
ggplot_urban_roads_plot(urban_data = r_MN, sites = s_MN, col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")

r_MN <- trim_raster(city = "MN", data_source = urbanization_cover_data)
basic_urban_site_plot(urban_data = r_MN, sites = s_MN)
ggplot_urban_pct_cover_plot(urban_data = r_MN, sites = s_MN, file_suffix = "urban_zone.jpg", guide_title = "% Urban")



s_PX <- trim_site_data(city = "PX")

r_PX <- trim_raster(city = "PX", data_source = urbanization_roads_data)
basic_urban_site_plot(urban_data = r_PX, sites = s_PX)
ggplot_urban_roads_plot(urban_data = r_PX, sites = s_PX, col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")

r_PX <- trim_raster(city = "PX", data_source = urbanization_cover_data)
basic_urban_site_plot(urban_data = r_PX, sites = s_PX)
ggplot_urban_pct_cover_plot(urban_data = r_PX, sites = s_PX, file_suffix = "urban_zone.jpg", guide_title = "% Urban")



s_LA <- trim_site_data(city = "LA")

r_LA <- trim_raster(city = "LA", data_source = urbanization_roads_data)
basic_urban_site_plot(urban_data = r_LA, sites = s_LA)
ggplot_urban_roads_plot(urban_data = r_LA, sites = s_LA, col_pal = NLCD_impervious_roads_color_key(), file_suffix = "urbanroads_ggplot.jpg")

r_LA <- trim_raster(city = "LA", data_source = urbanization_cover_data)
basic_urban_site_plot(urban_data = r_LA, sites = s_LA)
ggplot_urban_pct_cover_plot(urban_data = r_LA, sites = s_LA, file_suffix = "urban_zone.jpg", guide_title = "% Urban")
