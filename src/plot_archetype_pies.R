# Make plots of "pie charts" on map. Often used to make inferences about 
# population structure in space.
# https://stackoverflow.com/questions/45504311/plot-geom-scatterpie-on-a-geom-tile-plot

library(ggnewscale)
library(sp)
library(raster)
library(cowplot)
library(viridis)
library(scatterpie)


plot_pies_on_map <- function(city, spp, k_, out_file, width = 10, height = 10) {
  # Get map data including sites and % urban cover for background
  primary_dir <- here::here()
  map_dat <-
    readr::read_rds(paste0(
      primary_dir,
      "/data/macrosystems_urban_spatial_data.rds_",
      city
    ))
  urban_data <- map_dat[[2]]
  sites <- map_dat[[1]]
  
  # Need to have different figure buffers if the x values are negative
  # (both PX and LA are negative on the AEA projection)
  # figure_margins <-
  #   raster::extent(
  #     xmin(sites) * 0.999,
  #     xmax(sites) * 1.001,
  #     ymin(sites) * 0.999,
  #     ymax(sites) * 1.001
  #   )
  figure_margins <-
    raster::extent(xmin(sites) * 1.0001,
                   xmax(sites) * 0.9999,
                   ymin(sites) * 0.9999,
                   ymax(sites) * 1.0001)
  
  # Prepare map data
  cr <- crop(urban_data, figure_margins)
  coords <- xyFromCell(cr, seq_len(ncell(cr)))
  urban_df <- stack(as.data.frame(getValues(cr)))
  names(urban_df) <- c('value', 'variable')
  urban_df <- cbind(coords, urban_df)
  
  # The "off map" area is coded as >100 %, so reset it to zero
  urban_df$value[(urban_df$value > 100)] <- 0
  urban_df$value_factor <- as.factor(as.integer(urban_df$value))
  
  # Prepare pie / genetics data
  spp_city <- paste0(spp, "_", city)
  pies <-
    readr::read_csv(
      paste0(
        primary_dir,
        "/data/SNP_data/by-city/",
        spp_city,
        "/populations_20_pct/Q_Data_",
        spp_city,
        ".csv"
      )
    )
  # Merge genetic data and site coordinates
  pies <-
    pies %>% filter(k == k_) %>% pivot_wider(id_cols = c(sample, site_abbv, city, nlcd_urban_pct, k))
  pies <-
    pies %>% left_join(as_tibble(sites), by = c("city", "site_abbv")) %>% dplyr::select(site_abbv, lat, long, starts_with("X"))
  # Collapse proportions across samples
  pies <-
    pies %>% group_by(site_abbv, lat, long) %>% summarise(across(starts_with("X"), ~ sum(.x)))
  
  # Get k category column names
  # Reverse it so that it uses the other end of the color spectrum
  archetype_cols <-
    colnames(pies)[!(colnames(pies) %in% c("site_abbv", "city", "lat", "long"))]
  
  # Set colors for pies and other params
  my_pal <- rev(RColorBrewer::brewer.pal(n = 5, name = "RdYlBu"))
  guide_title = "% Urban"
  
  gg <-
    ggplot() +
    #ggtitle(label = as.character(as.factor(sites@data$city[1]))) +
    coord_fixed(1) + #prevents weird warping of map due to tile size
    theme_map() +
    geom_tile(data = urban_df,
              aes(
                x = x,
                y = y,
                fill = value_factor,
                color = value_factor
              )) +
    scale_fill_viridis(option = "A",
                       discrete = T,
                       guide = "none") +
    scale_color_viridis(option = "A",
                        discrete = T,
                        guide = "none") + #prevents white borders on tile
    theme(plot.title = element_text(
      hjust = 0.5,
      size = 30,
      vjust = 0
    )) +
    #guides(fill = guide_legend(title = guide_title),
    #color = guide_legend(title = guide_title)) +
    new_scale_fill() +
    geom_scatterpie(aes(x = long, y = lat, group = site_abbv),
                    data = pies,
                    cols = archetype_cols) +
    scale_fill_manual(values = c(my_pal), guide = "none")
  
  gg
  
  ggsave(out_file,
         dpi = "print",
         width = width,
         height = height)
}
