#' This script extracts the data for easily making maps with NLCD urbanized area data types
#' Specifically, it trims the data so it isn't so dang big and unwieldy
#' 
#'
#' The U.S. Geological Survey (USGS), in partnership with several federal agencies,
#' has developed and released four National Land Cover Database (NLCD) products over
#' the past two decades: NLCD 1992, 2001, 2006, and 2011. This one is for data from
#' 2016 and describes urban imperviousness.
#'
#' https://www.mrlc.gov/data/type/urban-imperviousness
#'
#' NLCD imperviousness products represent urban impervious surfaces as a percentage of
#' developed surface over every 30-meter pixel in the United States. NLCD 2016 updates all
#' previously released versions of impervious products for CONUS (NLCD 2001, NLCD 2006,
#' NLCD 2011) along with a new date of impervious surface for 2016. New for NLCD 2016 is
#' an impervious surface descriptor layer. This descriptor layer identifies types of roads,
#' core urban areas, and energy production sites for each impervious pixel to allow deeper
#' analysis of developed features.
#
#' https://www.mrlc.gov/data/nlcd-2016-developed-imperviousness-descriptor-conus
###########################################################################################

library(raster)
library(stringr)
library(tidyr)
library(readr)
library(dplyr)

data_projection_settings <- 
  "+proj=longlat +datum=WGS84 +no_defs"


trim_site_data <-
  function(city = "None", filter_SNP = TRUE) {
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
      as.data.frame(read.csv(file = "data/site_data.csv",
                             header = T))
    
    if(filter_SNP){
      # Drop any sites from which all samples are filtered out. There are ~5
      # sites that had no remaining samples, so it doesn't make sense to plot
      # them..
      final_samples <- rbind(
        read.delim("SNP_data/CD/popmap_CD.txt", sep = '\t', header = F),
        read.delim("SNP_data/DS/popmap_DS.txt", sep = '\t', header = F),
        read.delim("SNP_data/EC/popmap_EC.txt", sep = '\t', header = F),
        read.delim("SNP_data/LS/popmap_LS.txt", sep = '\t', header = F),
        read.delim("SNP_data/PA/popmap_PA.txt", sep = '\t', header = F),
        read.delim("SNP_data/TO/popmap_TO.txt", sep = '\t', header = F)
      )
      # Pull out city, site, management combination
      final_samples$site <-
        stringr::str_extract(final_samples$V1, "(?<=...)[:graph:]{2,15}(?=\\..)")
      final_sites_to_use <- unique(final_samples$site)
      
      # Combine in site level data so it can be matched
      sitedata <-
        unite(
          sitedata,
          "site_comb",
          c(city_abbv, site_abbv, management_type),
          sep = ".",
          remove = F
        ) %>% 
        dplyr::mutate(site_comb = str_replace(site_comb, " ", "_"))
      
      sitedata <-
        filter(sitedata,
               site_comb %in% c(final_sites_to_use)) %>% 
        left_join(final_samples %>% group_by(site) %>% count(),
                  by = c("site_comb" = "site"))
      
      nrow(sitedata)
    }
    
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
  function(city = "None") {
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
    # r_BA <- trim_spatial(city = "BA")
    # ggplot_urban_pct_cover_plot(urban_data = r_BA, sites = s_BA,
    #    file_suffix = "urban_zone.jpg", guide_title = "% Urban")
    
    data_source <- "spatial_data/NLCD_2016_Impervious/NLCD_2016_Impervious_L48_20190405.img"
    
    # Read spatial data - check for raster or shapefile
    if (stringr::str_sub(data_source, -4, -1) == ".img") {
      r <-
        raster::raster(data_source)
    } else if (stringr::str_sub(data_source, -4, -1) == ".shp") {
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


create_spatial_rds_files <- function(){
  # This function creates the trimmed spatial data as rds
  # Make sure you set the R_MAX_VSIZE=100Gb in ~/.Renviron
  
  spatial_dat <- 
    "spatial_data/trimmed_spatial"
  
  s_BA <- trim_site_data(city = "BA")
  r_BA <- trim_spatial(city = "BA")
  readr::write_rds(list(s_BA, r_BA), paste0(spatial_dat, "_BA.rds"))
  
  s_BO <- trim_site_data(city = "BO")
  r_BO <- trim_spatial(city = "BO")
  readr::write_rds(list(s_BO, r_BO), paste0(spatial_dat, "_BO.rds"))
  
  s_MN <- trim_site_data(city = "MN")
  r_MN <- trim_spatial(city = "MN")
  readr::write_rds(list(s_MN, r_MN), paste0(spatial_dat, "_MN.rds"))
  
  s_PX <- trim_site_data(city = "PX")
  r_PX <- trim_spatial(city = "PX")
  readr::write_rds(list(s_PX, r_PX), paste0(spatial_dat, "_PX.rds"))
  
  s_LA <- trim_site_data(city = "LA")
  r_LA <- trim_spatial(city = "LA")
  readr::write_rds(list(s_LA, r_LA), paste0(spatial_dat, "_LA.rds"))
}
