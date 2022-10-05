###########################################################################################
# Execution script for project pipeline

# https://sedac.ciesin.columbia.edu/
# https://www.worldpop.org/geodata/summary?id=22518

###########################################################################################
# Set working directory for the repository (should be the git repo)
library(here)
setwd(
  here() # Make sure to launch directly from the src dirctory in order for here function to
         # work. Can also put path to the repo here instead.
)

# Load source code
source("src/config.R")
source("src/utils.R")
source("src/NLCD_make_data.R")
source("src/NLCD_spatial_plots.R")

# Create spatial data
# Creates a csv of urban cover by site as well as spatial rds to aid in plotting
make_all_urban_cover_data()

# Make five plots containing urban cover and sampling sites by urban area
make_all_urban_cover_plots()

# # Generate redlining data by site
# write_site_level_redlining_to_csv(spatial_data = redlining_data)
