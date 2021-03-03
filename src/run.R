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
load("src/config.R")
load("src/utils.R")
load("src/NLCD_spatial_analysis.R")
load("src/redline_spatial_analysis.R")

# Generate urban percent cover by site (NLCD data)
make_site_urban_pct_csv()

# Make five plots containing urban cover and sampling sites by urban area
make_all_urban_cover_plots()

# Generate redlining data by site
write_site_level_redlining_to_csv(spatial_data = redlining_data)
