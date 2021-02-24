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
load("src/sample_map.R")

BA_urbancov <- get_envt_data_from_site(urban_data = r_BA, sites = s_BA)
BO_urbancov <- get_envt_data_from_site(urban_data = r_BO, sites = s_BO) #50% under 20
MN_urbancov <- get_envt_data_from_site(urban_data = r_MN, sites = s_MN) #43% under 20
PX_urbancov <- get_envt_data_from_site(urban_data = r_PX, sites = s_PX)
LA_urbancov <- get_envt_data_from_site(urban_data = r_LA, sites = s_LA)

d <- 
  rbind(BA_urbancov, BO_urbancov, MN_urbancov, PX_urbancov, LA_urbancov) %>%
  dplyr::select(city, urban_pct)

