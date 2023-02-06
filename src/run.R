###########################################################################################
# Execution script for project pipeline

# https://sedac.ciesin.columbia.edu/
# https://www.worldpop.org/geodata/summary?id=22518

###########################################################################################
# Set working directory for the repository (should be the git repo)
library(here)
setwd(
  here() # Make sure to launch directly from the src directory in order for here function to
         # work. Can also put path to the repo here instead.
)

# Load source code
source("src/config.R")
source("src/utils.R")
source("src/NLCD_make_data.R")
source("src/NLCD_spatial_plots.R")
source("src/archetype_analysis.R")
source("src/plot_archetypes_admixture.R")
source("src/plot_archetype_pies.R")

# Create spatial data
# Creates a csv of urban cover by site as well as spatial rds to aid in plotting
make_all_urban_cover_data()

# Make five plots containing urban cover and sampling sites by urban area
make_all_urban_cover_plots()

# # Generate redlining data by site
# write_site_level_redlining_to_csv(spatial_data = redlining_data)

# Run Archetype analysis - all metapopulations
run_archetype_analysis()

# Run Archetype analysis - by city and species for Baltimore only
run_archetype_analysis(spp = "CD", city = "BA")
run_archetype_analysis(spp = "DS", city = "BA")
run_archetype_analysis(spp = "EC", city = "BA")
run_archetype_analysis(spp = "LS", city = "BA")
run_archetype_analysis(spp = "PA", city = "BA")
run_archetype_analysis(spp = "TO", city = "BA")

# Run a few more
run_archetype_analysis(spp = "EC", city = "PX")
run_archetype_analysis(spp = "LS", city = "PX")
run_archetype_analysis(spp = "EC", city = "LA")
run_archetype_analysis(spp = "LS", city = "LA")

# Make Archetypal analysis plots - metapopulations
make_ADMIXTURE_plot(spp = "CD", 
                    k_list = c(3, 6),
                    out_file = "figures/genetics/CD_archetypes.png",
                    width = 10
)
make_ADMIXTURE_plot(spp = "DS", 
                    k_list = c(2, 4, 8),
                    out_file = "figures/genetics/DS_archetypes.png",
                    width = 10
)
make_ADMIXTURE_plot(spp = "EC", 
                    k_list = c(3, 6),
                    out_file = "figures/genetics/EC_archetypes.png",
                    width = 6
)
make_ADMIXTURE_plot(spp = "LS", 
                    k_list = c(5, 10),
                    out_file = "figures/genetics/LS_archetypes.png",
                    width = 10
)
make_ADMIXTURE_plot(spp = "PA", 
                    k_list = c(5, 10),
                    out_file = "figures/genetics/PA_archetypes.png",
                    width = 10
)
make_ADMIXTURE_plot(spp = "TO", 
                    k_list = c(5, 10),
                    out_file = "figures/genetics/TO_archetypes.png"
)


# Make Archetypal analysis plots - metapopulations but sorted by urbanness
make_ADMIXTURE_plot_by_urban(
  spp = "CD",
  k_list = c(2, 3, 4, 5),
  out_file = "figures/genetics/CD_archetypes_by_ULC.png",
  width = 10
)
make_ADMIXTURE_plot_by_urban(
  spp = "DS",
  k_list = c(2, 3, 4, 5),
  out_file = "figures/genetics/DS_archetypes_by_ULC.png",
  width = 10
)
make_ADMIXTURE_plot_by_urban(
  spp = "EC",
  k_list = c(2, 3, 4, 5),
  out_file = "figures/genetics/EC_archetypes_by_ULC.png",
  width = 6
)
make_ADMIXTURE_plot_by_urban(
  spp = "LS",
  k_list = c(2, 3, 4, 5),
  out_file = "figures/genetics/LS_archetypes_by_ULC.png",
  width = 10
)
make_ADMIXTURE_plot_by_urban(
  spp = "PA",
  k_list = c(2, 3, 4, 5),
  out_file = "figures/genetics/PA_archetypes_by_ULC.png",
  width = 10
)
make_ADMIXTURE_plot_by_urban(spp = "TO",
                             k_list = c(2, 3, 4, 5),
                             out_file = "figures/genetics/TO_archetypes_by_ULC.png")


# Make Archetypal analysis plots - by species for Baltimore only
make_ADMIXTURE_plot_by_City(spp = "CD",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/CD_BA_archetypes.png",
                            width = 4
                            )
make_ADMIXTURE_plot_by_City(spp = "DS",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/DS_BA_archetypes.png",
                            width = 4
)
make_ADMIXTURE_plot_by_City(spp = "EC",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/EC_BA_archetypes.png",
                            width = 3
)
make_ADMIXTURE_plot_by_City(spp = "LS",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/LS_BA_archetypes.png",
                            width = 3
)
make_ADMIXTURE_plot_by_City(spp = "PA",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/PA_BA_archetypes.png",
                            width = 4
)
make_ADMIXTURE_plot_by_City(spp = "TO",
                            city = "BA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/TO_BA_archetypes.png",
                            width = 5
)
make_ADMIXTURE_plot_by_City(spp = "EC",
                            city = "PX",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/EC_PX_archetypes.png",
                            width = 4
)
make_ADMIXTURE_plot_by_City(spp = "LS",
                            city = "PX",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/LS_PX_archetypes.png",
                            width = 5
)
make_ADMIXTURE_plot_by_City(spp = "EC",
                            city = "LA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/EC_LA_archetypes.png",
                            width = 4
)
make_ADMIXTURE_plot_by_City(spp = "LS",
                            city = "LA",
                            k_list = c(2, 3, 4, 5),
                            out_file = "figures/genetics/LS_LA_archetypes.png",
                            width = 5
)

# Plot pie graphs with archetypal analysis results for all k and species for Baltimore only
plot_pies_on_map(city = "BA", spp = "CD", k_ = 2, out_file = "figures/genetics/CD_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "CD", k_ = 3, out_file = "figures/genetics/CD_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "CD", k_ = 4, out_file = "figures/genetics/CD_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "CD", k_ = 5, out_file = "figures/genetics/CD_BA_pies_k5.png")

plot_pies_on_map(city = "BA", spp = "DS", k_ = 2, out_file = "figures/genetics/DS_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "DS", k_ = 3, out_file = "figures/genetics/DS_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "DS", k_ = 4, out_file = "figures/genetics/DS_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "DS", k_ = 5, out_file = "figures/genetics/DS_BA_pies_k5.png")

plot_pies_on_map(city = "BA", spp = "EC", k_ = 2, out_file = "figures/genetics/EC_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "EC", k_ = 3, out_file = "figures/genetics/EC_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "EC", k_ = 4, out_file = "figures/genetics/EC_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "EC", k_ = 5, out_file = "figures/genetics/EC_BA_pies_k5.png")

plot_pies_on_map(city = "BA", spp = "LS", k_ = 2, out_file = "figures/genetics/LS_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "LS", k_ = 3, out_file = "figures/genetics/LS_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "LS", k_ = 4, out_file = "figures/genetics/LS_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "LS", k_ = 5, out_file = "figures/genetics/LS_BA_pies_k5.png")

plot_pies_on_map(city = "BA", spp = "PA", k_ = 2, out_file = "figures/genetics/PA_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "PA", k_ = 3, out_file = "figures/genetics/PA_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "PA", k_ = 4, out_file = "figures/genetics/PA_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "PA", k_ = 5, out_file = "figures/genetics/PA_BA_pies_k5.png")

plot_pies_on_map(city = "BA", spp = "TO", k_ = 2, out_file = "figures/genetics/TO_BA_pies_k2.png")
plot_pies_on_map(city = "BA", spp = "TO", k_ = 3, out_file = "figures/genetics/TO_BA_pies_k3.png")
plot_pies_on_map(city = "BA", spp = "TO", k_ = 4, out_file = "figures/genetics/TO_BA_pies_k4.png")
plot_pies_on_map(city = "BA", spp = "TO", k_ = 5, out_file = "figures/genetics/TO_BA_pies_k5.png")

plot_pies_on_map(city = "PX", spp = "EC", k_ = 2, out_file = "figures/genetics/EC_PX_pies_k2.png")
plot_pies_on_map(city = "PX", spp = "EC", k_ = 3, out_file = "figures/genetics/EC_PX_pies_k3.png")
plot_pies_on_map(city = "PX", spp = "EC", k_ = 4, out_file = "figures/genetics/EC_PX_pies_k4.png")
plot_pies_on_map(city = "PX", spp = "EC", k_ = 5, out_file = "figures/genetics/EC_PX_pies_k5.png")

plot_pies_on_map(city = "PX", spp = "LS", k_ = 2, out_file = "figures/genetics/LS_PX_pies_k2.png")
plot_pies_on_map(city = "PX", spp = "LS", k_ = 3, out_file = "figures/genetics/LS_PX_pies_k3.png")
plot_pies_on_map(city = "PX", spp = "LS", k_ = 4, out_file = "figures/genetics/LS_PX_pies_k4.png")
plot_pies_on_map(city = "PX", spp = "LS", k_ = 5, out_file = "figures/genetics/LS_PX_pies_k5.png")

plot_pies_on_map(city = "LA", spp = "EC", k_ = 2, out_file = "figures/genetics/EC_LA_pies_k2.png")
plot_pies_on_map(city = "LA", spp = "EC", k_ = 3, out_file = "figures/genetics/EC_LA_pies_k3.png")
plot_pies_on_map(city = "LA", spp = "EC", k_ = 4, out_file = "figures/genetics/EC_LA_pies_k4.png")
plot_pies_on_map(city = "LA", spp = "EC", k_ = 5, out_file = "figures/genetics/EC_LA_pies_k5.png")

plot_pies_on_map(city = "LA", spp = "LS", k_ = 2, out_file = "figures/genetics/LS_LA_pies_k2.png")
plot_pies_on_map(city = "LA", spp = "LS", k_ = 3, out_file = "figures/genetics/LS_LA_pies_k3.png")
plot_pies_on_map(city = "LA", spp = "LS", k_ = 4, out_file = "figures/genetics/LS_LA_pies_k4.png")
plot_pies_on_map(city = "LA", spp = "LS", k_ = 5, out_file = "figures/genetics/LS_LA_pies_k5.png")
