require(graph4lg)
require(tidyverse)
# library(popgraph)
# library(gstudio)
# library(igraph)

do_network_plot <- function(spp_, city_) {
  # Load in genind format data
  gen_ <-
    readRDS(paste0("SNP_data/", spp_, "/", spp_, "_estimatedgeno_genind.rds"))
  gen_ <- gen_[pop(gen_) == city_]
  
  # Extract site names (Will use these as pops)
  pop(gen_) <-
    paste0(pop(gen_),
           ".",
           str_replace(rownames(gen_$tab), "......", "") %>% str_replace("..$", ""))
  
  site_info <-
    read_csv("data/site_data_urban_cov.csv") %>%
    dplyr::select(site_abbv,
                  city_abbv,
                  management_type,
                  nlcd_urban_pct,
                  lat,
                  long) %>%
    mutate(site_abbv = str_replace(site_abbv, " ", "_"))  %>%
    unite(
      "site_abbv",
      city_abbv,
      site_abbv,
      management_type,
      sep = ".",
      remove = F
    ) %>%
    filter(city_abbv == city_, site_abbv %in% pop(gen_)) %>%
    select(site_abbv, long, lat) %>%
    rename(ID = site_abbv, x = long, y = lat) %>%
    as.data.frame()
  
  # Euclidean genetic distance computed in the same way as with the function popgraph
  # from popgraph package, i.e. after a PCA and two SVD, among other computation steps
  # (option ‘dist=PG’). This distance differs from the conditional genetic distance
  # (cGD) computed from a population graph by summing genetic distances along shortest paths.
  mat_pg <- mat_gen_dist(x = gen_, dist = "PG")
  mat_pg <- reorder_mat(mat_pg,
                        order = as.character(rownames(mat_pg))[order(as.character(rownames(mat_pg)))])
  mat_pg
  
  mat_geo <- mat_geo_dist(
    data = site_info,
    ID = "ID",
    x = "x",
    y = "y",
    crds_type = "polar"
  )
  mat_geo <- reorder_mat(mat_geo,
                         order = as.character(rownames(mat_geo))[order(as.character(rownames(mat_geo)))])
  mat_geo
  
  scatter_dist(mat_pg, mat_geo,
               pts_col = "black")
  
  ibd <- mantel.randtest(dist(mat_pg), dist(mat_geo), nrepet = 9999)
  ibd_tidy <- data.frame(
    "Observation" = ibd$obs,
    "Hypothesis" = ibd$alter,
    "Reps" = ibd$rep,
    "Std.Obs" = ibd$expvar[1],
    "Expectation" = ibd$expvar[2],
    "Variance" = ibd$expvar[3],
    "p-value" = ibd$pvalue,
    "spp" = spp_,
    "city" = city_
  )
  
  # Connect to 3 (or some other number) nearest neighbors
  # If 'topo = 'knn”, a k-nearest neighbor graph whose links are weighted with
  # values from 'mat_w' is created.
  graph_k3 <- gen_graph_topo(
    mat_w = mat_pg,
    mat_topo = mat_pg,
    topo = "knn",
    k = 3
  )
  
  # If 'topo = 'mst”, the resulting graph will have the topology of a minimum
  # spanning tree. It means that the graph will not include any cycle (tree) and
  # it will be the subgraph with a tree topology with the minimum total links'
  # weight (based on 'mat_topo' values).
  graph_mst <- gen_graph_topo(mat_w = mat_pg,
                              mat_topo = mat_pg,
                              topo = "mst")
  
  p <- plot_graph_lg(
    graph = graph_mst,
    mode = "spatial",
    crds = site_info,
    link_width = "inv_w"
  )
  # Adjust label size
  p[["layers"]][[3]][["aes_params"]][["size"]] <- 3
  p
  
  ggsave(
    paste0(
      "figures/IBD/distance_genetic_network_bycity_",
      spp_,
      "_",
      city_,
      ".png"
    ),
    dpi = "print",
    height = 7,
    width = 7,
    units = "in"
  )
  
  return(list(p, ibd_tidy))
}


do_by_city_ibd <- function() {
  g01 <- do_network_plot("CD", "BA")
  g02 <- do_network_plot("CD", "LA")
  g03 <- do_network_plot("CD", "PX")
  
  g04 <- do_network_plot("DS", "BA")
  g05 <- do_network_plot("DS", "BO")
  g06 <- do_network_plot("DS", "MN")
  g07 <- do_network_plot("DS", "PX")
  
  g08 <- do_network_plot("EC", "BA")
  g09 <- do_network_plot("EC", "LA")
  g10 <- do_network_plot("EC", "PX")
  
  g11 <- do_network_plot("LS", "BA")
  g12 <- do_network_plot("LS", "BO")
  g13 <- do_network_plot("LS", "LA")
  g14 <- do_network_plot("LS", "MN")
  g15 <- do_network_plot("LS", "PX")
  
  g16 <- do_network_plot("PA", "BA")
  g17 <- do_network_plot("PA", "BO")
  g18 <- do_network_plot("PA", "LA")
  g19 <- do_network_plot("PA", "PX")
  
  g20 <- do_network_plot("TO", "BA")
  g21 <- do_network_plot("TO", "BO")
  g22 <- do_network_plot("TO", "LA")
  g23 <- do_network_plot("TO", "MN")
  g24 <- do_network_plot("TO", "PX")
  
  all_ibd_stats <- rbind(
    g01[[2]],
    g02[[2]],
    g03[[2]],
    
    g04[[2]],
    g05[[2]],
    g06[[2]],
    g07[[2]],
    
    g08[[2]],
    g09[[2]],
    g10[[2]],
    
    g11[[2]],
    g12[[2]],
    g13[[2]],
    g14[[2]],
    g15[[2]],
    
    g16[[2]],
    g17[[2]],
    g18[[2]],
    g19[[2]],
    
    g20[[2]],
    g21[[2]],
    g22[[2]],
    g23[[2]],
    g24[[2]]
  )
  
  readr::write_csv(all_ibd_stats, "output/isolation-by-distance-by-city-mantel-test.csv")
  
}
