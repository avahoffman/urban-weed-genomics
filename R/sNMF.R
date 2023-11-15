# sNMF alternative to Structure for validation

# BiocManager::install("LEA")
library(LEA)
library(ggplot2)
library(ggh4x) # facets on plots
library(cowplot)

convert_from_structure_for_sNMF <- function(){
  base_dir <- paste0(here::here(), "/SNP_data/")
  do_conversion <- function(input.file) {
    struct2geno(input.file = input.file, ploidy = 2, FORMAT = 2, extra.row = 0, extra.col = 2)
  }
  
  # EC - polyRAD
  input.file <- paste0(base_dir, "EC/EC_estimatedgeno_noheader.structure")
  do_conversion(input.file)
  # EC - SNPs from Stacks
  input.file <- paste0(base_dir, "EC/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # LS - polyRAD
  input.file <- paste0(base_dir, "LS/LS_estimatedgeno_noheader.structure")
  do_conversion(input.file)
  # LS - SNPs from Stacks
  input.file <- paste0(base_dir, "LS/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # DS - polyRAD
  # Note that I needed to use GenoDive to convert from n-ploidy lines per individual 
  # in the structure file to 2 lines per individual..
  input.file <- paste0(base_dir, "DS/DS_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # DS - SNPs from Stacks
  # Output from Stacks is already 2 lines per individual..
  input.file <- paste0(base_dir, "DS/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
}


do_sNMF_plot <- function(indata){
  gg <-
    ggplot(data = indata, aes(x = variable, y = value, fill = K_)) +
    geom_col(width = 1, color = NA) +
    facet_nested(~ city,
                 scales = "free_x",
                 space = "free",
                 switch = "both") +
    theme_classic() +
    # facetted_pos_scales(
    #   x = list(
    #     city == "Baltimore" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Baltimore"),]$x),
    #     city == "Boston" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Boston"),]$x),
    #     city == "Los Angeles" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Los Angeles"),]$x),
    #     city == "Minneapolis" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Minneapolis"),]$x),
    #     city == "*" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "*"),]$x),
    #     city == "Phoenix" ~ scale_x_discrete(position = "top", labels = x_lbl[(x_lbl$city == "Phoenix"),]$x)
    #   )
    # ) +
    #scale_fill_manual(values = c(my_pal)) +
  scale_fill_viridis_d(option = "E") +
    theme(
      legend.position = "none",
      # axis.text.x.top = element_text(
      #   angle = 90,
      #   hjust = 0,
      #   size = 5,
      #   vjust = 0.5
      # ),
      axis.text.x = element_blank(),
      panel.spacing = unit(0, "lines"),
      ggh4x.facet.nestline = element_line(linetype = 3),
      axis.ticks.length = unit(-1, "inch"),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0,5,5,5), "pt")
    )
  
  return(gg)
}


do_all_sNMF <- function(){
  base_dir <- paste0(here::here(), "/SNP_data/")
  
  # Do sNMF
  do_sNMF <- function(K_, infile_base){
    base_dir <- paste0(here::here(), "/SNP_data/")
    
    obj.snmf <- snmf(
      paste0(base_dir, infile_base, ".geno"),
      K = K_,
      alpha = 100,
      project = "new",
      entropy = TRUE
    )
    qmatrix <- Q(obj.snmf, K = K_)
    
    t_mat <- as.data.frame( t(qmatrix) )
    
    # Read in structure detail
    struct.detail <- read.table(paste0(base_dir, infile_base))[,1:2]
    ind_names <- struct.detail[rep(c(TRUE, FALSE), nrow(struct.detail)/2),]
    
    # Do some labeling according to structure output file
    colnames(t_mat) <- ind_names$V1
    t_mat$K_ <- rownames(t_mat)
    dat <- reshape2::melt(t_mat)
    
    # Pull out city data
    dat$city <- substr(dat$variable, 4, 5)
    
    # Reinstitute meaningful city names
    if (any(dat$city == "BA"))
      dat[dat$city == "BA", ]$city <- "Baltimore"
    if (any(dat$city == "BO"))
      dat[dat$city == "BO", ]$city <- "Boston"
    if (any(dat$city == "LA"))
      dat[dat$city == "LA", ]$city <- "Los Angeles"
    if (any(dat$city == "MN"))
      dat[dat$city == "MN", ]$city <- "Minneapolis"
    if (any(dat$city == "PX"))
      dat[dat$city == "PX", ]$city <- "Phoenix"
    
   
    return(dat) 
  }
  
  # EC
  sNMF_EC_K2 <- do_sNMF(K_ = 2, infile_base = "EC/EC_estimatedgeno_noheader.structure")
  p1 <- do_sNMF_plot(sNMF_EC_K2)
  #.LS
  sNMF_LS_K3 <- do_sNMF(K_ = 3, infile_base = "LS/LS_estimatedgeno_noheader.structure")
  p2 <- do_sNMF_plot(sNMF_LS_K3)
  # 
  
  mega_plot <- plot_grid(
    p1,
    p2,
    #align = 'h',
    #axis = "b",
    hjust = 0,
    ncol = 1
    #rel_widths = c(1, 0.2),
    # labels = c(
    #   "(a) Bermuda grass",
    #   "(b) crabgrass",
    #   "(c) horseweed",
    #   "(d) prickly lettuce",
    #   "(e) bluegrass",
    #   "(f) dandelion"
    # )
  )
  mega_plot

}