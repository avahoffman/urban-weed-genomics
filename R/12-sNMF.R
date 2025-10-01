# sNMF alternative to Structure for validation

# BiocManager::install("LEA")
library(LEA)
library(ggplot2)
library(ggh4x) # facets on plots
library(cowplot)

convert_from_structure_for_sNMF <- function(){
  # Convert structure output from PolyRAD and Stacks to .geno format, compatible with sNMF
  # Data must be diploid in struct2geno function
  # Tried other methods to keep data polyploid but was unsuccessful. 
  
  base_dir <- paste0(here::here(), "/SNP_data/")
  do_conversion <- function(input.file) {
    struct2geno(input.file = input.file, ploidy = 2, FORMAT = 2, extra.row = 0, extra.col = 2)
  }
  
  # EC - polyRAD
  input.file <- paste0(base_dir, "EC/EC_estimatedgeno_noheader.structure")
  do_conversion(input.file)
  # EC - SNPs from Stacks
  # Need to create a no header version of the file
  # Make sure Los Angeles is actually Los.Angeles (no spaces allowed)
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
  # First, open the header-containing structure file from polyRAD in GenoDive, e.g. DS_estimatedgeno.structure
  # Go to Data > Transform and leave all params the same (diploid by default)
  # Go to File > Export and select structure (it will save as a .gdv file)
  # Use a text editor to strip the header and save with .structure, e.g., DS_estimatedgeno_noheader_genodivefix.structure
  input.file <- paste0(base_dir, "DS/DS_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # DS - SNPs from Stacks
  # Output from Stacks is already 2 lines per individual..
  input.file <- paste0(base_dir, "DS/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # CD - polyRAD
  input.file <- paste0(base_dir, "CD/CD_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # CD - SNPs from Stacks
  input.file <- paste0(base_dir, "CD/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # PA - polyRAD
  input.file <- paste0(base_dir, "PA/PA_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # PA - SNPs from Stacks
  input.file <- paste0(base_dir, "PA/populations_20_pct/populations_noheader.structure")
  do_conversion(input.file)
  
  # TO - polyRAD
  input.file <- paste0(base_dir, "TO/TO_estimatedgeno_noheader_genodivefix.structure")
  do_conversion(input.file)
  # TO - SNPs from Stacks
  input.file <- paste0(base_dir, "TO/populations_20_pct/populations_noheader.structure")
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
      #ggh4x.facet.nestline = element_line(linetype = 3),
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


do_all_sNMF <- function(width = 12, height = 12){
  base_dir <- paste0(here::here(), "/SNP_data/")
  
  # Do sNMF
  do_sNMF <- function(K_, infile_base) {
    base_dir <- paste0(here::here(), "/SNP_data/")
    reps_ <- 1
    outfile <- paste0(base_dir,
                      gsub("\\.structure", "", infile_base),
                      "_K",
                      K_,
                      "_snmfout.rds")
    
    if (!file.exists(outfile)) {
      # entropy = TRUE computes the cross-entropy criterion
      obj.snmf <- snmf(
        paste0(base_dir, infile_base, ".geno"),
        K = 2:10,
        project = "new",
        repetitions = 10,
        entropy = TRUE,
        ploidy = 2, # data coerced to be diploid..
        CPU = 4
      )
      
      # plot cross-entropy criterion for all runs in the object
      png(file = paste0(base_dir, infile_base, ".geno.K_cross_entropy.png"))
      plot(obj.snmf, pch = 19, cex = 1.2, main = paste(substr(infile_base, 1, 2), "cross-entropy plot by K"))
      dev.off()
      
      # Choose the best run based on cross entropy
      best_run <- which.min(cross.entropy(obj.snmf, K = K_))
      qmatrix <- Q(obj.snmf, K = K_, run = best_run)
      t_mat <- as.data.frame(t(qmatrix))
      
      # Read in structure detail
      struct.detail <- read.table(paste0(base_dir, infile_base))[, 1:2]
      ind_names <-
        struct.detail[rep(c(TRUE, FALSE), nrow(struct.detail) / 2), ]
      
      # Do some labeling according to structure output file
      colnames(t_mat) <- ind_names$V1
      t_mat$K_ <- rownames(t_mat)
      dat <- reshape2::melt(t_mat)
      
      # Pull out city data
      dat$city <- substr(dat$variable, 4, 5)
      
      # Reinstitute meaningful city names
      if (any(dat$city == "BA"))
        dat[dat$city == "BA",]$city <- "Baltimore"
      if (any(dat$city == "BO"))
        dat[dat$city == "BO",]$city <- "Boston"
      if (any(dat$city == "LA"))
        dat[dat$city == "LA",]$city <- "Los Angeles"
      if (any(dat$city == "MN"))
        dat[dat$city == "MN",]$city <- "Minneapolis"
      if (any(dat$city == "PX"))
        dat[dat$city == "PX",]$city <- "Phoenix"
      
      saveRDS(dat, outfile)
    } else {
      dat <- readRDS(outfile)
    }
    
    return(dat)
  }
  
  # Using the same K that were determined by Structure
  # EC
  sNMF_EC_K2 <- do_sNMF(K_ = 4, infile_base = "EC/EC_estimatedgeno_noheader.structure")
  p1 <- do_sNMF_plot(sNMF_EC_K2)
  # LS
  sNMF_LS_K3 <- do_sNMF(K_ = 3, infile_base = "LS/LS_estimatedgeno_noheader.structure")
  p2 <- do_sNMF_plot(sNMF_LS_K3)
  # CD
  sNMF_CD_K3 <- do_sNMF(K_ = 3, infile_base = "CD/CD_estimatedgeno_noheader_genodivefix.structure")
  p3 <- do_sNMF_plot(sNMF_CD_K3)
  # DS
  sNMF_DS_K3 <- do_sNMF(K_ = 3, infile_base = "DS/DS_estimatedgeno_noheader_genodivefix.structure")
  p4 <- do_sNMF_plot(sNMF_DS_K3)
  # PA
  sNMF_PA_K4 <- do_sNMF(K_ = 4, infile_base = "PA/PA_estimatedgeno_noheader_genodivefix.structure")
  p5 <- do_sNMF_plot(sNMF_PA_K4)
  # TO
  sNMF_TO_K3 <- do_sNMF(K_ = 4, infile_base = "TO/TO_estimatedgeno_noheader_genodivefix.structure")
  p6 <- do_sNMF_plot(sNMF_TO_K3)
  
  mega_plot <- plot_grid(
    p3 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    p4 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    p1 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    p2 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    p5 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    p6 + theme(plot.margin = unit(c(20, 5, 5, 5), "pt")),
    hjust = 0,
    ncol = 1,
    label_fontface = "italic",
    labels = c(
      "C. dactylon   K = 3",
      "D. sanguinalis   K = 3",
      "E. canadensis   K = 4",
      "L. serriola   K = 3",
      "P. annua   K = 4",
      "T. officinale   K = 4"
    )
  )
  mega_plot
  setwd(here::here())
  ggsave(
    paste0("figures/Fig5_structure/sNMF_all.png"),
    dpi = "print",
    width = width,
    height = height
  )
  
  quick_theme <- function(plot_) {
    return(plot_ +
             theme(
               plot.margin = unit(c(15, 5, 5, 5), "pt"),
               strip.text.x = element_text(size = 6)
             ))
  }
  
  # Plot for the supplemental material
  supp_plot <- plot_grid(
    quick_theme(p3),
    quick_theme(p4),
    quick_theme(p1),
    quick_theme(p2),
    quick_theme(p5),
    quick_theme(p6),
    hjust = 0,
    ncol = 1,
    label_size = 10,
    label_fontface = "italic",
    labels = c(
      "C. dactylon   K = 3",
      "D. sanguinalis   K = 3",
      "E. canadensis   K = 4",
      "L. serriola   K = 3",
      "P. annua   K = 4",
      "T. officinale   K = 4"
    )
  )
  supp_plot
  return(supp_plot)

}