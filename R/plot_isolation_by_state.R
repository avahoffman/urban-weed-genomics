# A note on genetic distance methods:
# We don't really know the extent to which drift plays a role in these species.
# Using a geometric approach (such as Roger's distance) has fewer assumptions.
# https://www.uwyo.edu/dbmcd/molmark/lect06/lect6.html has a nice summary of
# this. Also leverages the genpop distance functions. See https://www.rdocumentation.org/packages/adegenet/versions/2.1.8/topics/dist.genpop
#
# Installation note:
# SNPRelate needs gfortran to compile - see https://www.r-bloggers.com/2021/03/gfortran-support-for-r-on-macos-2/
#
# It was also stupidly hard to figure out how to read contents of gdsn objects
# used in SNPRelate because it uses a pointer system. Here's how to read sample IDs:
# read.gdsn(index.gdsn(genofile, path="sample.id"))

library(adegenet)
library(poppr)
library(stringr)
# BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)
library(ggrepel)
library(SNPRelate)

do_ibs <- function(infile_path, genind_spp, species_name) {
  # SNPRelate package used here
  # https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#identity-by-state-analysis
  # 
  
  # Read in file / open session
  #set.seed(999)
  genofile <-
    snpgdsOpen(infile_path,
               readonly = T,
               allow.duplicate = T)
  
  # Want to ensure the same samples are used in the genind object as this SNPRelate object.
  genind_samp <- rownames(genind_spp$tab)
  
  ibs <- snpgdsIBS(genofile, num.thread = 2,
                   sample.id = genind_samp)
  
  # Assume zero distance if data is missing
  ibs$ibs <- dplyr::na_if(ibs$ibs, "NaN") %>% replace(is.na(.), 0)
  
  # Replace sample name with actual city name
  genind_spp <- fix_pop_names(genind_spp)
  
  # Create clustering
  ibs.hc <- snpgdsHCluster(ibs)
  
  # ggplot prep
  hc <- ibs.hc$hclust
 # clus <- cutree(hc, 3)
  #g <- split(names(clus), clus)
  p <- ggtree(hc)
  #clades <- sapply(g, function(n) MRCA(p, n))
  #p <- groupClade(p, clades, group_name='subtree')
  d <- data.frame(label = hc$labels,
                  city = pop(genind_spp))

  # colors_ <- RColorBrewer::brewer.pal(n = 5, name = "Set3")
  colors_ <- viridis::viridis(n = 5, option = "H", begin = 0.2)
  my_pal <- setNames(colors_, 
           c("Baltimore", "Boston", "Los Angeles", "Minneapolis", "Phoenix"))
  
  # Create plot
  gg <- p %<+% d +
    layout_dendrogram() +
    geom_tippoint(aes(fill=factor(city), x=x),
                  size=2, shape=21, color='black') +
    scale_fill_manual(values = my_pal) +
    guides(fill = guide_legend(title = "")) +
    theme(legend.position="none") +
    ggtitle(paste0("      ", species_name)) 
  
  # Save plot
  ggsave(
    paste0(paste0("figures/trees/tree_", species_name, ".png")),
    dpi = "print",
    width = 12,
    height = 3
  )
  
  return(gg)
}


fix_pop_names <- function(x) {
  # x is a genind object

  # Replace sample name with actual city name
  pop(x) <-
    str_replace(pop(x), "(..).BA.(.+)", "Baltimore")
  pop(x) <-
    str_replace(pop(x), "(..).BO.(.+)", "Boston")
  pop(x) <-
    str_replace(pop(x), "(..).LA.(.+)", "Los Angeles")
  pop(x) <-
    str_replace(pop(x), "(..).MN.(.+)", "Minneapolis")
  pop(x) <-
    str_replace(pop(x), "(..).PX.(.+)", "Phoenix")

  return(x)
}


plot_dendrograms <- function(width = 12, height = 12) {
  # Read in data
  gen_ds <-
    read.genepop("SNP_data/DS/populations_20_pct/populations.snps.gen")
  gen_ls <-
    read.genepop("SNP_data/LS/populations_20_pct/populations.snps.gen")
  gen_pa <-
    read.genepop("SNP_data/PA/populations_20_pct/populations.snps.gen")
  gen_to <-
    read.genepop("SNP_data/TO/populations_20_pct/populations.snps.gen")
  gen_ec <-
    read.genepop("SNP_data/EC/populations_20_pct/populations.snps.gen")
  gen_cd <-
    read.genepop("SNP_data/CD/populations_20_pct/populations.snps.gen")
  
  # Create GDS files from the vcf files (only needs to be done once)
  if (!file.exists("SNP_data/DS/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/DS/populations_20_pct/populations.snps.vcf",
    "SNP_data/DS/populations_20_pct/snps.gds"
  )}
  if (!file.exists("SNP_data/LS/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/LS/populations_20_pct/populations.snps.vcf",
    "SNP_data/LS/populations_20_pct/snps.gds"
  )}
  if (!file.exists("SNP_data/PA/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/PA/populations_20_pct/populations.snps.vcf",
    "SNP_data/PA/populations_20_pct/snps.gds"
  )}
  if (!file.exists("SNP_data/TO/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/TO/populations_20_pct/populations.snps.vcf",
    "SNP_data/TO/populations_20_pct/snps.gds"
  )}
  if (!file.exists("SNP_data/CD/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/CD/populations_20_pct/populations.snps.vcf",
    "SNP_data/CD/populations_20_pct/snps.gds"
  )}
  if (!file.exists("SNP_data/EC/populations_20_pct/snps.gds")){
  snpgdsVCF2GDS(
    "SNP_data/EC/populations_20_pct/populations.snps.vcf",
    "SNP_data/EC/populations_20_pct/snps.gds"
  )}
  
  # Do Isolation by state
  p1 <- do_ibs(infile_path = "SNP_data/CD/populations_20_pct/snps.gds",
         genind_spp = gen_cd,
         species_name = "Bermuda grass")
  
  p2 <- do_ibs(infile_path = "SNP_data/DS/populations_20_pct/snps.gds",
         genind_spp = gen_ds,
         species_name = "crabgrass")

  p3 <- do_ibs(infile_path = "SNP_data/EC/populations_20_pct/snps.gds",
         genind_spp = gen_ec,
         species_name = "horseweed")
  
  p4 <- do_ibs(infile_path = "SNP_data/LS/populations_20_pct/snps.gds",
         genind_spp = gen_ls,
         species_name = "prickly lettuce")
  
  p5 <- do_ibs(infile_path = "SNP_data/PA/populations_20_pct/snps.gds",
         genind_spp = gen_pa,
         species_name = "bluegrass")
  
  p6 <- do_ibs(infile_path = "SNP_data/TO/populations_20_pct/snps.gds",
         genind_spp = gen_to,
         species_name = "dandelion")
  
  legend <- get_legend(
    # create some space to the left of the legend
    p6 + theme(legend.position = "bottom",
               legend.direction = "horizontal")
  )
  
  mega_plot <- plot_grid(
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    legend,
    align = 'v',
    axis = "l",
    #hjust = -1,
    ncol = 1,
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  )
  mega_plot
  setwd(here())
  ggsave(
    paste0("figures/trees/tree_all.png"),
    dpi = "print",
    width = width,
    height = height
  )
  
}
