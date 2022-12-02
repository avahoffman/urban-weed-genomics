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
# BiocManager::install("ggtree")
library(ggtree)
library(ggrepel)
library(SNPRelate)


do_phylogeny <- function(genind_spp, method, title) {
  # Replace sample name with actual city name
  genind_spp <- fix_pop_names(genind_spp)
  
  # Create a bootstrapped tree
  regionaltree <- genind_spp %>%
    genind2genpop(pop = ~ Pop) %>%
    aboot(
      sample = 1000,
      distance = rogers.dist,
      cutoff = 0,
      quiet = F,
      tree = 'nj'
    )
  
  ## treat each pop individual
  localcls <- as.list(unique(genind_spp$pop))
  regionaltree <- groupOTU(regionaltree, localcls)
  
  # Create ggtree plot
  ggtree(regionaltree, layout = "rectangular") +
    geom_label_repel(
      aes(label = label),
      force = 0,
      nudge_x = 0,
      nudge_y = 0
    ) +
    theme_tree() +
    ggtitle(title)
  
  # Save file
  ggsave(
    file = paste("figures/trees/tree", runif(1), ".jpg", sep = ""),
    height = 4,
    width = 5
  )
  
  # Repeat with mean imputed for missing data
  genind_spp$tab <- tab(genind_spp, NA.method = "mean")
  
  # Create a bootstrapped tree
  regionaltree <- genind_spp %>%
    genind2genpop(pop = ~ Pop) %>%
    aboot(
      sample = 1000,
      distance = rogers.dist,
      cutoff = 0,
      quiet = F,
      tree = method
    )
  
  ## treat each pop individual
  localcls <- as.list(unique(genind_spp$pop))
  regionaltree <- groupOTU(regionaltree, localcls)
  
  # Create ggtree plot
  ggtree(regionaltree, layout = "rectangular") +
    geom_label_repel(
      aes(label = label),
      force = 0,
      nudge_x = 0,
      nudge_y = 0
    ) +
    theme_tree() +
    ggtitle(paste(title, " - mean imputed for missing data", sep = ""))
  
  # Save file
  ggsave(
    file = paste("figures/trees/tree", runif(1), ".jpg", sep = ""),
    height = 4,
    width = 5
  )
}


do_ibs <- function(infile_path, genind_spp, spp) {
  # SNPRelate package used here
  # https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#identity-by-state-analysis
  # 
  
  # Read in file / open session
  set.seed(999)
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
  rv2 <- snpgdsCutTree(ibs.hc, samp.group = pop(genind_spp))
  
  grDevices::jpeg(
    file = paste("figures/trees/tree", runif(1), ".jpg", sep = ""),
    width = 7,
    height = 4,
    units = "in",
    res = 1000
  )
  plot(
    rv2$dendrogram,
    leaflab = "none",
    main = paste0(spp, ": Isolation-by-State (SNPRelate)")
  )
  legend(
    "topright",
    legend = levels(pop(genind_spp)),
    col = 1:nlevels(pop(genind_spp)),
    pch = 19,
    ncol = 4
  )
  dev.off()
  snpgdsClose(genofile)
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


run_phylogenies <- function() {
  # Read in data
  gen_ds <-
    read.genepop("data/SNP_data/by-species/DS/populations_20_pct/populations.snps.gen")
  gen_ls <-
    read.genepop("data/SNP_data/by-species/LS/populations_20_pct/populations.snps.gen")
  gen_pa <-
    read.genepop("data/SNP_data/by-species/PA/populations_20_pct/populations.snps.gen")
  gen_to <-
    read.genepop("data/SNP_data/by-species/TO/populations_20_pct/populations.snps.gen")
  gen_ec <-
    read.genepop("data/SNP_data/by-species/EC/populations_20_pct/populations.snps.gen")
  gen_cd <-
    read.genepop("data/SNP_data/by-species/CD/populations_20_pct/populations.snps.gen")
  
  # Create GDS files from the vcf files (only needs to be done once)
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/DS/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/DS/populations_20_pct/snps.gds"
  # )
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/LS/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/LS/populations_20_pct/snps.gds"
  # )
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/PA/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/PA/populations_20_pct/snps.gds"
  # )
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/TO/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/TO/populations_20_pct/snps.gds"
  # )
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/CD/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/CD/populations_20_pct/snps.gds"
  # )
  # snpgdsVCF2GDS(
  #   "data/SNP_data/by-species/EC/populations_20_pct/populations.snps.vcf",
  #   "data/SNP_data/by-species/EC/populations_20_pct/snps.gds"
  # )
  
  # can't do Erigeron or Cynodon at the population level.. because need at least 4 cities.
  do_phylogeny(gen_ds, method = "nj", title = "Digitaria")
  do_phylogeny(gen_ls, method = "nj",title = "Lactuca")
  do_phylogeny(gen_pa[-c(141, 142)], method = "nj", title = "Poa") # Drop MN because there are only 2 samples
  do_phylogeny(gen_to, method = "nj", title = "Taraxicum")
  
  # Do Isolation by state
  do_ibs(infile_path = "data/SNP_data/by-species/DS/populations_20_pct/snps.gds",
         genind_spp = gen_ds,
         spp = "Digitaria")
  do_ibs(infile_path = "data/SNP_data/by-species/CD/populations_20_pct/snps.gds",
         genind_spp = gen_cd,
         spp = "Cynodon")
  do_ibs(infile_path = "data/SNP_data/by-species/EC/populations_20_pct/snps.gds",
         genind_spp = gen_ec,
         spp = "Erigeron")
  do_ibs(infile_path = "data/SNP_data/by-species/LS/populations_20_pct/snps.gds",
         genind_spp = gen_ls,
         spp = "Lactuca")
  do_ibs(infile_path = "data/SNP_data/by-species/PA/populations_20_pct/snps.gds",
         genind_spp = gen_pa,
         spp = "Poa")
  do_ibs(infile_path = "data/SNP_data/by-species/TO/populations_20_pct/snps.gds",
         genind_spp = gen_to,
         spp = "Taraxicum")
  
}
