#' Excellent resources for polyRAD:
#' Publication: https://academic.oup.com/g3journal/article/9/3/663/6026786
#' Documentation: https://cran.r-project.org/web/packages/polyRAD/polyRAD.pdf
#' Tutorial: https://lvclark.r-universe.dev/articles/polyRAD/polyRADtutorial.html
#' Tutorial: https://lvclark.r-universe.dev/articles/polyRADtutorials/population_genetics.html
#'  
library(polyRAD) # CRAN v2.0.0

setwd("/home/ahoffm31/scratch4-mavolio2/scratch_AH/polyRAD/PA_matches_for_polyrad/")

samps <- read.table("/scratch4/mavolio2/scratch_AH/polyRAD/popmap_PA_polyrad.txt", sep = "\t")
samps_crit <- ceiling(c(nrow(samps)*0.2)) # 20% threshold, round up to nearest integer

myStacks <- readStacks("catalog.alleles.tsv.gz",
                       ".",
                       version = 2,
                       taxaPloidy = 4,
                       min.ind.with.reads = samps_crit,
                       min.ind.with.minor.allele = 2,
                       readAlignmentData = F)

saveRDS(myStacks, "/scratch4/mavolio2/scratch_AH/polyRAD/PA_polyRADdata.rds")
