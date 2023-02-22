library(polyRAD)

setwd("/data/mavolio2/catalogs_by_species/PA_metapop_catalog/")

samps <- read.table("/scratch4/mavolio2/scratch_AH/popmap_PA_polyrad.txt", sep = "\t")
samps_crit <- ceiling(c(nrow(samps)*0.2)) # 20% threshold, round up to nearest integer

myStacks <- readStacks("catalog.alleles.tsv.gz",
                       ".",
                       version = 2,
                       taxaPloidy = 4,
                       min.ind.with.reads = samps_crit,
                       min.ind.with.minor.allele = 2,
                       readAlignmentData = F)

saveRDS(myStacks, "/scratch4/mavolio2/scratch_AH/PA_polyRADdata.rds")