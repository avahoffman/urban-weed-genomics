library(polyRAD) # CRAN v2.0.0

spp_ <- "DS"

setwd("/scratch4/mavolio2/scratch_AH/")
samps <- read.table(paste0("popmap_", spp_, "_polyrad.txt"), sep = "\t")

#####

polyRAD_dat <- readRDS(paste0(spp_, "_polyRADdata.rds"))
hh <- HindHe(polyRAD_dat )
samps$HindHe <- rowMeans(hh, na.rm = TRUE)[samps$V1]
samps$Depth <- rowSums(polyRAD_dat$locDepth)[samps$V1]
write.csv(samps, paste0(spp_, "_HindHe.csv"))

#####

od <- TestOverdispersion(polyRAD_dat, to_test = 1:12)
saveRDS(od, paste0(spp_, "_overdispersion.rds"))
