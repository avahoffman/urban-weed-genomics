library(polyRAD)

spp_ <- "PA"
inbreeding_ <- 0.3 # overdisp = 5

#####

mydata <- readRDS(paste0("/scratch4/mavolio2/scratch_AH/polyRAD/", spp_, "_polyRADdata.rds"))
od <- readRDS(paste0("/scratch4/mavolio2/scratch_AH/polyRAD/", spp_, "_overdispersion.rds"))
samps <- read.table(paste0(spp_, "_HindHe.csv"), sep = ",", header = T)

#####

hh <- HindHe(mydata)
hhByLoc <- colMeans(hh, na.rm = TRUE)

# Filter out reasonable frequencies
alfreq <- colMeans(mydata$depthRatio, na.rm = TRUE)
theseloci <- GetLoci(mydata)[mydata$alleles2loc[alfreq >= 0.05 & alfreq < 0.5]]
theseloci <- unique(theseloci)
hh_05 <- colMeans(hh[, theseloci], na.rm = TRUE)

print(mean(hh, na.rm = T))
print(mean(hh_05, na.rm = T))

expHindHe <- ExpectedHindHe(mydata, overdispersion = od$optimal, inbreeding = inbreeding_)
thresh_ <- quantile(expHindHe, probs = c(0.005, 0.995), na.rm = T)
keeploci <- names(hhByLoc)[ hhByLoc > thresh_[1] & hhByLoc < thresh_[2] ]
mydata2 <- SubsetByLocus(mydata, keeploci)
mydata2
saveRDS(mydata2, paste0("/scratch4/mavolio2/scratch_AH/polyRAD/", spp_, "_filtered_RADdata.rds"))

set.seed(420)
mydataPopStruct <- IteratePopStruct(mydata2, overdispersion = od$optimal)
mydata2_pca <- cbind(samps, mydataPopStruct$PCA[samps$V1,])
write.csv(mydata2_pca, paste0("/scratch4/mavolio2/scratch_AH/polyRAD/", spp_, "_IteratePopStructPCA.csv"))
saveRDS(mydataPopStruct, paste0("/scratch4/mavolio2/scratch_AH/polyRAD/", spp_, "_estimatedgeno_RADdata.rds"))

#####

