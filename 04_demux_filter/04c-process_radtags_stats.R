# Process-radtags----

# Plot process_radtags per-library and per-sample statistics to identify and
# remove outliers

## Import excel spreadsheet with data and rename
# Checking_process_radtags <- readr::read_csv("output/process_radtags-sample_output.csv")
processRadtags <- Checking_process_radtags

## Make into dataframe
processRadtags <-data.frame(processRadtags)

## Add a species and city identifier
processRadtags <- processRadtags %>%
    mutate(species = case_when (
      startsWith(Filename, "DS") ~ "DS",
      startsWith(Filename, "CD") ~ "CD",
      startsWith(Filename, "PA") ~ "PA",
      startsWith(Filename, "PA") ~ "PA",
      startsWith(Filename, "LS") ~ "LS",
      startsWith(Filename, "TO") ~ "TO",
      startsWith(Filename, "TE") ~ "TE"))

processRadtags <- processRadtags %>%
  mutate(city = case_when (
    str_detect(Filename, "BA") ~ "Baltimore",
    str_detect(Filename, "LA") ~ "LosAngeles",
    str_detect(Filename, "PX") ~ "Phoenix",
    str_detect(Filename, "BO") ~ "Boston",
    str_detect(Filename, "MN") ~ "Minneapolis"))

## Look at the proportion of each sublibrary that each sample makes up.
## I've calculated this by simply dividing the total reads for each sample by
## the total reads for each relevant sublibrary. I've plotted all samples here
## by species (not by sublibraries), as it doesn't matter exactly which sublibraries
## each sample is a part of but serves as a check as to whether a specific sample
## was sequenced proportionally low. 

### Create different datasets for each species
DS_processRadtags <- subset(processRadtags,species == "DS")
CD_processRadtags <- subset(processRadtags,species == "CD")
PA_processRadtags <- subset(processRadtags,species == "PA")
PA_processRadtags <- subset(processRadtags,species == "PA")
LS_processRadtags <- subset(processRadtags,species == "LS")
TO_processRadtags <- subset(processRadtags,species == "TO")
TE_processRadtags <- subset(processRadtags,species == "TE")


## Plots are shown for each species because variable limits may differ per species. 
## This code could potentially be looped so it is done automatically for each species. 

## _CD----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(CD_processRadtags$prop_sample_per_library)
sd(CD_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
CD_processRadtags <- CD_processRadtags[order(CD_processRadtags$prop_sample_per_library),]

ppi <- 300
png("CD_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,10,1.5,1.5) + 0.5)
dotchart(CD_processRadtags$prop_sample_per_library,
         labels = CD_processRadtags$Filename,
         xlab='CD: Proportion of Library',
         xlim=c(0,0.9),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(CD_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("CD_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=CD_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,0.9))
title(xlab='CD: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaCD
CD_processRadtags <- CD_processRadtags[order(CD_processRadtags$retained_reads),]

#### Descriptive statistics
summary(CD_processRadtags$retained_reads)
sd(CD_processRadtags$retained_reads)

#### Plot dotchart
png("CD_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(CD_processRadtags$retained_reads/1e6, 
         labels = CD_processRadtags$Filename,
         xlab='CD: Retained Reads (M)',
         xlim=c(0,60),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(CD_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("CD_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=CD_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,60))
title(xlab='CD: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
CD_processRadtags <- CD_processRadtags[order(CD_processRadtags$prop_reads_retained_per_sample),]

#### Descriptive statistics
summary(CD_processRadtags$prop_reads_retained_per_sample)
sd(CD_processRadtags$prop_reads_retained_per_sample)

#### Plot Dotchatt
png("CD_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(CD_processRadtags$prop_reads_retained_per_sample, 
         labels = CD_processRadtags$Filename,
         xlab='CD: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(CD_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("CD_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=CD_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='CD: Proportion of Retained Reads', line=2)
title(ylab='City', line=4.3)

.=dev.off()

## _DS----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(DS_processRadtags$prop_sample_per_library)
sd(DS_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
DS_processRadtags <- DS_processRadtags[order(DS_processRadtags$prop_sample_per_library),]

ppi <- 300
png("DS_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,10,1.5,1.5) + 0.5)
dotchart(DS_processRadtags$prop_sample_per_library,
         labels = DS_processRadtags$Filename,
         xlab='DS: Proportion of Library',
         xlim=c(0,0.9),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(DS_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("DS_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=DS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,0.9))
title(xlab='DS: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaDS
DS_processRadtags <- DS_processRadtags[order(DS_processRadtags$retained_reads),]

#### Descriptive statistics
summary(DS_processRadtags$retained_reads)
sd(DS_processRadtags$retained_reads)

#### Plot dotchart
png("DS_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(DS_processRadtags$retained_reads/1e6, 
         labels = DS_processRadtags$Filename,
         xlab='DS: Retained Reads (M)',
         xlim=c(0,60),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(DS_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("DS_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=DS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,60))
title(xlab='DS: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
DS_processRadtags <- DS_processRadtags[order(DS_processRadtags$prop_reads_retained_per_sample),]

#### Descriptive statistics
summary(DS_processRadtags$prop_reads_retained_per_sample)
sd(DS_processRadtags$prop_reads_retained_per_sample)

#### Plot Dotchatt
png("DS_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(DS_processRadtags$prop_reads_retained_per_sample, 
         labels = DS_processRadtags$Filename,
         xlab='DS: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(DS_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("DS_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=DS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='DS: Proportion of Retained Reads', line=2)
title(ylab='City', line=4.3)

.=dev.off()


## _EC----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(EC_processRadtags$prop_sample_per_library)
sd(EC_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
EC_processRadtags <- EC_processRadtags[order(EC_processRadtags$prop_sample_per_library),]

ppi <- 300
png("EC_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,10,1.5,1.5) + 0.5)
dotchart(EC_processRadtags$prop_sample_per_library,
         labels = EC_processRadtags$Filename,
         xlab='EC: Proportion of Library',
         xlim=c(0,0.7),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(EC_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("EC_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=EC_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,0.9))
title(xlab='EC: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaEC
EC_processRadtags <- EC_processRadtags[order(EC_processRadtags$retained_reads),]

#### Descriptive statistics
summary(EC_processRadtags$retained_reads)
sd(EC_processRadtags$retained_reads)

#### Plot dotchart
png("EC_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(EC_processRadtags$retained_reads/1e6, 
         labels = EC_processRadtags$Filename,
         xlab='EC: Retained Reads (M)',
         xlim=c(0,50),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(EC_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("EC_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=EC_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,50))
title(xlab='EC: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
EC_processRadtags <- EC_processRadtags[order(EC_processRadtags$prop_reads_retained_per_sample),]

#### Descriptive statistics
summary(EC_processRadtags$prop_reads_retained_per_sample)
sd(EC_processRadtags$prop_reads_retained_per_sample)

#### Plot Dotchatt
png("EC_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(EC_processRadtags$prop_reads_retained_per_sample, 
         labels = EC_processRadtags$Filename,
         xlab='EC: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(EC_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("EC_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=EC_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='EC: Proportion of Retained Reads', line=2)
title(ylab='City', line=4.3)

.=dev.off()


## _PA----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(PA_processRadtags$prop_sample_per_library)
sd(PA_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
PA_processRadtags <- PA_processRadtags[order(PA_processRadtags$prop_sample_per_library),]

ppi <- 300
png("PA_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,12,1.5,1.5) + 0.5)
dotchart(PA_processRadtags$prop_sample_per_library,
         labels = PA_processRadtags$Filename,
         xlab='PA: Proportion of Library',
         xlim=c(0,0.85),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(PA_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("PA_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=PA_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.85,
        ylim=c(0,0.9))
title(xlab='PA: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaPA
PA_processRadtags <- PA_processRadtags[order(PA_processRadtags$retained_reads),]

#### Descriptive statistics
summary(PA_processRadtags$retained_reads)
sd(PA_processRadtags$retained_reads)

#### Plot dotchart
png("PA_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(PA_processRadtags$retained_reads/1e6, 
         labels = PA_processRadtags$Filename,
         xlab='PA: Retained Reads (M)',
         xlim=c(0,50),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(PA_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("PA_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=PA_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,50))
title(xlab='PA: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
PA_processRadtags <- PA_processRadtags[order(PA_processRadtags$prop_reads_retained_per_sample),]

#### Descriptive statistics
summary(PA_processRadtags$prop_reads_retained_per_sample)
sd(PA_processRadtags$prop_reads_retained_per_sample)

#### Plot Dotchatt
png("PA_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(PA_processRadtags$prop_reads_retained_per_sample, 
         labels = PA_processRadtags$Filename,
         xlab='PA: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(PA_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("PA_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=PA_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='PA: Proportion of Retained Reads', line=2)
title(ylab='City', line=4.3)

.=dev.off()


## _LS----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(LS_processRadtags$prop_sample_per_library)
sd(LS_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
LS_processRadtags <- LS_processRadtags[order(LS_processRadtags$prop_sample_per_library),]

ppi <- 300
png("LS_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,12,1.5,1.5) + 0.5)
dotchart(LS_processRadtags$prop_sample_per_library,
         labels = LS_processRadtags$Filename,
         xlab='LS: Proportion of Library',
         xlim=c(0,0.7),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(LS_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("LS_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=LS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,0.9))
title(xlab='LS: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaLS
LS_processRadtags <- LS_processRadtags[order(LS_processRadtags$retained_reads),]

#### Descriptive statistics
summary(LS_processRadtags$retained_reads)
sd(LS_processRadtags$retained_reads)

#### Plot dotchart
png("LS_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(LS_processRadtags$retained_reads/1e6, 
         labels = LS_processRadtags$Filename,
         xlab='LS: Retained Reads (M)',
         xlim=c(0,50),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(LS_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("LS_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=LS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,50))
title(xlab='LS: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
LS_processRadtags <- LS_processRadtags[order(LS_processRadtags$prop_reads_retained_per_sample),]

# ID'd a big outlier- will remove
# LS.LA.MAR.U.1
LS_processRadtags <- dplyr::filter(LS_processRadtags, Filename != 'LS.LA.MAR.U.1')

#### Descriptive statistics
summary(LS_processRadtags$prop_reads_retained_per_sample)
sd(LS_processRadtags$prop_reads_retained_per_sample)

#### Plot Dotchatt
png("LS_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(LS_processRadtags$prop_reads_retained_per_sample, 
         labels = LS_processRadtags$Filename,
         xlab='LS: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(LS_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("LS_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=LS_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='LS: Proportion of Retained Reads', line=2)
title(ylab='City', line=4.3)

.=dev.off()


## _TO----

### __A. Proportion of Library each sample comprises----

#### Descriptive statistics
summary(TO_processRadtags$prop_sample_per_library)
sd(TO_processRadtags$prop_sample_per_library)

#### Sort the dataframe by the percent total
TO_processRadtags <- TO_processRadtags[order(TO_processRadtags$prop_sample_per_library),]

ppi <- 300
png("TO_Prop.Lib.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,12,1.5,1.5) + 0.5)
dotchart(TO_processRadtags$prop_sample_per_library,
         labels = TO_processRadtags$Filename,
         xlab='TO: Proportion of Library',
         xlim=c(0,0.7),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=1.5)
abline(v=mean(TO_processRadtags$prop_sample_per_library), lty=2, col='#b22222')

.=dev.off()

#### Proportion of library each sample comprises: Per-Population Boxplot
png("TO_Prop.Lib_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_sample_per_library~city, data=TO_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,0.9))
title(xlab='TO: Proportion of Library', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __B. Reads per-sample----

###$ Sort by retained reaTO
TO_processRadtags <- TO_processRadtags[order(TO_processRadtags$retained_reads),]

#### Descriptive statistics
summary(TO_processRadtags$retained_reads)
sd(TO_processRadtags$retained_reads)

#### Plot dotchart
png("TO_Reads.per.Samp.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,1,1) + 0.5)
dotchart(TO_processRadtags$retained_reads/1e6, 
         labels = TO_processRadtags$Filename,
         xlab='TO: Retained Reads (M)',
         xlim=c(0,30),
         cex = 0.28, 
         cex.axis=0.1,
         cex.lab=2)
abline(v=mean(TO_processRadtags$retained_reads/1e6), lty=2, col='#b22222')

.=dev.off()

#### Per-population boxplot
png("TO_Reads.Per.Sample_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(retained_reads/1e6~city, data=TO_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0,30))
title(xlab='TO: Retained Reads (M)', line=2)
title(ylab='City', line=4.3)

.=dev.off()

### __C. Percent retained----

#### Sort
TO_processRadtags <- TO_processRadtags[order(TO_processRadtags$prop_reads_retained_per_sample),]

#### Descriptive statistics
summary(TO_processRadtags$prop_reads_retained_per_sample)
sd(TO_processRadtags$prop_reads_retained_per_sample)

# ID'd a big outlier- will remove
# TO.BO.WL1.M.2
TO_processRadtags <- dplyr::filter(TO_processRadtags, Filename != 'TO.BO.WL1.M.2')


#### Plot Dotchatt
png("TO_Perc.Retained.png", width=5*ppi, height=10*ppi, res= 300)

par(mar=c(5,9,0.5,0.5) + 0.5)
dotchart(TO_processRadtags$prop_reads_retained_per_sample, 
         labels = TO_processRadtags$Filename,
         xlab='TO: Proportion of Retained Reads',
         xlim=c(0.99975,1),
         cex = 0.28, 
         cex.axis=0.05,
         cex.lab=2)
abline(v=mean(TO_processRadtags$prop_reads_retained_per_sample), lty=2, col='#b22222')

.=dev.off()

#### Per-population Boxplot
png("TO_Perc.Retained_Pop.png", width=5*ppi, height=5*ppi, res= 300)

par(mar=c(3,5,0.5,0.5) + 0.5)
boxplot(prop_reads_retained_per_sample~city, data=TO_processRadtags,
        las=1, 
        xlab='',
        ylab='',
        col=wes_palette("GrandBudapest1"),
        horizontal = TRUE,
        cex.axis=0.7,
        ylim=c(0.99975,1))
title(xlab='TO: Proportion of Retained Reads', line=2)
title(ylab='Population', line=4.3)

.=dev.off()


# _Outlier procedure----

## Samples were removed from downstream analysis by the following process:
## 1. First, if samples represented less than 1% of the library, they were removed.
##    These samples coincided with low read samples and can cause some population-level
##    issues downstream and so were removed.
## 2. Second, I tested out the threshold number for how many samples were removed when
##    not using samples with less than 1 M and 750,000 reads. 