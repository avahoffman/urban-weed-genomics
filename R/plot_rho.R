# library(adegenet)
# library(hierfstat)
# library(poppr)
library(reshape2)

read_rhos_dat <- function() {
  CD <- 
    read.csv("output/population_stats/genodive_output_rho_CD.csv")
  DS <-
    read.csv("output/population_stats/genodive_output_rho_DS.csv")
  EC <-
    read.csv("output/population_stats/genodive_output_rho_EC.csv")
  LS <-
    read.csv("output/population_stats/genodive_output_rho_LS.csv")
  PA <-
    read.csv("output/population_stats/genodive_output_rho_PA.csv")
  TO <-
    read.csv("output/population_stats/genodive_output_rho_TO.csv")
  
  rhos_ <- list(
    "CD" = CD,
    "DS" = DS,
    "EC" = EC,
    "LS" = LS,
    "PA" = PA,
    "TO" = TO
  )
  
  return(rhos_)
}


rhos <- read_rhos_dat()

rho_i <- rhos$EC
rownames(rho_i) <- rho_i$X 
rho_i <- rho_i[,-1]
rho_i[upper.tri(rho_i)] <- NA
rho_i[rho_i == 0] <- NA
rho_i$v1 <- rownames(rho_i)
rho_i <- melt(rho_i)

ggplot(rho_i, aes(v1, variable)) +
  theme_minimal() +
  geom_tile(aes(fill = value), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', limit = c(0,0.5)) +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'))