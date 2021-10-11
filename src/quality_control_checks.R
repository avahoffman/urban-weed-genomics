library(ggplot2)
library(ggthemes)


setwd("/Users/avahoffman/Dropbox/Research/Urban_evo/Macrosystems/urban-weed-genomics/")

df <-
  readxl::read_excel(path = "data/macrosystems_urban_leaf_samples_FINAL.xlsx",
                     sheet = "macrosystems_urban_leaf_samples",
                     na = "NA")


ggplot(data = df, aes(x=qubit_conc_ng_uL_in_storage)) +
  geom_density() + 
  xlim(0,400) +
  xlab("DNA concentration (ng/uL)") +
  theme_bw() +
  scale_colour_tableau()


ggsave(filename = "figures/DNA_concen.jpg", height = 3, width = 3)