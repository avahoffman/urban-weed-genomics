
spp_ <- "LS"

dat_ <- readr::read_delim(paste0("/Users/avahoffman/Dropbox/Research/Urban_evo/Macrosystems/urban-weed-genomics/output/AMOVA/", spp_, "_estimatedgeno.structure"), delim = "\t")
dat_ <-
  dat_ %>% 
  tidyr::separate(col = 1, into = c("spp", "city", "site", "mgmt", "rep"), "\\.", remove = F) %>% 
  tidyr::unite("site", c(city, site, mgmt), remove = F) %>% 
  dplyr::mutate("site" = as.numeric(factor(site)))

# Look at combos
dat_ %>% select(city, site) %>% distinct() %>% arrange(city, site) %>% print(n = 100)

dat_ <- dat_ %>% 
  select(-c(spp, city, mgmt, rep, `...2`))

write.table(
  dat_,
  sep = "\t",
  row.names = F,
  quote = F,
  file = paste0("/Users/avahoffman/Dropbox/Research/Urban_evo/Macrosystems/urban-weed-genomics/output/AMOVA/", spp_, "_estimatedgeno_bysite.structure")
)

dat