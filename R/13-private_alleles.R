calc_private_alleles <- function(in_data) {
  if (names(struct_file_LS)[1] != "ind") {
    stop("First column of Structure data input must be named 'ind' for 'individual'.")
  }
  if (names(struct_file_LS)[2] != "pop") {
    stop("Second column of Structure data input must be named 'pop' for population.")
  }
  
  get_private_alleles <- function(col) {
    alleles <- in_data %>% group_by(pop, !!sym(col)) %>% tally() %>% spread(!!sym(col), n, fill = 0)
    is_private <- colSums(alleles != 0) == 1
    priv_loci <- names(is_private[is_private])
    return(alleles %>% ungroup() %>% select(all_of(priv_loci)) %>% rowSums())
  }
  
  out <- sapply(names(in_data)[names(in_data) != "ind" &
                                 names(in_data) != "pop"], get_private_alleles)
  
  return(out)
}

calc_NON_private_alleles <- function(in_data) {
  if (names(struct_file_LS)[1] != "ind") {
    stop("First column of Structure data input must be named 'ind' for 'individual'.")
  }
  if (names(struct_file_LS)[2] != "pop") {
    stop("Second column of Structure data input must be named 'pop' for population.")
  }
  
  get_non_private_alleles <- function(col) {
    alleles <- in_data %>% group_by(pop, !!sym(col)) %>% tally() %>% spread(!!sym(col), n, fill = 0)
    isnot_private <- colSums(alleles != 0) != 1
    notpriv_loci <- names(isnot_private[isnot_private])
    return(alleles %>% ungroup() %>% select(all_of(notpriv_loci)) %>% rowSums())
  }
  
  out <- sapply(names(in_data)[names(in_data) != "ind" &
                                 names(in_data) != "pop"], get_non_private_alleles)
  
  return(out)
}


get_all_private_alleles <- function() {
  struct_file_CD <-
    read.table("SNP_data/zzz_private_alleles/CD_estimatedgeno.structure",
               header = T)
  PAs_CD <- calc_private_alleles(struct_file_CD)
  
  struct_file_DS <-
    read.table("SNP_data/zzz_private_alleles/DS_estimatedgeno.structure",
               header = T)
  PAs_DS <- calc_private_alleles(struct_file_DS)
  
  struct_file_EC <-
    read.table("SNP_data/zzz_private_alleles/EC_estimatedgeno.structure",
               header = T)
  PAs_EC <- calc_private_alleles(struct_file_EC)
  
  struct_file_LS <-
    read.table("SNP_data/zzz_private_alleles/LS_estimatedgeno.structure",
               header = T)
  PAs_LS <- calc_private_alleles(struct_file_LS)
  
  struct_file_PA <-
    read.table("SNP_data/zzz_private_alleles/PA_estimatedgeno.structure",
               header = T)
  PAs_PA <- calc_private_alleles(struct_file_PA)
  
  struct_file_TO <-
    read.table("SNP_data/zzz_private_alleles/TO_estimatedgeno.structure",
               header = T)
  PAs_TO <- calc_private_alleles(struct_file_TO)
  
  par(mfrow = c(3, 2))
  image(PAs_CD)
  image(PAs_DS)
  image(PAs_EC)
  image(PAs_LS)
  image(PAs_PA)
  image(PAs_TO)
  
  NPAs_CD <- calc_NON_private_alleles(struct_file_CD)
  NPAs_DS <- calc_NON_private_alleles(struct_file_DS)
  NPAs_EC <- calc_NON_private_alleles(struct_file_EC)
  NPAs_LS <- calc_NON_private_alleles(struct_file_LS)
  NPAs_PA <- calc_NON_private_alleles(struct_file_PA)
  NPAs_TO <- calc_NON_private_alleles(struct_file_TO)
  
  PAs_all <-
    rbind(
      tibble(
        spp = "CD",
        city = c("BA", "LA", "PX"),
        PA = PAs_CD %>% rowSums(),
        NPA = NPAs_CD %>% rowSums()
      ),
      tibble(
        spp = "DS",
        city = c("BA", "BO", "MN", "PX"),
        PA = PAs_DS %>% rowSums(),
        NPA = NPAs_DS %>% rowSums()
      ),
      tibble(
        spp = "EC",
        city = c("BA", "LA", "PX"),
        PA = PAs_EC %>% rowSums(),
        NPA = NPAs_EC %>% rowSums()
      ),
      tibble(
        spp = "LS",
        city = c("BA", "BO", "LA", "MN", "PX"),
        PA = PAs_LS %>% rowSums(),
        NPA = NPAs_LS %>% rowSums()
      ),
      tibble(
        spp = "PA",
        city = c("BA", "BO", "LA", "MN", "PX"),
        PA = PAs_PA %>% rowSums(),
        NPA = NPAs_PA %>% rowSums()
      ),
      tibble(
        spp = "TO",
        city = c("BA", "BO", "LA", "MN", "PX"),
        PA = PAs_TO %>% rowSums(),
        NPA = NPAs_TO %>% rowSums()
      )
    )
  
  PAs_all <- PAs_all %>% mutate(Percent_private_alleles = PA / (PA + NPA) *
                                  100)
  write_csv(PAs_all, "output/population_stats/private_alleles.csv")
}
