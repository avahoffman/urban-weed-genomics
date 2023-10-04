library(readr)

calc_Pr_structure<-function(K) {
  # see https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf
  # K <- c(-4356, -3983, -3982, -3983, -4006)
  # calc_Pr_structure(K)
  # Should yield
  # [1] 2.159605e-163  2.119416e-01  5.761169e-01  2.119416e-01  2.174919e-11
  K_adj <- sapply(K, function(x) exp(x - mean(K)))
  Pr <- sapply(K, function(x) exp(x- mean(K)) / (sum(K_adj)))
  return(Pr)
}


Pr_plot <- function(spp_, bayes = F) {
  values <- read_csv("output/structure/structure_k_Pr.csv")
  values_filt <- values[(values$species == spp_),]
  if(bayes) {
    values_filt[,3] <- calc_Pr_structure(unlist(values_filt[,3]))
  } else {
    values_filt[,3] <- unlist(values_filt[,3])
  }
  if(spp_ == "LS"){
    # VERY low value for K = 5
    values_filt <- values_filt[1:4,]
  }
  ggplot(values_filt, aes(x = K, y = `lnPr(X|K)`)) +
    geom_point() +
    labs(y = if(bayes) "Pr(K|X)" else "lnPr(X|K)") + 
    geom_line(lty = 3) +
    theme_classic()
}

