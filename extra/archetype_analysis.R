#' This script performs Archetypal analysis
#' See Gimbernat-Mayol 2022: https://doi.org/10.1371/journal.pcbi.1010301
#' Requires the archetypal-analysis package from python be installed to run.
#'
#' Clone from github here:
#' https://github.com/AI-sandbox/archetypal-analysis


run_archetype_analysis <-
  function(spp = c("CD", "DS", "EC", "LS", "PA", "TO"),
           k_vals = 3) {
    for (spp_ in spp) {
      primary_dir <- here::here()
      
      # Set the working directory so that archetypal analysis can be easily run
      temp_dir <-
        paste0(primary_dir,
               "/SNP_data/",
               spp_,
               "/populations_20_pct")
      setwd(temp_dir)
      
      # Run archetypal analysis
      for (k_ in k_vals) {
        system2('archetypal-analysis',
                paste0('-i populations.snps.vcf -o AA -k ', k_))
      }
      
      setwd(here::here())
    }
  }
