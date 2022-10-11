# This script performs Archetypal analysis
# See Gimbernat-Mayol 2022: https://doi.org/10.1371/journal.pcbi.1010301
# Requires the archetypal-analysis package from python be installed to run.


run_archetype_analysis <- function(
    spp = c("CD", "DS", "EC", "LS", "PA", "TO"),
    city = NULL
) {
  for (spp_ in spp) {
    primary_dir <- here::here()
    
    # Set the working directory so that archetypal analysis can be easily run
    # Dir names are different depending on if the catalog was created for the
    # metapopulation or by speces-city combination.
    if (is.null(city)) {
      temp_dir <-
        paste0(primary_dir,
               "/data/SNP_data/by-species/",
               spp_,
               "/populations_20_pct")
      setwd(temp_dir)
    } else {
      temp_dir <-
        paste0(
          primary_dir,
          "/data/SNP_data/by-city/",
          spp_,
          "_",
          city,
          "/populations_20_pct"
        )
      setwd(temp_dir)
    }
    
    # Run archetypal analysis
    system2('archetypal-analysis',
            '-i populations.snps.vcf -o AA -k 2')
    system2('archetypal-analysis',
            '-i populations.snps.vcf -o AA -k 3')
    system2('archetypal-analysis',
            '-i populations.snps.vcf -o AA -k 4')
    system2('archetypal-analysis',
            '-i populations.snps.vcf -o AA -k 5')
  }
}
