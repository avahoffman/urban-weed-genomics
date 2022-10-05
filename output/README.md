This document is intended to help understand the function of the files contained within the output directory 
of urban-weed-genomics. This folder contains all related output files for the urban-weed-genomics project. 
File descriptions are found below. 

#### `clone_filter_out.csv`

This file contains output from [Step 3 - Remove PCR Clones](/bash/README.md/#step-3---remove-pcr-clones).
It includes the sublibrary name (i.e., file), the number of input reads, the number of output reads that
remain after PCR replicates were removed, and the number of PCR replicates that were identified. 

#### `cstacks-metapop-catalog_samples-included.csv`

This file contains the list of samples that were used to create the Metapopulation catalog for each species.
A subset of samples were chosen randomly from each population (i.e., city) and based on samples that had a
similiar percent coverage and number of retained reads.

#### ` parameter_opt-data.csv`

This file contains the results of [Step 5a - Run `denovo_map.sh`](/bash/README.md/#step-5a---run-denovo_mapsh);
where iterations were run on a [subset of samples](/output/parameter_opt-samples_used_in_popmap.csv) 
to decide an optimal set of parameters to use throughout the stacks pipeline. 
Results from these iterations include the:

1. `no.SNPs` = the number of SNPs identified using each iteration, and
2. `no.r60.loci` = the number of loci identified and shared across 60% of the samples. 

The change in R60 loci, as depicted in this [figure](/figures/parameter_opt_Figures.pdf) 
was used to decide the [optimal set of parameters for each species](/output/parameter_opt-final_params_per_species.csv).

#### `parameter_opt-final_params_per_species.csv`

This file contains the optimal parameters used for each species to create catalogs and stacks.

#### `parameter_opt-samples_used_in_popmap.csv`

The subset of samples for each species that were used to identify optimal parameters are listed here.

#### `process_radtags-discarded_samples.csv`

Samples identified and discarded as part of [Step 4d - Identify low-coverage and low-quality samples from](#step-4d---identify-low-coverage-and-low-quality-samples-from)
are listed here.

#### `process_radtags-kept_samples.csv`

Samples kept for downstream analysis after [Step 4d - Identify low-coverage and low-quality samples from](#step-4d---identify-low-coverage-and-low-quality-samples-from)
are listed here.

#### `process_radtags-library_output.csv`

This files contains the library-wide statistics and output from [Step 4c - Assess the raw, processed, and cleaned data](#step-4c---assess-the-raw-processed-and-cleaned-data).

#### `process_radtags-sample_output.csv`

This files contains the per-sample statistics and output from [Step 4c - Assess the raw, processed, and cleaned data](#step-4c---assess-the-raw-processed-and-cleaned-data).

#### `ustacks-discarded_samples.csv`

Samples identified and discarded as part of [Step 5b - Run `ustacks`](#step-5b---run-ustacks) are listed here.

#### `ustacks-kept_samples.csv`

Samples kept for downstream analysis after [Step 5b - Run `ustacks`](#step-5b---run-ustacks) are listed here.

