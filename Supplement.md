- [About](#about)
- [Introduction](#introduction)
  - [0.1 Raw Data File Naming -
    Sublibraries](#11-raw-data-file-naming---sublibraries)
  - [0.2 A Note on File Transfers](#12-a-note-on-file-transfers)
  - [0.3 A Note on Species Names](#a-note-on-species-names)
- [1 File Transfer](#21-transfer-files)
- [2 File Concatenation and Stacks
  Installation](#22-concatenate-files-and-install-stacks)
  - [2.1 Concatenate Files for each
    Sublibrary](#221-concatenate-files-for-each-sublibrary)
  - [2.2 Download and Install Stacks](#222-download-and-install-stacks)
- [3 PCR Clone Removal](#23-remove-pcr-clones)
  - [3.1 Run PCR Clone Removal
    Script](#231-run-pcr-clone-removal-script)
  - [3.2 Parse PCR Clone Removal
    Results](#232-parse-pcr-clone-removal-results)
- [4 Sample Demultiplexing and
  Filtering](#step-4---demultiplexing-and-sample-filtering)
  - [4.1 Demultiplex and Filter](#step-4a---demultiplex-and-filter)
  - [4.2 Organize files](#step-4b---organize-files)
  - [4.3 Assess the raw, processed, and cleaned
    data](#step-4c---assess-the-raw-processed-and-cleaned-data)
  - [4.4 Identify low-coverage and low-quality samples
    from](#step-4d---identify-low-coverage-and-low-quality-samples-from)
- [5 Stacks: Metapopulation Catalog Building and Parameter
  Search](#step-5---metapopulation-catalog-building-and-parameter-search)
  - [5.1 Run `denovo_map.sh`](#run-denovo_map.sh)
  - [5.2 Run `ustacks`](#step-5b---run-ustacks)
  - [5.3 Correct File Names](#step-5c---correct-file-names)
  - [5.4 Choose catalog
    samples/files](#step-5d---choose-catalog-samplesfiles)
- [6 Stacks: Metapopulation catalog with
  `cstacks`](#step-6---metapopulation-catalog-with-cstacks)
- [7 Stacks: Metapopulation locus matching with
  `sstacks`](#step-7---metapopulation-locus-matching-with-sstacks)
- [8 Genotype probabilities with `polyRAD`](#step-8---polyrad)
  - [8.1 Make `RADdata` object](#make-raddata-object)
  - [8.2 Calculate overdispersion](#calculate-overdispersion)
  - [8.3 Estimate genotypes](#estimate-genotypes)
  - [8.4 Final filter and file cleanup](#final-filter-and-file-cleanup)
- [9 `Structure` Analysis](#structure-analysis)
  - [9.1 Running Structure on `polyRAD`
    output](#running-structure-on-polyrad-output)
- [10 Conceptual Figure](#conceptual-figure)
- [11 NLCD Data and Site Plots](#nlcd-data-and-site-plots)
  - [11.1 Preparing NLCD Data](#preparing-nlcd-data)
  - [11.2 Climate normals](#climate-normals)
  - [11.3 Maps of sampling locations](#maps-of-sampling-locations)
- [12 Principal components analysis &
  plots](#principal-components-analysis-plots)
- [13 Genetic structure using Structure and
  sNMF](#genetic-structure-using-structure-and-snmf)
  - [13.1 `Structure` optimal K](#structure-optimal-k)
  - [13.2 Plotting Structure output](#plotting-structure-output)
  - [13.3 Validation of Structure results with
    sNMF](#validation-of-structure-results-with-snmf)
- [14 Genetic Differentiation](#genetic-differentiation)
  - [14.1 Among city - Jost’s D, $G_{ST}$,
    $F_{ST}$](#among-city---josts-d-g_st-f_st)
  - [14.2 Among city - pairwise $\rho$](#among-city---pairwise-rho)
  - [14.3 Within city - allelic
    richness](#within-city---allelic-richness)
  - [14.4 Within city - $F_{IS}$ (homozygosity within
    population)](#within-city---f_is-homozygosity-within-population)
  - [14.5 Within city - $\bar{r}_d$ - Linkage
    disequilibrium](#within-city---barr_d---linkage-disequilibrium)
  - [14.6 Within city - private alleles](#within-city---private-alleles)
  - [14.7 AMOVA](#amova)
- [15 `SessionInfo()`](#sessioninfo)
- [Appendix](#appendix-books)
  - [15.1 File Organization](#file-organization-bookmark_tabs)
  - [15.2 Aspera Transfer File Names](#aspera-transfer-file-names)
  - [15.3 `clone_filter` File Names](#clone_filter-file-names)

<!-- To nicely format this markdown doc, I recommend using
Sublime Text editor. Highlight the text in question and select Edit >
Wrap > Wrap paragraph at 72 characters. --> 
<!-- Allow "S" labeling on Figs and Tables  -->
<!-- End -->

# About

This is supplementary material intended to be paired with files on
GitHub at the repository
<https://github.com/avahoffman/urban-weed-genomics>.

# Introduction

In this experiment, we used quaddRAD library prep to prepare the sample
DNA. This means that there were both two unique outer barcodes (typical
Illumina barcodes) *AND* two unique inner barcodes (random barcode bases
inside the adapters) for each sample - over 1700 to be exact!

The sequencing facility demultiplexes samples based on the outer
barcodes (typically called 5nn and i7nn). Once this is done, each file
still contains a mix of the inner barcodes. We will refer to these as
“sublibraries” because they are sort of halfway demultiplexed. We
separate them out bioinformatically later.

## 0.1 Raw Data File Naming - Sublibraries

Here’s a bit of information on the file name convention. The typical raw
file looks like this:

`AMH_macro_1_1_12px_S1_L001_R1_001.fastq.gz`

- These are author initials and “macro” stands for “Macrosystems”. These
  are on every file.

  `AMH_macro`

- The first number is the *i5nn* barcode for the given sublibrary. We
  know all these samples have a i5nn barcode “1”, so that narrows down
  what they can be. The second number is the *i7nn* barcode for the
  given sublibrary. We know all these samples have a i7nn barcode “1”,
  so that further narrows down what they can be.

  `1_1`

- This refers to how many samples are in the sublibrary. “12px” means
  12-plexed, or 12 samples. In other words, we will use the inner
  barcodes to further distinguish 12 unique samples in this sublibrary.

  `12px`

- This is a unique sublibrary name. S1 = 1 i5nn and 1 i7nn.

  `S1`

- This means this particular file came from lane 1 of the NovaSeq. There
  are four lanes. All samples should appear across all four lanes.

  `L001`

- This is the first (R1) of two paired-end reads (R1 and R2).

  `R1`

- The last part doesn’t mean anything - it was just added automatically
  before the file suffix (`fastq.gz`)

  `001.fastq.gz`

## 0.2 A Note on File Transfers

There are three main systems at play for file transfer: the local
machine, the sequencing facility’s (GRCF) Aspera server, and
[MARCC](https://www.marcc.jhu.edu/). The Aspera server is where the data
were/are stored immediately after sequencing. MARCC is where we plan to
do preprocessing and analysis. Scripts and text files are easy for me to
edit on my local machine. We used [Globus](https://app.globus.org) to
transfer these small files from my local machine to MARCC.

Midway through this analyses, we transitioned to another cluster, JHU’s
Rockfish. Scripts below, with the exception of file transfer from the
Aspera server, should reflect the new filesystem, though you will have
to adjust the file paths accordingly.

![File transfer schematic](figures/file_transfer.jpg)  

## 0.3 A Note on Species Names

Throughout this study, we examined 6 species. Sometimes we used
abbreviations for easier file naming. These are:

1.  *Cynodon dactylon*, “CD”, or Bermuda grass

2.  *Digitaria sanguinalis*, “DS”, or crabgrass

3.  *Erigeron canadensis*, “EC”, or horseweed

4.  *Lactuca serriola*, “LS”, or prickly lettuce

5.  *Poa annua*, “PA”, or bluegrass

6.  *Taraxacum officinale*, “TO”, or dandelion

# 1 File Transfer

Referred to through files as “Step 1”. Files can be found in the
`01_transfer_files/` directory.

This directory contains files named in this convention:
`01-aspera_transfer_n.txt`. These are text files containing the *names*
of `fastq.gz` files that we wanted to transfer from the sequencing
facility’s Aspera server to the computing cluster
([MARCC](https://www.marcc.jhu.edu/)). This was to maximize ease of
transferring only certain files over at once, since transferring could
take a long time. We definitely did this piecemeal. Possible file names
shown in [Aspera Transfer File Names](#aspera-transfer-file-names).
There are multiple of these files so that we could parallelize (replace
n with the correct number in the command used below). This text file
will need to be uploaded to your scratch directory in MARCC.

Files were then transferred using the following commands. Before
starting, make sure you are in a data transfer node. Then, load the
aspera module. Alternatively, you can install the Aspera transfer
software and use that.

    module load aspera

Initiate the transfer from within your scratch directory:

    ascp -T -l8G -i /software/apps/aspera/3.9.1/etc/asperaweb_id_dsa.openssh
    --file-list=01-aspera_transfer_n.txt
    --mode=recv --user=<aspera-user> --host=<aspera-IP> /scratch/users/<me>@jhu.edu

# 2 File Concatenation and Stacks Installation

Referred to through files as “Step 2”. Files can be found in the
`02_concatenate_and_check/` directory.

## 2.1 Concatenate Files for each Sublibrary

**Step 2a**. We ran my samples across the whole flow cell of the
NovaSeq, so results came in 8 files for each demultiplexed sublibrary (4
lanes \* paired reads). For example, for sublibrary 1_1, we’d see the
following 8 files:

    AMH_macro_1_1_12px_S1_L001_R1_001.fastq.gz
    AMH_macro_1_1_12px_S1_L001_R2_001.fastq.gz
    AMH_macro_1_1_12px_S1_L002_R1_001.fastq.gz
    AMH_macro_1_1_12px_S1_L002_R2_001.fastq.gz    
    AMH_macro_1_1_12px_S1_L003_R1_001.fastq.gz
    AMH_macro_1_1_12px_S1_L003_R2_001.fastq.gz
    AMH_macro_1_1_12px_S1_L004_R1_001.fastq.gz
    AMH_macro_1_1_12px_S1_L004_R2_001.fastq.gz

The `02_concatendate_and_check/02-concat_files_across4lanes.sh` script
finds all files in the working directory with the name pattern
`*_L001_*.fastq.gz` and then concatenates across lanes 001, 002, 003,
and 004 so they can be managed further. The “L001” part of the filename
is then eliminated. For example the 8 files above would become:

    AMH_macro_1_1_12px_S1_R1.fastq.gz
    AMH_macro_1_1_12px_S1_R2.fastq.gz

Rockfish uses [slurm](https://slurm.schedmd.com/overview.html) to manage
jobs. To run the script, use the `sbatch` command. For example:

    sbatch ~/code/02-concat_files_across4lanes.sh

This command will run the script from within the current directory, but
will look for and pull the script from the code directory. This will
concatenate all files within the current directory that match the loop
pattern.

## 2.2 Download and Install Stacks

**Step 2b**. On Rockfish,
[Stacks](https://catchenlab.life.illinois.edu/stacks/) will need to be
downloaded to each user’s code directory. Stacks, and software in
general, should be compiled in an interactive mode or loaded via module.
For more information on interactive mode, see `interact --usage`.

    interact -p debug -g 1 -n 1 -c 1
    module load gcc

Now download Stacks. We used version 2.60.

    wget http://catchenlab.life.illinois.edu/stacks/source/stacks-2.60.tar.gz
    tar xfvz stacks-2.60.tar.gz

Next, go into the stacks-2.60 directory and run the following commands:

    ./configure --prefix=/home/<your_username>/code4-<PI_username>
    make
    make install
    export PATH=$PATH:/home/<your_username>/code4-<PI_username>/stacks-2.60

The filesystem patterns on your cluster might be different, and you
should change these file paths accordingly.

# 3 PCR Clone Removal

Referred to through files as “Step 3”. Files can be found in the
`03_clone_filter/` directory.

## 3.1 Run PCR Clone Removal Script

**Step 3a**. The `03-clone_filter.sh` script runs `clone_filter` from
[Stacks](https://catchenlab.life.illinois.edu/stacks/). The program was
run with options `--inline_inline --oligo_len_1 4 --oligo_len_2 4`. The
`--oligo_len_x 4` options indicate the 4-base pair degenerate sequence
was included on the outside of the barcodes for detecting PCR
duplicates. The script uses the file name prefixes listed for each
single sub-pooled library in `03-clone_filter_file_names.txt` and loops
to run `clone_filter` on all of them. Possible file names shown in
[`clone_filter` File Names](#clone_filter-file-names).

## 3.2 Parse PCR Clone Removal Results

**Step 3b**. If you want to extract descriptive statistics from the
`clone_filter` output, you can use the `03.5-parse_clone_filter.py`
script to do so. It can be run on your local terminal after transferring
the `clone_filter.out` logs to your local computer.

``` r
source("03_clone_filter/examine_clones.R")
make_cloneplot()
```

<figure>
<img src="Supplement_files/figure-gfm/PCR_clones-1.png"
alt="PCR clone removal statistics" />
<figcaption aria-hidden="true">PCR clone removal statistics</figcaption>
</figure>

# 4 Sample Demultiplexing and Filtering

Files can be found in the `04_demux_filter/` directory.

## 4.1 Demultiplex and Filter

*Step 4a*. The `04-process_radtags.sh` script runs `process_radtags`
from [Stacks](https://catchenlab.life.illinois.edu/stacks/). The program
was run with options
`-c -q --inline_inline --renz_1 pstI --renz_2 mspI --rescue --disable_rad_check`.
The script uses the same file prefixes as [Step 3 -
`03-clone_filter.sh`](#step-3---remove-pcr-clones). Each sub-pooled
library has a forward and reverse read file that was filtered in the
previous step. Like the [above section](#step-3---03-clone_filtersh),
the script uses the file name prefixes listed for each single sub-pooled
library in `04-process_radtags_file_names.txt` and loops to run
`process_radtags` on all of them. Possible file names shown in
[`clone_filter` File Names](#clone_filter-file-names).

Each sub-pooled library also has a demultiplexing file (`04-demux/`
directory) that contains the sample names and inner(i5 and i7) barcodes.
For example, the sublibrary 1_1, we’d see the following barcode file:

    ATCACG  AGTCAA  DS.BA.PIK.U.1
    CGATGT  AGTTCC  DS.BA.PIK.U.2
    TTAGGC  ATGTCA  DS.BA.PIK.U.3
    TGACCA  CCGTCC  DS.BA.PIK.U.4
    ACAGTG  GTCCGC  DS.BA.PIK.U.5
    GCCAAT  GTGAAA  DS.BA.DHI.U.1
    CAGATC  GTGGCC  DS.BA.DHI.U.2
    ACTTGA  GTTTCG  DS.BA.DHI.U.3
    GATCAG  CGTACG  DS.BA.DHI.U.4
    TAGCTT  GAGTGG  DS.BA.DHI.U.5
    GGCTAC  ACTGAT  DS.BA.GA.U.1
    CTTGTA  ATTCCT  DS.BA.GA.U.2

The ‘process_radtags’ command will demultiplex the data by separating
out each sublibrary into the individual samples. It will then clean the
data, and will remove low quality reads and discard reads where a
barcode was not found.

## 4.2 Organize files

*Step 4b*. In a new directory, make sure the files are organized by
species. In the `process_radtags` script, we specified that files be
sent to `~/scratch/demux/*sublibrary_name*` (reasoning for this is in
[Step 4c](#step-4c---assess-the-raw-processed-and-cleaned-data)), but
files should manually be organized into species folders (i.e.,
`~/scratch/demux/*SPP*`) after `process_radtags` is performed. For
example, the file “DS.MN.L01-DS.M.1.1.fq.gz” should be sent to the
`~/scratch/demux/DS` directory.

Note: this is not automated at this point but it would be nice to
automate the file moving process so it’s not forgotten at this point.

## 4.3 Assess the raw, processed, and cleaned data

*Step 4c*. In the script for [Step 4](#step-4---demultiplex-and-filter),
we have specified that a new output folder be created for each
sublibrary. The output folder is where all sample files and the log file
will be dumped for each sublibrary. It is important to specify a
different output folder if you have multiple sublibraries because we
will be assessing the output log for each sublibrary individually (and
otherwise, the log is overwritten when the script loops to a new
sublibrary).

The utility `stacks-dist-extract` can be used to extract data from the
log file. First, we examined the library-wide statistics to identify
sublibraries where barcodes may have been misentered or where sequencing
error may have occurred. We used:

     stacks-dist-extract process_radtags.log total_raw_read_counts
     

to pull out data on the total number of sequences, the number of
low-quality reads, whether barcodes were found or not, and the total
number of retained reads per sublibary. Look over these to make sure
there are no outliers or sublibraries that need to be checked and rerun.

Next, we used:

     stacks-dist-extract process_radtags.log per_barcode_raw_read_counts
     

to analyze how well each sample performed. There are three important
statistics to consider for each sample.

1.  *The proportion of reads per sample for each sublibrary* indicates
    the proportion that each individual was processed and sequenced
    within the overall library. This is important to consider as cases
    where a single sample dominates the sublibrary may indicate
    contamination. (Field `prop_sample_per_library`).

2.  *The number of reads retained for each sample* can be an indicator
    of coverage. It is most likely a good idea to remove samples with a
    very low number of reads. Where you decide to place the cutoff for
    low coverage samples is dependent on your dataset. For example, a
    threshold of 1 million reads is often used but this is not
    universal. (Field `retained_reads`).

3.  *The proportion of reads retained for each sample* can also indicate
    low-quality samples and will give an idea of the variation in
    coverage across samples. (Field `prop_reads_retained_per_sample`).

Output for sublibraries for this step are summarized in
[`process_radtags-library_output.csv`](output/process_radtags-sample_output.csv).

Output for individual samples for this step are summarized in
[`process_radtags-sample_output.csv`](output/process_radtags-sample_output.csv).

The script `04c-process_radtags_stats.R` was used to create many plots
for easily assessing each statistic. Output from this step can be found
in `figures/process_radtags/` where figures are organized by species.

The script `04c-radtags_filter_summary.R` summarizes the filtering
results from all samples.

``` r
source("04_demux_filter/04c-radtags_filter_summary.R")
make_filterplot()
```

<figure>
<img src="Supplement_files/figure-gfm/radtag_stats-1.png"
alt="RAD tag processing statistics" />
<figcaption aria-hidden="true">RAD tag processing
statistics</figcaption>
</figure>

## 4.4 Identify low-coverage and low-quality samples from

*Step 4d*. Using the same output log and the above statistics, we
removed low-coverage and low-quality samples that may skew downstream
analyses.

Samples were identified and removed via the following procedure:

1.  First, samples that represented less than **1% of the sequenced
    sublibrary** were identified and removed. These samples correlate to
    low-read and low-coverage samples.

2.  Next, a threshold of **1 million retained reads per sample** was
    used to remove any remaining low-read samples. Low-read samples
    correlate to low coverage and will lack enough raw reads to
    contribute to downstream analyses.

Good/kept samples are listed in
[`process_radtags-kept_samples.csv`](output/process_radtags-kept_samples.csv).

Discarded samples are listed in
[`process_radtags-discarded_samples.csv`](output/process_radtags-discarded_samples.csv).

``` r
source("04_demux_filter/04c-radtags_filter_summary.R")
make_manual_discard_plot()
```

<figure>
<img src="Supplement_files/figure-gfm/radtag_manual_filter-1.png"
alt="RAD tag manual filtering summary" />
<figcaption aria-hidden="true">RAD tag manual filtering
summary</figcaption>
</figure>

Note: At this point, we started using Stacks 2.62 for its
multi-threading capabilities. Functionality of the previous steps should
be the same, however.

# 5 Stacks: Metapopulation Catalog Building and Parameter Search

Files can be found in the `05_ustacks_and_params/` directory.

Going forward, when we use the term **metapopulation**, we are referring
to the collection of all samples within species among all cities where
the species was present.

It is important to conduct preliminary analyses that will identify an
optimal set of parameters for the dataset (see [Step
5a](#step-5a---denovo_mapsh)). Following the parameter optimization, the
program `ustacks` can be run to generate a catalog of loci.

## 5.1 Run `denovo_map.sh`

*Step 5a*. Stack assembly will differ based on several different aspects
of the dataset(such as the study species, the RAD-seq method used,
and/or the quality and quantity of DNA used). So it is important to use
parameters that will maximize the amount of biological data obtained
from stacks.

There are three main parameters to consider when doing this:

1.  *m* = controls the minimum number of raw reads required to form a
    stack(implemented in `ustacks`)

2.  *M* = controls the number of mismatches between stacks to to merge
    them into a putative locus (implemented in `ustacks`)

3.  *n* = controls the number of mismatches allowed between stacks to
    merge into the catalog (implemented in `cstacks`)

There are two main ways to optimize parameterization:

1.  an iterative method were you sequentially change each parameter
    while keeping the other parameters fixed (described in *Paris et al.
    2017*), or

2.  an iterative method were you sequentially change the values of *M*
    and *n*(keeping *M* = *n*) while fixing *m* = 3, and then test *m* =
    2, 4 once the optimal *M* = *n* is determined(described in *Rochette
    and Catchen 2017*, *Catchen 2020*).

We performed the second method and used the `denovo_map.sh` script to
run the `denovo_map.pl` command to perform iterations. This script
requires that we first choose a subset of samples to run the iterations
on. The samples should be representative of the overall dataset; meaning
they should include all populations and have similar read coverage
numbers. Read coverage numbers can be assessed by looking at the
descriptive statistics produced from [Step
4c](#step-4c---assess-the-raw-processed-and-cleaned-data).

Place these samples in a text file (`popmap_test_samples.txt`) with the
name of the sample and specify that all samples belong to the same
population. For example, `popmap_test_samples.txt` should look like…

    DS.BA.GA.U.1    A
    DS.PX.BUF.M.5   A
    DS.B0.HC4.M.1   A
    ...

It is important to have all representative samples treated as one
population because you will assess outputs found across 80% of the
individuals. The script will read this text file from the `--popmap`
argument.

The script also requires that you specify an output directory after
`-o`. This should be unique to the parameter you are testing… for
example, if you are testing *M* = 3, then you could make a subdirectory
labeled `stacks.M3` where all outputs from `denovo_map.sh` will be
placed. Otherwise, for each iteration, the outputs will be overwritten
and you will lose the log from the previous iteration. The
`denovo_map.sh` script also requires that you direct it toward where
your samples are stored, which is your directory built in [Step
4b](#step-4b---organize-files). Make sure to run the
`--min-samples-per-pop 0.80` argument.

To decide which parameters to use, examine the following from each
iteration:

1.  the average sample coverage: This is obtained from the summary log
    in the `ustacks` section of `denovo_map.log`. If samples have a
    coverage \<10x, you will have to rethink the parameters you use
    here.

2.  the number of assembled loci shared by 80% of samples: This can be
    found in the `haplotypes.tsv` by counting the number of loci:
    `cat populations.haplotypes.tsv | grep -v ^"#" | wc -l`

3.  the number of polymorphic loci shared by 80% of samples: This can be
    found in `populations.sumstats.tsv` or by counting
    `populations.hapstats.tsv`:
    `cat populations.hapstats.tsv | grep -v "^#" | wc -l`

4.  the number of SNPs per locus shared by 80% of samples: found in
    `denovo_map.log` or by counting the number of SNPs in
    `populations.sumstats.tsv`:
    `populations.sumstats.tsv | grep -v ^"#" | wc -l`

The script `05a-param_opt-figures_script.R` was used to create plots for
assessing the change in shared loci across parameter iterations.

Based on this optimization step, we used the following parameters:

| Species | M (locus mismatches) | n (catalog mismatches) | m (minimum reads) |
|:--------|---------------------:|-----------------------:|------------------:|
| CD      |                    8 |                      8 |                 3 |
| DS      |                   10 |                     10 |                 3 |
| EC      |                    8 |                      8 |                 3 |
| LS      |                    7 |                      7 |                 3 |
| PA      |                    5 |                      5 |                 3 |
| TO      |                    6 |                      6 |                 3 |

Final parameter optimization values for the Stacks pipeline.

## 5.2 Run `ustacks`

*Step 5b*. `ustacks` builds *de novo* loci in each individual sample. We
have designed the `ustacks` script so that the process requires three
files:

- `05-ustacks_n.sh` : the shell script that executes `ustacks`
- `05-ustacks_id_n.txt` : the sample ID number
- `05-ustacks_samples_n.txt` : the sample names that correspond to the
  sample IDs

The sample ID should be derived from the `order_id` column(first column)
on the master spreadsheet. It is unique (1-1736) across all of the
samples.

The sample name is the corresponding name for each sample ID in the
spreadsheet. E.g., sample ID “9” corresponds to sample name
“DS.BA.DHI.U.4”. Sample naming convention is
species.city.site.management_type.replicate_plant.

`05-ustacks_n.sh` should have an out_directory (`-o` option) that will
be used for all samples (e.g., `stacks/ustacks`). Files can be processed
piecemeal into this directory. There should be three files for every
sample in the output directory:

- `<samplename>.alleles.tsv.gz`
- `<samplename>.snps.tsv.gz`
- `<samplename>.tags.tsv.gz`

Multiple versions of the `05-ustacks_n.sh` script can be run in parallel
(simply replace n in the three files above with the correct number).

A small number of samples (13) were discarded at this stage as the
`ustacks` tool was unable to form any primary stacks corresponding to
loci. See
[output/ustacks-discarded_samples.csv](output/ustacks-discarded_samples.csv).

| ustacks discarded | Retained reads | Proportion of sub-library |
|:------------------|---------------:|--------------------------:|
| no                |       10075192 |                 0.1608982 |
| yes               |        7490510 |                 0.1230658 |

Summary of samples discarded at the ustacks step of the Stacks pipeline.

## 5.3 Correct File Names

*Step 5c*. This step contains a script `05b-fix_filenames.sh` which uses
some simple regex to fix filenames that are output in previous steps.
Stacks adds an extra “1” at some point at the end of the sample name
which is not meaningful. The following files:

- DS.MN.L02-DS.M.3.1.alleles.tsv.gz
- DS.MN.L03-DS.U.2.1.tags.tsv.gz
- DS.MN.L09-DS.U.1.1.snps.tsv.gz

become:

- DS.MN.L02-DS.M.3.alleles.tsv.gz
- DS.MN.L03-DS.U.2.tags.tsv.gz
- DS.MN.L09-DS.U.1.snps.tsv.gz

The script currently gives some strange log output, so it can probably
be optimized/improved. The script should be run from the directory where
the changes need to be made. Files that have already been fixed will not
be changed.

## 5.4 Choose catalog samples/files

*Step 5d*. In the next step, we will choose the files we want to go into
the catalog. This involves a few steps:

1.  Create a meaningful directory name. This could be the date (e.g.,
    `stacks_22_01_25`).

2.  Copy the `ustacks` output for all of the files you want to use in
    the reference from Step 5b. Remember this includes three files per
    sample. So if you have 20 samples you want to include in the
    reference catalog, you will transfer 3 x 20 = 60 files into the
    meaningful directory name. The three files per sample should follow
    this convention:

- `<samplename>.alleles.tsv.gz`
- `<samplename>.snps.tsv.gz`
- `<samplename>.tags.tsv.gz`

3.  Remember the meaningful directory name. You will need it in Step 6.

# 6 Stacks: Metapopulation catalog with `cstacks`

Files can be found in the `06_cstacks/` directory.

`cstacks` builds the locus catalog from all the samples specified. The
accompanying script, `cstacks_SPECIES.sh` is relatively simple since it
points to the directory containing all the sample files. It follows this
format to point to that directory:

    cstacks -P ~/directory ...

Make sure that you use the meaningful directory from Step 5c and that
you have copied all the relevant files over. Otherwise this causes
[problems
downstream](https://groups.google.com/g/stacks-users/c/q3YYPprmYnU/m/cH5RB5KwBQAJ).
For example, you might edit the code to point to
`~/scratch/stacks/stacks_22_01_25`.

    cstacks -P ~/scratch/stacks/stacks_22_01_25 ...

The tricky thing is ensuring enough compute memory to run the entire
process successfully. There is probably space to optimize this process.

The `cstacks` method uses a “population map” file, which in this project
is `cstacks_popmap_SPECIES.txt`. This file specifies which samples to
build the catalog from and categorizes them into your ‘populations’, or
in this case, cities using two tab-delimited columns, e.g.:

    DS.BA.GA.U.1    Baltimore
    DS.BA.GA.U.2    Baltimore
    DS.BA.GA.U.3    Baltimore
    DS.BA.GA.U.4    Baltimore
    DS.BA.GA.U.5    Baltimore
    ...

Make sure the samples in this file correspond to the input files located
in e.g., `~/scratch/stacks/stacks_22_01_25`.

`cstacks` builds three files for use in all your samples (in this
pipeline run), mirroring the sample files output
by[`ustacks`](#step-5---ustacks):

- `catalog.alleles.tsv.gz`
- `catalog.snps.tsv.gz`
- `catalog.tags.tsv.gz`

| Sample           | Species | City |
|:-----------------|:--------|:-----|
| DS.BA.PIK.U.1    | DS      | BA   |
| DS.BA.GA.U.4     | DS      | BA   |
| DS.BA.LH-1.M.4   | DS      | BA   |
| DS.BA.LH-3.M.1   | DS      | BA   |
| DS.BA.WB.U.2     | DS      | BA   |
| DS.BA.LL-4.M.5   | DS      | BA   |
| DS.BA.LH-2.M.5   | DS      | BA   |
| DS.BA.TRC.U.3    | DS      | BA   |
| DS.BA.W3.M.2     | DS      | BA   |
| DS.BA.RG-1.M.1   | DS      | BA   |
| DS.BA.LL-3.M.3   | DS      | BA   |
| DS.BA.RG-2.M.4   | DS      | BA   |
| DS.BO.HC1.M.3    | DS      | BO   |
| DS.BO.HC4.M.5    | DS      | BO   |
| DS.BO.LC1.M.3    | DS      | BO   |
| DS.BO.LC2.M.2    | DS      | BO   |
| DS.BO.LC3.M.5    | DS      | BO   |
| DS.BO.WL1.M.2    | DS      | BO   |
| DS.BO.WL2.M.1    | DS      | BO   |
| DS.BO.WL3.M.5    | DS      | BO   |
| DS.BO.I4.U.1     | DS      | BO   |
| DS.BO.R1.U.4     | DS      | BO   |
| DS.BO.R2.U.2     | DS      | BO   |
| DS.BO.R4.U.4     | DS      | BO   |
| DS.MN.L05-DS.M.3 | DS      | MN   |
| DS.MN.L09-DS.M.3 | DS      | MN   |
| DS.MN.L11-DS.M.1 | DS      | MN   |
| DS.MN.L02-DS.U.1 | DS      | MN   |
| DS.MN.L02-DS.M.4 | DS      | MN   |
| DS.MN.L03-DS.U.3 | DS      | MN   |
| DS.MN.L04-DS.U.5 | DS      | MN   |
| DS.MN.L06-DS.U.3 | DS      | MN   |
| DS.MN.L07-DS.U.3 | DS      | MN   |
| DS.MN.L09-DS.U.3 | DS      | MN   |
| DS.MN.L11-DS.U.1 | DS      | MN   |
| DS.MN.L11-DS.U.5 | DS      | MN   |
| DS.PX.BUF.M.1    | DS      | PX   |
| DS.PX.PIE.M.2    | DS      | PX   |
| DS.PX.ALA.M.1    | DS      | PX   |
| DS.PX.MTN.M.6    | DS      | PX   |
| DS.PX.LAP.M.3    | DS      | PX   |
| DS.PX.NUE.M.4    | DS      | PX   |
| DS.PX.WES.M.2    | DS      | PX   |
| DS.PX.DF1.M.1    | DS      | PX   |
| DS.PX.ENC.M.1    | DS      | PX   |
| DS.PX.DOW.M.1    | DS      | PX   |
| DS.PX.DOW.M.4    | DS      | PX   |
| DS.PX.DF2.M.3    | DS      | PX   |
| CD.BA.LA.U.2     | CD      | BA   |
| CD.BA.TRC.U.3    | CD      | BA   |
| CD.BA.WGP.M.2    | CD      | BA   |
| CD.BA.LH-2.M.2   | CD      | BA   |
| CD.BA.LL-4.M.1   | CD      | BA   |
| CD.BA.PIK.U.2    | CD      | BA   |
| CD.BA.WB.U.2     | CD      | BA   |
| CD.BA.CP.U.4     | CD      | BA   |
| CD.BA.FH.U.1     | CD      | BA   |
| CD.BA.PSP.M.4    | CD      | BA   |
| CD.BA.AA.U.4     | CD      | BA   |
| CD.BA.RG-1.M.2   | CD      | BA   |
| CD.BA.W3.M.3     | CD      | BA   |
| CD.BA.GA.U.3     | CD      | BA   |
| CD.BA.WBO.U.5    | CD      | BA   |
| CD.LA.WHI.M.3    | CD      | LA   |
| CD.LA.SEP.M.3    | CD      | LA   |
| CD.LA.SEP.M.4    | CD      | LA   |
| CD.LA.ROS.M.5    | CD      | LA   |
| CD.LA.MR2.M.2    | CD      | LA   |
| CD.LA.ALL.M.2    | CD      | LA   |
| CD.LA.ALL.M.5    | CD      | LA   |
| CD.LA.VAL.M.5    | CD      | LA   |
| CD.LA.HAR.M.4    | CD      | LA   |
| CD.LA.LUB.M.3    | CD      | LA   |
| CD.LA.GLO.M.4    | CD      | LA   |
| CD.LA.ZOO.M.3    | CD      | LA   |
| CD.LA.NWH.M.5    | CD      | LA   |
| CD.LA.KIN.M.3    | CD      | LA   |
| CD.LA.KIN.M.5    | CD      | LA   |
| CD.PX.CAM.U.5    | CD      | PX   |
| CD.PX.MON.U.5    | CD      | PX   |
| CD.PX.PKW.U.5    | CD      | PX   |
| CD.PX.LAP.M.4    | CD      | PX   |
| CD.PX.NES.U.4    | CD      | PX   |
| CD.PX.PAL.M.3    | CD      | PX   |
| CD.PX.ASU.M.1    | CD      | PX   |
| CD.PX.NUE.M.5    | CD      | PX   |
| CD.PX.WES.M.3    | CD      | PX   |
| CD.PX.MAN.M.4    | CD      | PX   |
| CD.PX.CLA.M.3    | CD      | PX   |
| CD.PX.DF1.M.5    | CD      | PX   |
| CD.PX.COY.M.5    | CD      | PX   |
| CD.PX.RPC.M.3    | CD      | PX   |
| CD.PX.ENC.M.2    | CD      | PX   |
| EC.BA.LH-2.M.2   | EC      | BA   |
| EC.BA.WBO.U.4    | EC      | BA   |
| EC.BA.WB.U.5     | EC      | BA   |
| EC.BA.FH.U.3     | EC      | BA   |
| EC.BA.CP.U.2     | EC      | BA   |
| EC.BA.TRC.U.3    | EC      | BA   |
| EC.BA.LL-4.M.4   | EC      | BA   |
| EC.BA.WB.U.1     | EC      | BA   |
| EC.BA.PIK.U.5    | EC      | BA   |
| EC.BA.PSP.M.4    | EC      | BA   |
| EC.BA.GA.U.2     | EC      | BA   |
| EC.BA.LL-3.M.3   | EC      | BA   |
| EC.BA.ML.U.1     | EC      | BA   |
| EC.BA.TRC.U.5    | EC      | BA   |
| EC.BA.ML.U.3     | EC      | BA   |
| EC.LA.SGB.U.2    | EC      | LA   |
| EC.LA.SGB.U.5    | EC      | LA   |
| EC.LA.DUR.U.2    | EC      | LA   |
| EC.LA.HOW.U.2    | EC      | LA   |
| EC.LA.SAN.U.2    | EC      | LA   |
| EC.LA.VER.U.1    | EC      | LA   |
| EC.LA.VER.U.4    | EC      | LA   |
| EC.LA.VB2.U.4    | EC      | LA   |
| EC.LA.AC2.U.2    | EC      | LA   |
| EC.LA.AC1.U.1    | EC      | LA   |
| EC.LA.VB1.U.1    | EC      | LA   |
| EC.LA.VB1.U.3    | EC      | LA   |
| EC.LA.SGR.U.4    | EC      | LA   |
| EC.LA.SGR.U.5    | EC      | LA   |
| EC.LA.HOW.U.3    | EC      | LA   |
| EC.PX.BUF.M.1    | EC      | PX   |
| EC.PX.BUF.M.3    | EC      | PX   |
| EC.PX.ALA.M.3    | EC      | PX   |
| EC.PX.MTN.M.2    | EC      | PX   |
| EC.PX.WES.M.1    | EC      | PX   |
| EC.PX.WES.M.2    | EC      | PX   |
| EC.PX.MAN.M.1    | EC      | PX   |
| EC.PX.CLA.M.1    | EC      | PX   |
| EC.PX.PSC.M.1    | EC      | PX   |
| EC.PX.DF1.M.1    | EC      | PX   |
| EC.PX.DOW.M.1    | EC      | PX   |
| EC.PX.DOW.M.2    | EC      | PX   |
| EC.PX.COY.M.2    | EC      | PX   |
| EC.PX.COY.M.3    | EC      | PX   |
| EC.PX.ALA.M.5    | EC      | PX   |
| LS.BA.WB.U.1     | LS      | BA   |
| LS.BA.WB.U.2     | LS      | BA   |
| LS.BA.DHI.U.2    | LS      | BA   |
| LS.BA.GA.U.1     | LS      | BA   |
| LS.BA.PIK.U.3    | LS      | BA   |
| LS.BA.PIK.U.5    | LS      | BA   |
| LS.BA.CP.U.2     | LS      | BA   |
| LS.BA.ML.U.2     | LS      | BA   |
| LS.BA.WBO.U.3    | LS      | BA   |
| LS.BO.WL3.M.4    | LS      | BO   |
| LS.BO.I1.U.1     | LS      | BO   |
| LS.BO.I2.U.1     | LS      | BO   |
| LS.BO.WL2.M.2    | LS      | BO   |
| LS.BO.R1.U.2     | LS      | BO   |
| LS.BO.R2.U.4     | LS      | BO   |
| LS.BO.R3.U.3     | LS      | BO   |
| LS.BO.HC4.M.3    | LS      | BO   |
| LS.BO.LC4.M.2    | LS      | BO   |
| LS.LA.VET.M.4    | LS      | LA   |
| LS.LA.SSV.M.1    | LS      | LA   |
| LS.LA.NAV.M.4    | LS      | LA   |
| LS.LA.SHO.M.2    | LS      | LA   |
| LS.LA.WES.M.3    | LS      | LA   |
| LS.LA.GLO.M.3    | LS      | LA   |
| LS.LA.HOW.U.5    | LS      | LA   |
| LS.LA.SAN.U.2    | LS      | LA   |
| LS.LA.ARR.U.2    | LS      | LA   |
| LS.MN.L06-LS.U.2 | LS      | MN   |
| LS.MN.L06-LS.U.5 | LS      | MN   |
| LS.MN.L07-LS.U.4 | LS      | MN   |
| LS.MN.L08-LS.U.5 | LS      | MN   |
| LS.MN.L09-LS.U.3 | LS      | MN   |
| LS.MN.L01-LS.M.4 | LS      | MN   |
| LS.MN.L01-LS.U.3 | LS      | MN   |
| LS.MN.L02-LS.U.1 | LS      | MN   |
| LS.MN.L05-LS.U.2 | LS      | MN   |
| LS.PX.MON.U.2    | LS      | PX   |
| LS.PX.PKW.U.5    | LS      | PX   |
| LS.PX.PIE.M.4    | LS      | PX   |
| LS.PX.ALA.M.3    | LS      | PX   |
| LS.PX.PAL.M.3    | LS      | PX   |
| LS.PX.MAN.M.2    | LS      | PX   |
| LS.PX.NUE.M.1    | LS      | PX   |
| LS.PX.ENC.M.4    | LS      | PX   |
| LS.PX.COY.M.3    | LS      | PX   |
| PA.BA.PIK.U.1    | PA      | BA   |
| PA.BA.LH-3.M.2   | PA      | BA   |
| PA.BA.LH-3.M.3   | PA      | BA   |
| PA.BA.WB.U.1     | PA      | BA   |
| PA.BA.AA.U.1     | PA      | BA   |
| PA.BA.WGP.M.3    | PA      | BA   |
| PA.BA.LL-4.M.3   | PA      | BA   |
| PA.BA.LA.U.2     | PA      | BA   |
| PA.BA.LH-2.M.2   | PA      | BA   |
| PA.BA.W3.M.3     | PA      | BA   |
| PA.BA.RG-1.M.2   | PA      | BA   |
| PA.BA.LL-3.M.5   | PA      | BA   |
| PA.BO.I2.U.3     | PA      | BO   |
| PA.BO.HC1.M.4    | PA      | BO   |
| PA.BO.R3.U.2     | PA      | BO   |
| PA.BO.HC4.M.5    | PA      | BO   |
| PA.BO.R4.U.2     | PA      | BO   |
| PA.BO.WL2.M.5    | PA      | BO   |
| PA.BO.WL4.M.4    | PA      | BO   |
| PA.BO.LC4.M.4    | PA      | BO   |
| PA.BO.HC2.M.1    | PA      | BO   |
| PA.BO.R1.U.2     | PA      | BO   |
| PA.BO.WL1.M.1    | PA      | BO   |
| PA.BO.I1.U.5     | PA      | BO   |
| PA.LA.ALL.M.5    | PA      | LA   |
| PA.LA.SEP.M.1    | PA      | LA   |
| PA.LA.SEP.M.5    | PA      | LA   |
| PA.LA.WHI.M.2    | PA      | LA   |
| PA.LA.ROS.M.5    | PA      | LA   |
| PA.LA.LUB.M.2    | PA      | LA   |
| PA.LA.GLO.M.2    | PA      | LA   |
| PA.LA.ZOO.M.4    | PA      | LA   |
| PA.LA.ZOO.M.5    | PA      | LA   |
| PA.LA.NWH.M.2    | PA      | LA   |
| PA.LA.KIN.M.4    | PA      | LA   |
| PA.LA.POP.M.4    | PA      | LA   |
| PA.PX.BUF.M.3    | PA      | PX   |
| PA.PX.PIE.M.4    | PA      | PX   |
| PA.PX.LAP.M.5    | PA      | PX   |
| PA.PX.ALA.M.1    | PA      | PX   |
| PA.PX.PAP.M.2    | PA      | PX   |
| PA.PX.PAP.M.5    | PA      | PX   |
| PA.PX.DF1.M.2    | PA      | PX   |
| PA.PX.RPP.U.3    | PA      | PX   |
| PA.PX.ENC.M.4    | PA      | PX   |
| PA.PX.ENC.M.5    | PA      | PX   |
| PA.PX.COY.M.1    | PA      | PX   |
| PA.PX.BUF.M.2    | PA      | PX   |
| TO.BA.WBO.U.4    | TO      | BA   |
| TO.BA.CP.U.1     | TO      | BA   |
| TO.BA.FH.U.1     | TO      | BA   |
| TO.BA.LH-3.M.4   | TO      | BA   |
| TO.BA.WGP.M.3    | TO      | BA   |
| TO.BA.GA.U.4     | TO      | BA   |
| TO.BA.PIK.U.4    | TO      | BA   |
| TO.BA.PSP.M.1    | TO      | BA   |
| TO.BA.RG-2.M.2   | TO      | BA   |
| TO.BO.HC1.M.4    | TO      | BO   |
| TO.BO.HC2.M.5    | TO      | BO   |
| TO.BO.HC3.M.1    | TO      | BO   |
| TO.BO.HC4.M.5    | TO      | BO   |
| TO.BO.LC1.M.1    | TO      | BO   |
| TO.BO.LC2.M.5    | TO      | BO   |
| TO.BO.LC3.M.1    | TO      | BO   |
| TO.BO.WL2.M.1    | TO      | BO   |
| TO.BO.I2.U.3     | TO      | BO   |
| TO.LA.WHI.M.5    | TO      | LA   |
| TO.LA.HAR.M.4    | TO      | LA   |
| TO.LA.MR1.M.1    | TO      | LA   |
| TO.LA.GLO.M.5    | TO      | LA   |
| TO.LA.ZOO.M.1    | TO      | LA   |
| TO.LA.NWH.M.4    | TO      | LA   |
| TO.LA.VNS.M.2    | TO      | LA   |
| TO.LA.PEP.M.5    | TO      | LA   |
| TO.LA.COM.M.4    | TO      | LA   |
| TO.MN.L11-TO.M.3 | TO      | MN   |
| TO.MN.L02-TO.U.1 | TO      | MN   |
| TO.MN.L04-TO.U.1 | TO      | MN   |
| TO.MN.L06-TO.U.2 | TO      | MN   |
| TO.MN.L08-TO.U.5 | TO      | MN   |
| TO.MN.L09-TO.U.2 | TO      | MN   |
| TO.MN.L11-TO.U.3 | TO      | MN   |
| TO.MN.L05-TO.M.5 | TO      | MN   |
| TO.MN.L08-TO.M.5 | TO      | MN   |
| TO.PX.BUF.M.1    | TO      | PX   |
| TO.PX.ALA.M.2    | TO      | PX   |
| TO.PX.LAP.M.4    | TO      | PX   |
| TO.PX.WES.M.1    | TO      | PX   |
| TO.PX.CLA.M.1    | TO      | PX   |
| TO.PX.DF1.M.1    | TO      | PX   |
| TO.PX.DF2.M.1    | TO      | PX   |
| TO.PX.COY.M.1    | TO      | PX   |
| TO.PX.COY.M.6    | TO      | PX   |

Subset of samples used in SNP catalog creation.

# 7 Stacks: Metapopulation locus matching with `sstacks`

Files can be found in the `07_sstacks/` directory.

All samples in the population (or all samples you want to include in the
analysis) are matched against the catalog produced in
[`cstacks`](#step-6---cstacks) with `sstacks`, run in script
`stacks_SPECIES.sh` and `stacks_SPECIES_additional.sh`. It runs off of
the samples based in the output directory *and* the listed samples in
`sstacks_samples_SPECIES.txt` and
`sstacks_samples_SPECIES_additional.txt` (respectively), so make sure
all your files (sample and catalog, etc.) are there and match.
`sstacks_samples_SPECIES.txt` takes the form:

    DS.BA.GA.U.1
    DS.BA.GA.U.2
    DS.BA.GA.U.3
    DS.BA.GA.U.4
    DS.BA.GA.U.5
    ...

There should be a new file produced at this step for every sample in the
output directory:

- `<samplename>.matches.tsv.gz`

A small number of samples generated very few matches to the catalog
(such as only 4 loci matching, obviously not enough to draw any
conclusions) and therefore aren’t used in the next step. See
[output/sstacks-discarded_samples.csv](output/sstacks-discarded_samples.csv).

# 8 Genotype probabilities with `polyRAD`

Files can be found in the `08_polyRAD/` directory.

## 8.1 Make `RADdata` object

We used the [polyRAD package](https://github.com/lvclark/polyRAD) to
call genotypes because many of our species are polyploid or have
historical genome duplication. PolyRAD takes the catalog output
(`catalog.alleles.tsv.gz`) and accompanying matches to the catalog
(e.g., `CD.BA.AA.U.1.matches.tsv.gz`) to create genotype likelihoods for
species with diploidy and/or polyploidy.

We used the catalog and match files to create a RADdata object class in
R for each species. We ran this on the Rockfish compute cluster at Johns
Hopkins University, with the `make_polyRAD_<spp>.R` script doing the
brunt of the work. The R script was wrapped by `polyrad_make_<spp>.sh`
to submit the script to the SLURM scheduler.

*Relevant Parameters:*

- `min.ind.with.reads` was set to 20% of samples. This means we
  discarded any loci not found in at least 20% of samples for each
  species.
- `min.ind.with.minor.allele` was set to `2`. This means a locus must
  have at least this many samples with reads for the minor allele in
  order to be retained.

*Requires:*

- `popmap_<spp>_polyrad.txt`, a list of samples and population
- output from sstacks

*Outputs:*

- `<spp>_polyRADdata.rds`, RDS object (the RADdata object)

## 8.2 Calculate overdispersion

Next, we calculated overdispersion using the
`polyRAD_overdispersion_<spp>.R` script, wrapped by
`polyrad_overd_<spp>.sh` to submit the script to the SLURM scheduler.

Requires:

- `popmap_<spp>_polyrad.txt`, a list of samples and population
- `<spp>_polyRADdata.rds`, RDS object (the RADdata object) output from
  the previous step

Outputs:

- `<spp>_overdispersion.rds`, RDS object (the overdispersion test
  output)

## 8.3 Estimate genotypes

Next, we calculated filtered loci based on the expected Hind/He
statistic and estimated population structure/genotypes using the
`polyRAD_filter_<spp>.R` script, wrapped by `polyrad_filt_<spp>.sh` to
submit the script to the SLURM scheduler.

We used the table in [this
tutorial](https://lvclark.r-universe.dev/articles/polyRAD/polyRADtutorial.html),
which estimated an inbreeding based on the ploidy, optimal
overdispersion value, and mean Hind/He. These values are hardcoded in
`polyRAD_filter_<spp>.R`.

Requires:

- `popmap_<spp>_polyrad.txt`, a list of samples and population
- `<spp>_polyRADdata.rds`, RDS object (the RADdata object) output from
  the previous step
- `<spp>_overdispersion.rds`, RDS object (the overdispersion test
  output) output from the previous step

Outputs:

- `<spp>_filtered_RADdata.rds`, RDS object (RADdata object filtered for
  appropriate Hind/He)
- `<spp>_IteratePopStructPCA.csv`, data output from the genotype
  estimate PCA, suitable for plotting
- `<spp>_estimatedgeno_RADdata.rds`, RDS object (RADdata object with
  genotype estimates)

## 8.4 Final filter and file cleanup

The output `<spp>_estimatedgeno_RADdata.rds` needs to be converted to
genind and structure format for further analysis and steps. There is a
little cleanup involved so the population information is retained. For
example, Structure needs the population identity to be an integer, not a
string. This set of functions can be run on a laptop.

At this stage, we also visually assessed the $H_{ind}/H_e$ statistic
versus the locus depth (see `check_coverage` inside the
`convert_genomics.R` script). We removed the following samples from
further analysis:

| Sample         |
|:---------------|
| CD.BA.PSP.M.1  |
| CD.BA.DHI.U.2  |
| CD.BA.DHI.U.3  |
| CD.BA.RG-1.M.5 |
| CD.BA.RG-1.M.4 |
| DS.BO.WL1.M.4  |
| DS.BO.I1.U.3   |
| EC.BO.R4.U.1   |
| LS.BO.HC2.M.5  |
| LS.BO.LC4.M.3  |
| LS.BO.R2.U.4   |
| LS.BO.R2.U.1   |
| PA.BA.LH-3.M.4 |
| PA.BA.AA.U.3   |
| PA.BA.AA.U.4   |
| PA.PX.RPP.U.2  |
| PA.BO.HC2.M.4  |
| PA.PX.RPP.U.1  |
| TO.BA.TRC.U.1  |
| TO.BA.TRC.U.3  |
| TO.BO.R4.U.1   |
| TO.BA.TRC.U.2  |
| TO.BO.R4.U.2   |
| TO.BO.R2.U.2   |

Subset of samples discarded after genotype estimation using polyRAD.

``` r
source("08_polyRAD/convert_genomics.R")
```

``` r
convert_all()
```

# 9 `Structure` Analysis

Files, inlcuding model parameters, can be found in the `09_structure/`
directory.

Structure documentation can be found
[here](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_doc.pdf).

## 9.1 Running Structure on `polyRAD` output

`polyRAD` outputs genotype probabilities in a format suitable for
Structure. These files were named as:

    CD_estimatedgeno.structure
    DS_estimatedgeno.structure
    EC_estimatedgeno.structure
    LS_estimatedgeno.structure
    PA_estimatedgeno.structure
    TO_estimatedgeno.structure

We ran all species using a naive approach (not using prior information)
with $K={1,2,3,4,5}$ (`MAXPOPS` argument). To search for the most
appropriate K, We ran Structure through 5 replicate runs for each
combination of species and K, with 10000 iterations discarded as burn-in
and retained 20000 iterations. These runs created files that look like:

    structure_out_CD1_naive_f         // K = 1, rep 1
    structure_out_CD1_naive_rep2_f    // K = 1, rep 2
    structure_out_CD1_naive_rep3_f    // K = 1, rep 3
    structure_out_CD1_naive_rep4_f    // K = 1, rep 4
    structure_out_CD1_naive_rep5_f    // K = 1, rep 5
    structure_out_CD2_naive_f         // K = 2, rep 1
    structure_out_CD2_naive_rep2_f    // K = 2, rep 2
    structure_out_CD2_naive_rep3_f    // K = 2, rep 3
    ...

Within each species, we compressed the result files for all K and reps
and submitted to [Structure
Harvester](https://taylor0.biology.ucla.edu/structureHarvester/) to
choose the optimal K using the Delta-K method (see
<https://link.springer.com/article/10.1007/s12686-011-9548-7>). Once the
optimal K was selected per species, we re-ran Structure using a greater
number of iterations (100000) for final output and plotting.

# 10 Conceptual Figure

We made a conceptual figure to help readers understand PCA patterns we
might expect.

``` r
source("R/10-Fig1-conceptual_fig.R")
```

``` r
make_fig1()
```

# 11 NLCD Data and Site Plots

NLCD is used in Fig. 2 and IBE analysis, described below.

From the USGS:

> The U.S. Geological Survey (USGS), in partnership with several federal
> agencies, has developed and released four National Land Cover Database
> (NLCD) products over the past two decades: NLCD 1992, 2001, 2006, and
> 2011. This one is for data from 2016 and describes urban
> imperviousness.

> <https://www.mrlc.gov/data/type/urban-imperviousness>

> NLCD imperviousness products represent urban impervious surfaces as a
> percentage of developed surface over every 30-meter pixel in the
> United States. NLCD 2016 updates all previously released versions of
> impervious products for CONUS (NLCD 2001, NLCD 2006, NLCD 2011) along
> with a new date of impervious surface for 2016. New for NLCD 2016 is
> an impervious surface descriptor layer. This descriptor layer
> identifies types of roads, core urban areas, and energy production
> sites for each impervious pixel to allow deeper analysis of developed
> features.

## 11.1 Preparing NLCD Data

First, we trimmed the large data. This makes a smaller `.rds` file for
each city.

``` r
source("R/10-trim_NLCD_spatial_data.R")
```

``` r
create_spatial_rds_files()
create_spatial_rds_files(spp = "CD")
```

## 11.2 Climate normals

We obtained climate normals for plotting from
<https://www.ncei.noaa.gov/access/us-climate-normals>. We used the
latest 30-year period (1991-2020): Most recent standard climatological
period (2021 release); which is recommended for most purposes.

## 11.3 Maps of sampling locations

Next, we made plots for each city’s sampling locations. Note that these
only include sites that had viable polymorphic loci.

``` r
source("R/10-Fig2-plot_map_of_samples.R")
```

``` r
# Plot Fig 2
make_all_urban_site_plots_with_clim_normals()
```

# 12 Principal components analysis & plots

The following creates PCA plots from polyRAD data.

``` r
source("R/11-plot_pca.R")
```

``` r
# Plot Fig 4
make_pca_city_and_pctimp_all()
```

In addition to coloring points by city, we also colored points by NLCD
Urban Cover %. Note that in the pdf version of this document, the figure
might appear on the next pages.

``` r
make_pca_pctimp_only_all()
```

<figure>
<img src="Supplement_files/figure-gfm/plot_pca_3-1.png"
alt="PCA colored by % Impervious surface." />
<figcaption aria-hidden="true">PCA colored by % Impervious
surface.</figcaption>
</figure>

# 13 Genetic structure using Structure and sNMF

## 13.1 `Structure` optimal K

Within each species, we compressed the result files for all K and reps
and submitted to [Structure
Harvester](https://taylor0.biology.ucla.edu/structureHarvester/) to
choose the optimal K using the Delta-K method (see
<https://link.springer.com/article/10.1007/s12686-011-9548-7>).

The results were:

CD: K=3 DS: K=3 EC: K=2 LS: K=3 PA: K=4 TO: K=3

``` r
# This file contains output from various K from Structure..
read_csv("output/structure/structure_k_Pr.csv")
```

The code below generates plots of various K (e.g., K={1-5}) vs
likelihood, but did not end up being used in the manuscript.

``` r
source("R/12-structure_k.R")
```

## 13.2 Plotting Structure output

The code below generates plots for Structure results.

``` r
source("R/12-plot_structure.R")
```

``` r
make_structure_multi_plot()
```

## 13.3 Validation of Structure results with sNMF

We ran sNMF as an alternative to Structure to validate the results. We
coerced all polyploid data to diploid data to make the file types
compatible with the sNMF function in R. The snmf() function computes an
entropy criterion that evaluates the quality of fit of the statistical
model to the data by using a cross-validation technique. We plotted the
cross-entropy criterion for K=\[2:10\] for all species. Using the best
K, we then selected the best of 10 runs in each K using the
`which.min()` function.

``` r
source("R/12-sNMF.R")
```

The following runs sNMF and generates the figure. .

``` r
do_all_sNMF()
```

<figure>
<img src="Supplement_files/figure-gfm/snmf_2-1.png"
alt="Ancestry coefficients obtained using snmf(). As with the Structure analysis, E. canadensis (horseweed) and L. serriola (prickly lettuce) appear to have the most population structure. D. sanguinalis, E. canadensis, and L. serriola from Phoenix appear unique. In general, sNMF produced larger K for most species, which will create more sensitivity to admixture." />
<figcaption aria-hidden="true">Ancestry coefficients obtained using
<code>snmf()</code>. As with the Structure analysis, <em>E.
canadensis</em> (horseweed) and <em>L. serriola</em> (prickly lettuce)
appear to have the most population structure. <em>D. sanguinalis</em>,
<em>E. canadensis</em>, and <em>L. serriola</em> from Phoenix appear
unique. In general, sNMF produced larger K for most species, which will
create more sensitivity to admixture.</figcaption>
</figure>

# 14 Genetic Differentiation

## 14.1 Among city - Jost’s D, $G_{ST}$, $F_{ST}$

We used `polyrad::calcPopDiff()` to calculate continental population
statistics for each species.

``` r
source("R/13-calc_popdiff_stats.R")
```

``` r
do_all_continental_stats()
```

``` r
# CD as an example
read.csv("output/population_stats/popdiff_stats_CD.csv")
```

    ##   X statistic      value
    ## 1 1     JostD 0.30579677
    ## 2 2       Gst 0.02735719
    ## 3 3       Fst 0.02812163

## 14.2 Among city - pairwise $\rho$

We used [GenoDive v. 3.0.6](https://doi.org/10.1111/1755-0998.13145) to
calculate pairise $\rho$ (rho) among cities within species. Note that
there is a p-value correction for testing multiple cities (species are
treated as independent, however).

This can be run in GenoDive by selecting Pairwise Differentiation from
the Analysis menu and selecting the “rho” statistic from the dropdown.

We used the following script to clean up the results.

``` r
source("R/13-rho.R")
compile_rho_table()
```

| Species | City1 | City2 | rho | p-value | adjusted p-value |
|:---|:---|:---|---:|---:|:---|
| CD | PX | BA | 0.050 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| CD | LA | BA | 0.046 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| CD | PX | LA | 0.015 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| DS | MN | BA | 0.031 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| DS | BO | BA | 0.018 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| DS | PX | BA | 0.012 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| DS | MN | BO | 0.007 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| DS | PX | BO | -0.002 | 0.875 | <span style="     ">0.955</span> |
| DS | PX | MN | -0.002 | 0.955 | <span style="     ">0.955</span> |
| EC | PX | BA | 0.098 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| EC | PX | LA | 0.087 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| EC | LA | BA | 0.038 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.001</span> |
| LS | PX | BA | 0.077 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | PX | MN | 0.069 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | PX | LA | 0.061 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | PX | BO | 0.056 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | MN | LA | 0.039 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | LA | BA | 0.038 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | BO | BA | 0.032 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | MN | BO | 0.021 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | LA | BO | 0.010 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0011</span> |
| LS | MN | BA | 0.009 | 0.002 | <span style=" font-weight: bold;   text-decoration: underline; ">0.002</span> |
| PA | PX | BO | 0.028 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| PA | LA | BO | 0.024 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| PA | PX | BA | 0.015 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| PA | LA | BA | 0.011 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0015</span> |
| PA | BO | BA | 0.008 | 0.002 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0024</span> |
| PA | PX | LA | -0.002 | 0.972 | <span style="     ">0.972</span> |
| TO | PX | BA | 0.023 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | PX | MN | 0.015 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | PX | BO | 0.013 | 0.002 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0025</span> |
| TO | LA | BA | 0.011 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | LA | BO | 0.009 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | MN | LA | 0.009 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | PX | LA | 0.009 | 0.027 | <span style=" font-weight: bold;   text-decoration: underline; ">0.03</span> |
| TO | BO | BA | 0.008 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | MN | BO | 0.008 | 0.001 | <span style=" font-weight: bold;   text-decoration: underline; ">0.0014</span> |
| TO | MN | BA | 0.001 | 0.098 | <span style="     ">0.098</span> |

Rho statistics for pairwise comparison between cities.

## 14.3 Within city - allelic richness

We used [GenoDive v. 3.0.6](https://doi.org/10.1111/1755-0998.13145) to
calculate several additional statistics.

This can be run in GenoDive by selecting Analysis \> Genetic Diversity,
and selecting “Calculate indices separately for every population” and
selecting “Correct for unknown dosage of alleles” for the polyploid
species.

- Num: Number of alleles
- Eff_num: Effective number of alleles
- Ho: Observed Heterozygosity
- Hs: Heterozygosity within populations
- Gis: Inbreeding coefficient

``` r
head(read.csv("output/population_stats/genodive_genetic_diversity.csv"))
```

    ##   spp city    Num Eff_num    Ho    Hs   Gis
    ## 1  CD   BA  8.242   4.530 0.627 0.748 0.163
    ## 2  CD   LA 12.664   6.470 0.683 0.827 0.174
    ## 3  CD   PX 14.192   6.978 0.672 0.830 0.191
    ## 4  DS   BA 10.611   5.069 0.637 0.768 0.171
    ## 5  DS   BO 10.588   5.230 0.612 0.778 0.214
    ## 6  DS   MN 11.696   5.233 0.623 0.771 0.192

## 14.4 Within city - $F_{IS}$ (homozygosity within population)

We used [GenoDive v. 3.0.6](https://doi.org/10.1111/1755-0998.13145) to
calculate $F_{IS}$. This gives a good estimate of whether there are more
homozygotes than expected (positive number) or more heterozygotes than
expected (negative number). Notably, GenoDive accommodates polyploids
and reduces the bias on $F_{IS}$ by performing a permutation test. By
default, there are 999 permutations.

This can be run in GenoDive by selecting Analysis \> Hardy-Weinberg \>
Heterozygosity-based (Nei) method.

``` r
head(read.csv("output/population_stats/genodive_output_Fis.csv"))
```

    ##   Species Population  n   Fis
    ## 1      CD         BA 55 0.166
    ## 2      CD         LA 48 0.186
    ## 3      CD         PX 82 0.200
    ## 4      CD    Overall NA 0.187
    ## 5      DS         BA 55 0.208
    ## 6      DS         BO 52 0.252

## 14.5 Within city - $\bar{r}_d$ - Linkage disequilibrium

We used `poppr::ia()` to calculate the standardized index of association
of loci in the dataset ($\bar{r}_d$ or `rbarD`). We use the standardized
index of association to avoid the influence of different sample sizes,
as described by [Agapow and Burt
2001](https://doi.org/10.1046/j.1471-8278.2000.00014.x).

When `p.rD` is small (\<0.05) and rbarD is (relatively) higher, that is
a sign that the population could be in linkage disequilibrium.

An interesting note from the documentation:

> It has been widely used as a tool to detect clonal reproduction within
> populations. Populations whose members are undergoing sexual
> reproduction, whether it be selfing or out-crossing, will produce
> gametes via meiosis, and thus have a chance to shuffle alleles in the
> next generation. Populations whose members are undergoing clonal
> reproduction, however, generally do so via mitosis. This means that
> the most likely mechanism for a change in genotype is via mutation.
> The rate of mutation varies from species to species, but it is rarely
> sufficiently high to approximate a random shuffling of alleles. The
> index of association is a calculation based on the ratio of the
> variance of the raw number of differences between individuals and the
> sum of those variances over each locus. You can also think of it as
> the observed variance over the expected variance.

There is a nice description
[here](https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html).

``` r
source("R/13-rbarD.R")
calc_rbarD()
```

``` r
head(read.csv("output/population_stats/rbarD.csv"))
```

    ##   spp city  n       Ia  p.Ia     rbarD  p.rD
    ## 1  CD   BA 55 664.7655 0.001 0.2950534 0.001
    ## 2  CD   LA 48 470.5913 0.001 0.2070064 0.001
    ## 3  CD   PX 82 634.5787 0.001 0.2792566 0.001
    ## 4  DS   BA 55 557.7334 0.001 0.2123881 0.001
    ## 5  DS   BO 52 896.0906 0.001 0.3398873 0.001
    ## 6  DS   MN 81 578.2494 0.001 0.2192197 0.001

## 14.6 Within city - private alleles

We used a custom script to calculate % private alleles as the percentage
of total alleles unique to a particular city.

``` r
source("R/13-private_alleles.R")
get_all_private_alleles()
```

## 14.7 AMOVA

We performed hierarchical analysis of molecular variance (AMOVA; using
GenoDive 3.06) based on the Rho-statistics, which is based on a Ploidy
independent Infinite Allele Model. AMOVA is under the “Analysis” menu.
Structure format files typically have only one level of population
grouping. To ensure nestedness, we added a “City” level manually by
going to “Data \> Population grouping ..” menu and adding “City”. When
running the AMOVA, we selected “Advanced” and selected Individual nested
within Population nested within City. We used 999 permutations.

For bluegrass (PA), Minneapolis (MN) had a small sample size (2
individuals). We re-ran the AMOVA with these samples excluded but found
similar results.

| Species | Source of Variation | Nested in | SSD | d.f. | MS | Var-comp | %Var | F-value | P-value |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|:---|
| CD | Within Population | – | 34441.923 | 136 | 253.249 | 253.249 | 0.235 | 0.765 | – |
| CD | Among Population | City | 137243.842 | 46 | 2983.562 | 732.212 | 0.678 | 0.743 | 0.001 |
| CD | Among City | – | 18722.905 | 2 | 9361.452 | 94.155 | 0.087 | 0.087 | 0.001 |
| DS | Within Population | – | 39458.473 | 155 | 254.571 | 254.571 | 0.224 | 0.776 | – |
| DS | Among Population | City | 149988.556 | 65 | 2307.516 | 638.534 | 0.561 | 0.715 | 0.001 |
| DS | Among City | – | 47958.337 | 3 | 15986.112 | 244.451 | 0.215 | 0.215 | 0.001 |
| EC | Within Population | – | 15832.596 | 73 | 216.885 | 216.885 | 0.302 | 0.698 | – |
| EC | Among Population | City | 44957.413 | 31 | 1450.239 | 391.507 | 0.546 | 0.644 | 0.001 |
| EC | Among City | – | 10368.036 | 2 | 5184.018 | 108.704 | 0.152 | 0.152 | 0.001 |
| LS | Within Population | – | 10889.549 | 120 | 90.746 | 90.746 | 0.308 | 0.692 | – |
| LS | Among Population | City | 32081.032 | 59 | 543.746 | 159.914 | 0.544 | 0.638 | 0.001 |
| LS | Among City | – | 8721.721 | 4 | 2180.430 | 43.569 | 0.148 | 0.148 | 0.001 |
| PA - no MN | Within Population | – | 9129.284 | 128 | 71.323 | 71.323 | 0.230 | 0.770 | – |
| PA - no MN | Among Population | City | 34474.530 | 44 | 783.512 | 194.429 | 0.628 | 0.732 | 0.001 |
| PA - no MN | Among City | – | 7968.939 | 3 | 2656.313 | 43.730 | 0.141 | 0.141 | 0.002 |
| PA - w/ MN | Within Population | – | 10167.015 | 129 | 78.814 | 78.814 | 0.278 | 0.722 | – |
| PA - w/ MN | Among Population | City | 30922.681 | 44 | 702.788 | 168.366 | 0.593 | 0.681 | 0.001 |
| PA - w/ MN | Among City | – | 6999.878 | 4 | 1749.970 | 36.750 | 0.129 | 0.129 | 0.002 |
| TO | Within Population | – | 15727.071 | 162 | 97.081 | 97.081 | 0.314 | 0.686 | – |
| TO | Among Population | City | 51744.430 | 71 | 728.795 | 202.502 | 0.655 | 0.676 | 0.001 |
| TO | Among City | – | 4740.183 | 4 | 1185.046 | 9.469 | 0.031 | 0.031 | 0.001 |

AMOVA Statistics.

The following code plots the figure in the main manuscript.

``` r
source("R/14-AMOVA.R")
make_amova_plot()
```

<!-- ## Isolation by distance and environment -->
<!-- We used [multiple matrix regression with randomization (MMRR)](https://doi.org/10.1111/evo.12134) to determine the relative contributions of isolation by distance (i.e., an association between genetic and geographic distances) and isolation by environment (i.e., an association between genetic and environmental distances). -->
<!-- ### Genetic distance -->
<!-- We calculated the genetic dissimilarity matrix (among sites) using the [Cavalli-Sforza](https://doi.org/10.2307/2406616) (Chord) distance metric in GenoDive. While this can also be done using `adegenet`, we don't want to make assumptions about ploidy. We used the "*_estimatedgeno_sitesaspops.structure" files so that sites (not cities) were treated here as populations. -->
<!-- ### Geographic distance -->
<!-- We took the traditional approach to creating a geographic dissimilarity matrix (based on latitude and longitude) using euclidean distance. -->
<!-- ### Environmental data and distance -->
<!-- Environmental variables include the monthly averages in the middle of the day for: -->
<!-- - air temperature at 5cm above ground -->
<!-- - air temperature at 1.2m above ground -->
<!-- - soil temperature at 2.5cm below ground -->
<!-- - RH (relative humidity) at 5cm above ground -->
<!-- - RH at 1.2m above ground -->
<!-- Variables were extracted from historic datasets and modeled using a microclimate model. More information can be found on the [NicheMapR website](https://mrke.github.io/) (how the model works, what variables can be manipulated and what you can model, vignettes for running models in R). -->
<!-- This method was chosen because it takes data from global datasets (you can use both historic and current or pick specific years) but then accounts for site-specific variables (we can change the % shade, the slope or aspect of the landscape, and it considers elevation, average cloud cover, etc.). [Here’s the list](https://mrke.github.io/models/MicroClimate-Models) of all the different models/datasets we’re able to can pull from. It's meant for mechanistic niche modeling. -->
<!-- Variables in the file `site_data_DUC_environvars.csv` are all for the monthly averages at noon (12pm - hottest part of the day!) and are extreme. In other words, they are maximums. -->
<!-- Note that [this Stack Overflow post](https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r) is helpful with installing `NicheMapR`. -->
<!-- ```{r} -->
<!-- # devtools::install_github('mrke/NicheMapR') -->
<!-- # library(NicheMapR) -->
<!-- # -->
<!-- # test_site_coords <- c(sites[1,]$lat, sites[1,]$long) -->
<!-- # test_distance_to_city_center_km <- sites[1,]$distance_to_city_center_km -->
<!-- # micros_ <- micro_usa(loc = test_site_coords) -->
<!-- # -->
<!-- # loc <- c(-89.40, 43.07) -->
<!-- # micro <- micro_global(loc = loc) -->
<!-- ``` -->
<!-- Environmental distance was generated the same way as geographic distance above (euclidean distance). These are the four environmental variables mentioned in the main manuscript, although more environmental variables are present in the raw data: % Urban cover, Distance to city center, April soil temperature, and July soil temperature. In the raw data these appear as `nlcd_urban_pct`, `distance_to_city_center`, `soiltemp_2.5cm_Apr_12pm`, `soiltemp_2.5cm_Jul_12pm`. -->
<!-- ### Analysis -->
<!-- Note that for this analysis, we treated each *sampling site* as a distinct location. There would not be enough power to do a distance matrix among 3-5 cities. Code for generating matrices, running the MMRR, and generating figures can be found in the source code below. -->
<!-- ```{r mmrr_1, message=FALSE} -->
<!-- source("R/IBD_IBE_MMRR.R") -->
<!-- ``` -->
<!-- ```{r mmrr_2, eval = FALSE, message=FALSE} -->
<!-- make_mmrr_plot() -->
<!-- ``` -->
<!-- Below are the results of the MMRR, with all cities in the same model. Species are treated as independent, separate models. -->
<!-- ```{r mmrr_3a, echo = FALSE, message=FALSE} -->
<!-- table_ <- readr::read_csv("output/MMRR/MMRR_total_overall_models.csv") %>% -->
<!--   mutate(p = kableExtra::cell_spec( -->
<!--     `F p-value:`, -->
<!--     bold = ifelse(as.numeric(`F p-value:`) < 0.05, TRUE, FALSE), -->
<!--     underline = ifelse(as.numeric(`F p-value:`) < 0.05, TRUE, FALSE) -->
<!--   )) -->
<!-- knitr::kable( -->
<!--   table_, -->
<!--   format = "simple", -->
<!--   caption = "Overall model p-values from MMRR." -->
<!-- ) -->
<!-- ``` -->
<!-- ```{r mmrr_3b, echo = FALSE, message=FALSE} -->
<!-- table_ <- -->
<!--   readr::read_csv("output/MMRR/MMRR_total.csv") %>% -->
<!--   mutate(p = kableExtra::cell_spec( -->
<!--     p, -->
<!--     bold = ifelse(as.numeric(p) < 0.05, TRUE, FALSE), -->
<!--     underline = ifelse(as.numeric(p) < 0.05, TRUE, FALSE) -->
<!--   )) -->
<!-- knitr::kable( -->
<!--   table_, -->
<!--   format = "simple", -->
<!--   caption = "Full parameter estimates and statistics by species from running 9999 permutations ('Reps') via MMRR." -->
<!-- ) -->
<!-- ``` -->
<!-- We also repeated this within city. -->
<!-- ```{r mmrr_4a, echo = FALSE, message=FALSE} -->
<!-- table_ <- readr::read_csv("output/MMRR/MMRR_splitbycity_overall_models.csv") %>% -->
<!--   mutate(`p adjusted` = kableExtra::cell_spec(`p adjusted`, -->
<!--   bold = ifelse(as.numeric(`p adjusted`) < 0.05, TRUE, FALSE), -->
<!--   underline = ifelse(as.numeric(`p adjusted`) < 0.05, TRUE, FALSE) -->
<!-- )) -->
<!-- knitr::kable( -->
<!--   table_, -->
<!--   format = "simple", -->
<!--   caption = "Overall model p-values from MMRR, subset by each city. Adjusted p-values are corrected using the Benjamini and Hochberg method within each species' results." -->
<!-- ) -->
<!-- ``` -->
<!-- ```{r mmrr_4b, echo = FALSE, message=FALSE} -->
<!-- table_ <- readr::read_csv("output/MMRR/MMRR_splitbycity.csv") %>% -->
<!--   mutate(p = kableExtra::cell_spec( -->
<!--     p, -->
<!--     bold = ifelse(as.numeric(p) < 0.05, TRUE, FALSE), -->
<!--     underline = ifelse(as.numeric(p) < 0.05, TRUE, FALSE) -->
<!--   )) -->
<!-- knitr::kable( -->
<!--   table_, -->
<!--   format = "simple", -->
<!--   caption = "Parameter estimates and statistics from running 9999 permutations ('Reps') via MMRR, subset by each city. Note that p-values are not adjusted for multiple testing." -->
<!-- ) -->
<!-- ``` -->
<!-- ## Sites per city per species -->
<!-- Use `get_unique_site_post_genotyping.R`. -->

# 15 `SessionInfo()`

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.4.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] polysat_1.7-7     LEA_3.16.0        ggh4x_0.2.8       here_1.0.1       
    ##  [5] lubridate_1.9.3   forcats_1.0.0     purrr_1.0.2       tibble_3.2.1     
    ##  [9] tidyverse_2.0.0   viridis_0.6.5     viridisLite_0.4.2 raster_3.6-26    
    ## [13] sp_2.1-4          cowplot_1.2.0     stringr_1.5.1     readr_2.1.5      
    ## [17] polyRAD_2.0.0     dplyr_1.1.4       magrittr_2.0.3    tidyr_1.3.1      
    ## [21] ggplot2_4.0.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-4    gtable_0.3.6       xfun_0.52          lattice_0.22-6    
    ##  [5] tzdb_0.4.0         vctrs_0.6.5        tools_4.4.2        generics_0.1.3    
    ##  [9] parallel_4.4.2     fansi_1.0.6        highr_0.11         pkgconfig_2.0.3   
    ## [13] RColorBrewer_1.1-3 S7_0.2.0           lifecycle_1.0.4    compiler_4.4.2    
    ## [17] farver_2.1.2       textshaping_0.4.0  terra_1.8-60       codetools_0.2-20  
    ## [21] htmltools_0.5.8.1  yaml_2.3.10        pillar_1.9.0       crayon_1.5.2      
    ## [25] commonmark_1.9.1   tidyselect_1.2.1   digest_0.6.35      stringi_1.8.4     
    ## [29] labeling_0.4.3     rprojroot_2.0.4    fastmap_1.2.0      grid_4.4.2        
    ## [33] cli_3.6.3          dichromat_2.0-0.1  utf8_1.2.4         withr_3.0.2       
    ## [37] scales_1.4.0       bit64_4.0.5        timechange_0.3.0   rmarkdown_2.27    
    ## [41] ggtext_0.1.2       bit_4.0.5          gridExtra_2.3      ragg_1.3.2        
    ## [45] hms_1.1.3          kableExtra_1.4.0   evaluate_1.0.5     knitr_1.47        
    ## [49] markdown_1.13      rlang_1.1.4        gridtext_0.1.5     Rcpp_1.0.12       
    ## [53] glue_1.8.0         xml2_1.3.6         formatR_1.14       svglite_2.1.3     
    ## [57] rstudioapi_0.16.0  vroom_1.6.5        R6_2.5.1           systemfonts_1.1.0

# Appendix

## 15.1 File Organization

All data files for the Macrosystems project are permanently stored under
Meghan Avolio’s group resources in the Johns Hopkins University Rockfish
computing cluster. Files are stored under the ‘data’ directory under the
following subdirectories:

- **`01-raw_data`**: This folder contains the raw, unprocessed data
  files that were obtained directly from the sequencing server. There
  are eight `fastq.gz` files per sublibrary that correspond to the four
  sequencing lanes for each read direction.

- **`02-concatenated_data`**: This folder contains the concatenated,
  unprocessed files for each sublibrary (i.e., the files containing the
  sequences for each lane were combined to create one file per read
  direction).

- **`03-pcr_filtered_data`**: Here, you will find the resulting data
  files from the `clone_filter` program, where pcr replicates/clones
  have been removed from the raw sequences. There are two `fq.gz` files
  per sublibrary.

- **`04-process_radtags`**: This folder contains various subdirectories
  that correspond to the `process_radtags` program that demultiplexes
  and cleans the data. The `demux_txt_files` folder contains the .txt
  files used to identify barcodes and separate out the individual
  samples from each sublibrary. The resulting data files from the
  `process-radtags` program are separated by individual and can be found
  in the relevant species folder(i.e., CD, DS, EC, LS, PA, TE, TO). Each
  individual sample has four data files; `sampleID.1.fq.gz` and
  `sampleID.2.fq.gz` correspond to the forward and reverse reads for
  each sample and `sampleID.rem,1/2.fq.gz` contain the remainder reads
  that were cleaned and removed from the data sequence.

- **`05-ustacks-denovo_data`**: This folder contains species
  subdirectories that store the resulting data files from the `ustacks`
  program for each individual. There are three files per individual;
  `sampleID.allelles.tsv.gz`, `sampleID.snps.tsv.gz`, and,
  `sampleID.tags.tsv.gz`. These files should be permanently stored here
  and copied to a new directory for any new catalogs and/or when a group
  of samples are being aligned to a new catalog.

- **`catalogs_by_city`**: For any given species within city, there is
  likely to be a slightly different set of SNPs compared to the whole
  metapopulation of five cities. We examined 24 sets of species-city
  combinations. These catalogs are permanently stored here.

- **`catalogs_by_species`**: Metapopulation catalogs are stored within
  this folder for each species. The metapopulation catalog was created
  using samples from all populations to create a national catalog.

  Some notes about catalog directories:

  - Catalogs contain three files; `catalog.alleles.tsv.gz`,
    `catalog.snps.tsv.gz`, and `catalog.tafs.tsv.gz`. If you would like
    to use the catalog on a new project, you will need to copy all three
    files to a new project folder.
  - You can determine which individuals were used to create the catalog
    by looking at the `cstacks_popmap.txt` found within each folder.
    Specifically for the metapopulation catalogs, this information is
    also found in the
    [cstacks-metapop-catalog_samples-included.csv](/output/cstacks-metapop-catalog_samples-included.csv)
  - You can determine which individuals were subsequently aligned to the
    catalog and used in the subsequent stacks analysis by looking at the
    `popmap*.txt` found within each folder.
  - Each folder also contains the relevant ustacks and stacks pipeline
    scripts and output files (i.e., from `cstacks`, `gstacks`, `stacks`,
    `tsv2bam`, and `populations`),

## 15.2 Aspera Transfer File Names

See
[data/aspera_transfer_file_names.csv](data/aspera_transfer_file_names.csv).
Preview:

``` r
readLines("data/aspera_transfer_file_names.csv", 10)
```

    ##  [1] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L001_R1_001.fastq.gz" 
    ##  [2] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L001_R2_001.fastq.gz" 
    ##  [3] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L002_R1_001.fastq.gz" 
    ##  [4] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L002_R2_001.fastq.gz" 
    ##  [5] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L003_R1_001.fastq.gz" 
    ##  [6] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L003_R2_001.fastq.gz" 
    ##  [7] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L004_R1_001.fastq.gz" 
    ##  [8] "/Hoffman_macrosystems/AMH_macro_1_1_12px_S1_L004_R2_001.fastq.gz" 
    ##  [9] "/Hoffman_macrosystems/AMH_macro_1_10_8px_S10_L001_R1_001.fastq.gz"
    ## [10] "/Hoffman_macrosystems/AMH_macro_1_10_8px_S10_L001_R2_001.fastq.gz"

## 15.3 `clone_filter` File Names

See
[data/clone_filter_file_names.csv](data/clone_filter_file_names.csv).
Preview:

``` r
readLines("data/clone_filter_file_names.csv", 10)
```

    ##  [1] "AMH_macro_1_1_12px_S1"  "AMH_macro_1_10_8px_S10" "AMH_macro_1_11_8px_S11"
    ##  [4] "AMH_macro_1_12_8px_S12" "AMH_macro_1_13_8px_S13" "AMH_macro_1_14_8px_S14"
    ##  [7] "AMH_macro_1_2_12px_S2"  "AMH_macro_1_3_12px_S3"  "AMH_macro_1_4_12px_S4" 
    ## [10] "AMH_macro_1_5_8px_S5"
