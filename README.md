# Genet assignment sensitivity analyses using ϕ and SH in trembling aspen (*Populus tremuloides*)

This repository contains the R code used to evaluate the sensitivity of genet assignments to alternative thresholds of pairwise relatedness (ϕ) and shared heterozygosity (SH) in trembling aspen (*Populus tremuloides*), as reported in the manuscript and its Supporting Information.

## Analyses

The R script performs the following analyses:

- Comparison of pairwise relatedness (ϕ) distributions within and between genets
- Comparison of shared heterozygosity (SH) distributions within and between genets
- Sensitivity analyses using alternative ϕ and SH thresholds
- Evaluation of the proportion of within-genet and between-genet comparisons exceeding each threshold
- Identification of samples not assigned to multi-sample genets
- Comparison of singleton samples with individuals assigned to multi-sample genets
- Pairwise comparisons among singleton samples
- Verification of the final number of identified genets
- Generation of summary tables and supplementary results

These analyses were conducted to evaluate the consistency of genet assignments and examine the genetic differentiation between samples assigned to the same or different genets.

## R script

The main analysis script is:

`genet_assignment_phi_SH_sensitivity.R`

## Dataset

The analyses were conducted using:

- 421 trembling aspen samples
- 25 multi-sample genets comprising 413 samples
- 8 singleton genets
- 33 identified genets in total

Pairwise comparisons were used to evaluate genetic similarity within and between identified genets.

## Software and R packages

The analyses were performed using R and require the following packages:

- openxlsx
- dplyr
- tidyr

The packages can be installed in R using:

```r
install.packages(c(
  "openxlsx",
  "dplyr",
  "tidyr"
))
```

## Input files

The analyses require the following input files:

**Genet assignments**

`clonal_groups.xlsx`

This file contains the sample identifiers and their corresponding genet assignments.

The script uses the following columns:

- `Sample`: Sample identifier
- `Clone_Group`: Assigned genet

**Shared heterozygosity (SH)**

`SH_index_table.xlsx`

This file contains pairwise shared heterozygosity values.

The script uses the following columns:

- `X1`: First sample identifier
- `X2`: Second sample identifier
- `SH`: Shared heterozygosity index

**Pairwise relatedness (ϕ)**

`relatedness_output.xlsx`

This file contains pairwise relatedness estimates.

The script uses the following columns:

- `INDV1`: First sample identifier
- `INDV2`: Second sample identifier
- `RELATEDNESS_PHI`: Pairwise relatedness estimate (ϕ)

## Sensitivity analyses

### Pairwise relatedness (ϕ)

The following thresholds were evaluated:

- ϕ ≥ 0.350
- ϕ ≥ 0.375
- ϕ ≥ 0.400
- ϕ ≥ 0.425
- ϕ ≥ 0.450
- ϕ ≥ 0.475
- ϕ ≥ 0.490

For each threshold, the script calculates:

- Number of within-genet comparisons exceeding the threshold
- Percentage of within-genet comparisons exceeding the threshold
- Number of between-genet comparisons exceeding the threshold
- Percentage of between-genet comparisons exceeding the threshold

### Shared heterozygosity (SH)

The following thresholds were evaluated:

- SH ≥ 0.800
- SH ≥ 0.850
- SH ≥ 0.900
- SH ≥ 0.925
- SH ≥ 0.950
- SH ≥ 0.975
- SH ≥ 0.980

The same sensitivity calculations were performed for SH.

## Singleton genet verification

Eight samples were not assigned to any of the 25 multi-sample genets.

These samples were examined using:

- Maximum pairwise relatedness (ϕ) with samples assigned to multi-sample genets
- Maximum shared heterozygosity (SH) with samples assigned to multi-sample genets
- Pairwise comparisons among the eight singleton samples

The analyses evaluate whether these samples reach the selected high-similarity thresholds:

- ϕ ≥ 0.45
- SH ≥ 0.90

The results provide an additional assessment of the distinction between singleton samples and the identified multi-sample genets.

## Output files

The script generates the following Excel file:

`genet_assignment_phi_SH_sensitivity.xlsx`

The workbook contains six worksheets:

- `Sensitivity_table`: Comparison of alternative ϕ and SH thresholds
- `PHI_distribution`: Summary statistics for pairwise relatedness
- `SH_distribution`: Summary statistics for shared heterozygosity
- `Singleton_summary`: Maximum ϕ and SH values for each singleton sample
- `Singleton_pairs_PHI`: Pairwise relatedness among singleton samples
- `Singleton_pairs_SH`: Shared heterozygosity among singleton samples

## Reproducibility

The R code corresponds to the genet assignment sensitivity analyses described in the manuscript and Supporting Information.

The script uses previously generated pairwise relatedness estimates, shared heterozygosity values, and genet assignments.

Before running the script, the user must update the input directory specified by `base_dir` to match the location of the input files.

The script generates an Excel workbook containing the sensitivity analyses and singleton verification results.

## Code availability

The R code used for genet assignment sensitivity analyses based on pairwise relatedness (ϕ) and shared heterozygosity (SH) is publicly available in this repository.
