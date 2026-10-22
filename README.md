# pccf_exposure_comparison

Code, aggregate results and a portable R function for a national comparison of
non-random, population-weighted exposure assignment from residential postal codes with
PCCF+ random allocation, using 12.4 million addresses in the National Address Register as
the reference.

> Peters PA, McGaughey T. Population-Weighted Exposure Assignment from Residential Postal
> Codes: a National Comparison of Methods. *Preprints* 2026 (posted 22 October 2026).
> https://doi.org/10.20944/preprints202610.2150.v1

## Description

Statistics Canada's Postal Code Conversion File Plus (PCCF+) assigns a residential postal
code with several candidate locations to one location drawn at random in proportion to
population. The study compared four assignments for 12,410,638 residential National Address
Register (NAR) addresses:

- **A**: the population-weighted mean of the exposure over all PCCF+ candidate locations of
  the postal code;
- **B**: the exposure at the population-weighted centroid of the candidates;
- **C**: PCCF+ Version 8D random allocation; and
- **D**: the exposure at the address (reference).

Exposures were growing-season greenness (MODIS Enhanced Vegetation Index, EVI) and quintiles
of the 2021 Material and Social Deprivation Index (MSDI). The weighted mean had the same
greenness bias as random allocation and a 13.5% lower mean absolute error (20% lower for
rural postal codes). The modal candidate quintile reduced material deprivation quintile
misclassification from 12.2% to 10.2% overall and from 46.2% to 37.0% for rural postal
codes. The weighted centroid had the smallest positional error but did not improve
greenness.

## Contents

| Folder | Content |
|---|---|
| `R/` | **Portable function** `pcweight.R` (methods A and B for any dataset with postal codes) and a worked example, `pcweight_example.R` |
| `scripts/` | Analysis scripts `sw1` to `sw8`, which reproduce every table and figure in the paper |
| `outputs/` | Aggregate results reported in the paper: manuscript and supplementary tables, figures, and the step outputs they are built from |
| `paper/` | The preprint as posted on Preprints.org (PDF) |
| `data/` | Output folder of the scripts (not tracked) |
| `config.R` | All file paths, relative to the repository root and the validation checkout |

## Using the portable function

`R/pcweight.R` assigns exposures to six-character postal codes by methods A and B, using
the candidate locations and weights of a licensed PCCF+ 8D installation. Inputs are a table
of records with postal codes (optionally with the `Link_Source` of a PCCF+ run), an exposure
layer as a `terra` raster (value of the containing cell) or `sf` points (value of the nearest
point), and optionally DA attributes such as MSDI scores and quintiles and the DA boundary
file. It returns one row per postal code with the link source, the number and effective
number of candidates, the weighted mean, standard deviation and share missing of each
exposure, the weighted mean score, its quintile, the quintile probabilities and the modal
quintile, and the weighted centroid with its exposure and DA values.

```r
source("R/pcweight.R")
pccf <- pcw_read_pccf("<PCCF+ 8D folder>/DATA_R")
res  <- pcweight(records, pccf, pcode = "PCODE", exposure = evi_raster,
                 da_attr = msdi, da_scores = c(mat = "SCOREMAT", soc = "SCORESOC"),
                 da_quintiles = c(mat = "QuintMat", soc = "QuintSoc"), da_boundaries = da)
```

See `R/pcweight_example.R`. The function needs R with `data.table` and `sf` (and `terra`
for raster exposures). When no PCCF+ run is available, the link source is inferred with the
matching sequence of the PCCF+ 8D residential routine; this agreed with PCCF+ for all
754,656 postal codes in the study. The function reproduced the analysis pipeline exactly for
every postal code of Prince Edward Island (`scripts/sw7_pcweight_check.R`).

## Reproducing the paper

The analysis extends the PCCF+/NAR validation study and reuses its data and outputs:

1. **Run the validation pipeline.** Clone
   [pccf_nar_validation](https://github.com/unitepangaea/pccf_nar_validation) next to this
   repository (or set `VALIDATION` in `config_local.R`) and run its pipeline. That
   repository explains how to obtain the public data and where to place the licensed PCCF+
   files.
2. **Install R packages** (R 4.6 or later): `Rscript install_packages.R`
3. **Run the analysis**: `Rscript run_all.R`, or run the scripts one at a time in the order
   listed there. All scripts are run from the repository root, write to `data/derived/` and
   stop with an error if any of their checks fails. Each script can first be run for one
   province, e.g. `SW_PROV=11 Rscript scripts/sw1_candidates.R`.

Requirements: about 5 GB of free disk space and 16 GB of memory. To change paths, copy
`config_local_example.R` to `config_local.R`.

| Script | Step | Paper |
|---|---|---|
| `sw1_candidates.R` | Candidate locations and weights for every postal code | 2.3, Table 1 |
| `sw2_candidate_exposures.R` | EVI and MSDI at candidates; weighted centroids and their DA | 2.4 |
| `sw3_postal_values.R` | Methods A and B per postal code (both weight sets) | 2.4 |
| `sw4_compare.R` | Comparison with the address reference, Bayesian bootstrap, cell-level models | 2.5–2.6, 3.1–3.5 |
| `sw5_sensitivity.R` | Modal quintile rule and sensitivity analyses 1–6 | 2.7, 3.5–3.6 |
| `sw6_figures_tables.R` | Tables, Figures 1–4, Supplementary Tables S1–S12 | 3 |
| `sw7_pcweight_check.R` | Checks of `R/pcweight.R` against the pipeline; inferred link sources | 2.8, 3.7 |
| `sw8_manuscript_tables.R` | Manuscript tables as Markdown | 3 |

The versions of the outputs used in the paper are in `outputs/`. The bootstrap uses fixed
random seeds.

## Data and licences

- Code: MIT licence (`LICENSE`).
- Aggregate tables and figures in `outputs/`: CC BY 4.0 (`LICENSE-DATA.md`).
- Input data remain under their own terms (see the validation repository). PCCF+ is used
  under licence from Statistics Canada; neither the PCCF+ files nor any address-,
  postal-code- or DA-level result derived from them is published here.

## Citation

See `CITATION.cff`, or cite the preprint above (https://doi.org/10.20944/preprints202610.2150.v1).
The posted preprint is included as [`paper/preprints202610.2150.v1.pdf`](paper/preprints202610.2150.v1.pdf).
The study builds on the validation of PCCF+ against the NAR
(https://doi.org/10.20944/preprints202610.1204.v1; code at
https://github.com/unitepangaea/pccf_nar_validation).
