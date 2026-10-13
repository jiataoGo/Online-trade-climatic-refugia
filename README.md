# Plant trade and climatic refugia: GLMM analysis and figures

R scripts and input files for fitting the GLMMs underlying Figures 2 and 3, selecting models, performing statistical tests, and reproducing Figures 1–3 and Supplementary Figures S1–S2.

## Repository structure

- `code/`: the complete GLMM analysis and one R script per figure.
- `data/`: trade records, map layers, species suitability rasters, and cached Figure S2 data.
- `models/`: the full GLMM, candidate-model selection table, and selected GLMM, saved as RDS objects.
- `results/`: figure outputs in PNG and PDF format, with supporting tables.

The trade dataset contains 21,351 records representing 748 taxa: 472 non-naturalized, 232 naturalized non-invasive, and 44 invasive taxa. The raster collection contains 757 taxa; Figure S2 selects the 748 taxa present in the trade dataset.

## Requirements

The figures were developed using R 4.4.1. Install the required packages in R:

```r
install.packages(c(
  "tidyverse", "sf", "terra", "ggsci", "ggpubr", "patchwork",
  "lme4", "MuMIn", "performance", "car", "ggeffects", "ragg"
))
```

PDF export requires an R installation with Cairo support. Figure S2 uses Arial; font availability can affect appearance.

## Rebuild the GLMM analysis from raw records

Set the working directory to the repository root. Run:

```r
source("code/GLMM_analysis.R")
source("code/Figure_2.R")
source("code/Figure_3.R")
```

`GLMM_analysis.R` reads the trade CSV, aggregates records into 3,333 groups, constructs the original orthogonal contrasts, and standardizes latitude and longitude. It fits a binomial GLMM with a logit link and a shop-level random intercept, then fits the hierarchical candidate models and ranks them by AICc. It uses the original `bobyqa` optimizer and iteration limit, with up to six workers for model selection.

This script does not read precomputed RDS files. It saves newly fitted full and selected models and the candidate-model table to `models/`. It also writes model data, coefficients, diagnostics, performance metrics, likelihood-ratio tests, and R session information to `results/`.

The supplied RDS files are the original saved analysis results. The reconstructed model input was checked against the original model data and matched, and the candidate set contains the same 803 models. A full rerun of model selection was interrupted before completion; end-to-end numerical equivalence of the new workflow has not yet been verified.

Each selected interaction is tested by dropping it from the selected model. Main effects are tested by dropping them individually from an additive model, following the original analysis. The exported `AIC_full` and `AIC_reduced` columns identify the corresponding models explicitly. Model selection uses AICc; these test-table columns report ordinary AIC.

## Reproduce figures from saved models

Set the working directory to the repository root, which contains `code`, `data`, `models`, and `results`. Run any script independently in R:

```r
source("code/Figure_1.R")
source("code/Figure_2.R")
source("code/Figure_3.R")
source("code/Figure_S1.R")
source("code/Figure_S2.R")
```

Each script writes its outputs to `results/`, replacing existing files with the same names. Figures 2 and 3 read the selected model in `models/`: run the GLMM workflow above first to use a fresh fit, or use the supplied RDS to reproduce plots without repeating model selection. RDS files store R objects, not the code that created them; the fitting code is provided separately in `GLMM_analysis.R`.

Figure S2 reads `data/plot_data.rds` by default and checks that the trade CSV matches the cached data. To rebuild its summaries from the supplied suitability rasters, run the following command from the repository root in a terminal:

```sh
Rscript code/Figure_S2.R --recalculate
```

This also refreshes the cache and writes summary tables and richness rasters to `results/`.
