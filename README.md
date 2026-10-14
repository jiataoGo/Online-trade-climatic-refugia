# Plant trade and climatic refugia

`code/`: six R scripts; `data/`: inputs; `models/`: saved models and selection results; `results/`: Figures 1-3 and S1-S2 (PNG/PDF).

Use R 4.4.1 with Cairo support. Set the working directory to this folder.

```r
install.packages(c("tidyverse", "sf", "terra", "ggsci", "ggpubr", "patchwork",
                   "lme4", "MuMIn", "performance", "car", "ggeffects", "ragg"))

# Draw all five figures using the supplied models and S2 cache.
for (f in c("1", "2", "3", "S1", "S2")) {
  source(paste0("code/Figure_", f, ".R"))
}
```

To refit models from the trade data, run `source("code/GLMM_analysis.R")` before drawing Figures 2-3; this replaces the saved models. To rebuild S2 from source rasters, run `Rscript code/Figure_S2.R --recalculate` in a terminal.

Figures 1-3 and S1 use Times New Roman (`serif`); S2 uses Arial.
