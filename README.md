# CompRisksBenchmaRk <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/VallejosGroup/CompRisksBenchmark/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/VallejosGroup/CompRisksBenchmark/actions)
<!-- badges: end -->

A unified R framework for **fitting, predicting, and benchmarking competing risks survival models** using nested cross-validation.

---

## Overview

`CompRisksBenchmaRk` provides:

- **A model registry** — register any competing risks model once and reuse it across datasets and workflows.
- **Built-in model wrappers** — Fine-Gray (FGR), penalised Fine-Gray (FGRP), cause-specific Cox PH (csCPH), and random survival forests (RSF).
- **Nested cross-validation engine** — unbiased hyperparameter selection and performance estimation from pre-built split directories.
- **Comprehensive metrics** — time-varying Brier score, integrated Brier score, AUC, C-index (pec & SurvMetrics), calibration plots (pseudo-value + LOESS), and observed/expected ratios.
- **Clinical utility** — decision curve analysis with competing risks.
- **Missing data** — median/mode and MICE imputation integrated into the CV pipeline.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("VallejosGroup/CompRisksBenchmark")
```

### Required dependencies

```r
install.packages(c(
  "survival", "riskRegression", "prodlim", "cmprsk",
  "ggplot2", "dplyr", "tidyr", "scales",
  "arrow", "jsonlite", "rsample", "tibble",
  "data.table", "foreach", "abind"
))
```

### Optional dependencies (for specific models / features)

| Package | Used by |
|---|---|
| `randomForestSRC` | RSF model |
| `fastcmprsk` | FGRP model |
| `pec` | C-index (pec) metric |
| `mice` | MICE imputation |

---

## Quick start

### 1 — Simulate data and create nested splits

```r
library(CompRisksBenchmaRk)

# Simulate a competing risks dataset (2 causes)
df <- sim_cmprks(n = 2000, seed = 42)
head(df)

# Evaluation time grid
times <- quantile(df$Tobs[df$cause != 0], probs = seq(0.1, 0.9, by = 0.1))

# Write stratified nested CV splits to disk
create_nested_splits(
  df          = df,
  time_col    = "Tobs",
  event_col   = "cause",
  out_dir     = "BenchResults",
  outer_folds = 5,
  inner_folds = 3,
  seed        = 123,
  times       = times
)
```

### 2 — Run nested CV for a model without tuning

```r
results_fgr <- nested_cv_from_bench(
  out_dir   = "BenchResults",
  model_key = "FGR",
  seed      = 123,
  verbose   = TRUE
)
```

### 3 — Run nested CV for a model with hyperparameter tuning

```r
results_rsf <- nested_cv_from_bench(
  out_dir   = "BenchResults",
  model_key = "RSF",
  grid      = list(
    ntree     = c(200L, 500L),
    nodesize  = c(5L, 15L),
    splitrule = "logrankCR"
  ),
  seed    = 123,
  verbose = TRUE
)
```

### 4 — Summarise and plot results

```r
all_results <- list(FGR = results_fgr, RSF = results_rsf)
times_vec   <- results_fgr[[1]]$times

tab <- summarize_out_of_sample(all_results, num_causes = 2,
                                times = times_vec)
agg <- aggregate_out_of_sample(tab, metrics = c("bs", "ibs", "auc"))
plot_out_of_sample(agg, metrics = c("bs", "auc"), cause = 1)
```

### 5 — Register a custom model

```r
register_cr_model(
  key = "my_model",
  fit = function(x, y_time, y_event, args = list()) {
    # ... train your model ...
    list(my_fit = ..., causes = sort(unique(y_event[y_event != 0])))
  },
  predict_cif = function(fit_obj, x_new, times) {
    # return array [n, K, Tm]
  },
  info = function() list(
    name         = "My custom model",
    supports     = "CIF",
    needs_tuning = FALSE,
    default_grid = function() tibble::tibble()
  )
)

list_cr_models()   # check it's registered
```

---

## Package structure

```
CompRisksBenchmaRk/
├── R/
│   ├── CompRisksBenchmaRk-package.R   # package-level documentation
│   ├── registry.R                     # model registry
│   ├── utils.R                        # data prep, imputation, CV helpers
│   ├── metrics.R                      # Brier, AUC, calibration, DCA
│   ├── nested_cv.R                    # nested CV engine
│   ├── model_FGR.R                    # Fine-Gray wrapper
│   ├── model_FGRP.R                   # Penalised Fine-Gray wrapper
│   ├── model_csCPH.R                  # Cause-specific Cox PH wrapper
│   └── model_RSF.R                    # Random survival forest wrapper
├── tests/
│   └── testthat/
│       └── test-core.R
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

## Citation

If you use this package, please cite the original benchmark study:

> VallejosGroup (2024). *CompRisksBenchmark: a benchmarking framework for competing risks models*. GitHub. https://github.com/VallejosGroup/CompRisksBenchmark

---

## Licence

MIT © VallejosGroup
