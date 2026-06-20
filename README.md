# cmice — Constrained Multiple Imputation by Chained Equations

[![R-CMD-check](https://github.com/nguyenthephuong1/cmice/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/nguyenthephuong1/cmice/actions/workflows/R-CMD-check.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

R implementation of **CMICE** (Constrained MICE), a generalisable modification to the Fully Conditional Specification (FCS) multiple imputation framework that enforces observation-specific structural constraints during the multinomial sampling step. Designed for cancer-registry histology imputation where nonspecific morphology codes (e.g., ICD-O-3 8010 "carcinoma, NOS"; 8046 "non-small-cell carcinoma") restrict the admissible specific subtypes by diagnostic classification rules.

Compatible with the standard `mice()` / `with()` / `pool()` pipeline.

## Install

```r
# install.packages("remotes")
remotes::install_github("nguyenthephuong1/cmice")
```

## Quick start

> The example below uses the **mice** package (a suggested, not required,
> dependency of cmice) for the imputation pipeline, and dplyr for the
> illustrative data preparation. Install them first if needed:
> `install.packages(c("mice", "dplyr"))`.

```r
library(mice)
library(cmice)

# 1. Decompose the source variable by NA pattern (one column per code)
df_imp <- df %>%
  mutate(
    hist1 = case_when(hist_cat == "8010" ~ NA_character_,
                      hist_cat == "8046" ~ "NA2",
                      TRUE               ~ as.character(hist_cat)),
    hist2 = case_when(hist_cat == "8046" ~ NA_character_,
                      hist_cat == "8010" ~ "NA1",
                      TRUE               ~ as.character(hist_cat)),
    hist1 = factor(hist1),
    hist2 = factor(hist2)
  )

# 2. Define constrained imputation rules (named-category interface,
#    recommended; emits a clear error if a label is misspelled).
mice.impute.cons_na1 <- function(y, ry, x, wy = NULL, ...) {
  cmice::mice.impute.cons_named(y, ry, x, wy, exclude = c("SMC", "NA2"), ...)
}
mice.impute.cons_na2 <- function(y, ry, x, wy = NULL, ...) {
  cmice::mice.impute.cons_named(y, ry, x, wy,
                                exclude = c("SMC", "LAC", "NA1"), ...)
}

# 3. Run mice() with these custom methods
imp <- mice(df_imp,
            method = c(hist1 = "cons_na1", hist2 = "cons_na2"),
            m = 20, maxit = 5, seed = 166)
```

## Why named-category interface?

The legacy positional interface (`constrain = c(6, "NA2")`) silently coerces all elements to character via R's `c()` rules and fails silently if a positional index is wrong. The named interface validates labels against `levels(y)` at call time and raises an informative error on misspelling:

```r
mice.impute.cons_named(y, ry, x, exclude = c("SMc"))
# Error: `exclude` contains labels not in levels(y): "SMc".
# Available levels: "ADC", "SQC", "LAC", "SMC".
```

## Citation

If you use `cmice` in published work, please cite:

> Nguyen, P. T. (2026). Constrained Multiple Imputation for Nonspecific Histological Type in Cancer Registries Data: Methods and Practical Guidance. *Statistics in Medicine* (under review).

The R package itself is archived at Zenodo: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) (placeholder; resolved DOI after first archived release).

## License

MIT (see LICENSE file).
