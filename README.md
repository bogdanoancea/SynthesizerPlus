# SynthesizerPlus

[![DOI](https://zenodo.org/badge/1405547507.svg)](https://doi.org/10.5281/zenodo.23193076)

Synthetic data generation for R with copulas, built-in quality and
disclosure-risk evaluation, ggplot2 visualisation and multi-format I/O.

*SynthesizerPlus* learns the marginal distribution of every column and the
dependence between columns (Gaussian or Student-*t* copula) from a real data
set, and generates new records with the same statistical structure. It
handles numeric, integer, categorical, ordinal, logical, date and date-time
columns, missing values and their patterns, stratified models, and
univariate or multivariate time series.

## Installation

```r
# development version
remotes::install_github("bogdanoancea/SynthesizerPlus", build_vignettes = TRUE)
```

Optional packages enable extra file formats: `arrow` or `nanoparquet`
(Parquet), `arrow` (Feather), `hdf5r` (HDF5), `jsonlite` (JSON/NDJSON),
`readxl` with `writexl` or `openxlsx` (Excel), `haven` (SPSS, Stata, SAS),
`fst`, and `data.table` (fast CSV).

## Quick start

```r
library(SynthesizerPlus)

fit <- fit_synthesizer(airquality, by = "Month")   # learn
syn <- generate(fit, n = 1000, seed = 1)           # generate

cmp <- compare_synthetic(airquality, syn, seed = 1) # evaluate
cmp
plot(cmp, type = "marginals")
plot(cmp, type = "association")

write_data(syn, "synthetic.parquet")               # deliver
```

## Features

| Area | Functions |
|---|---|
| Synthesis | `fit_synthesizer()`, `generate()`, `synthesize()`, `simulate()` |
| Disclosure control | `disclosure_risk()` (replicated uniques, correct attribution probability), `dcr()`, simulated `census_sim` data |
| Closeness control | `generate(closeness = )`, `perturbation()`, `calibrate_closeness()` |
| ML data augmentation | `augment_data()`, `ts_windows()` |
| Time series | `fit_synthesizer()` on `ts`: copula-VAR, block and stationary bootstrap |
| Large data | `generate_to_file()` (chunked, CSV / Parquet / HDF5 / ...) |
| Multivariate distributions | `r_mvnorm()`, `r_mvt()`, `r_mvlnorm()`, `r_mvskewnorm()`, `r_dirichlet()`, `r_mvmixture()`, `r_copula()` (Gaussian, t, Clayton, Gumbel, Frank), `r_mvdist()`, `make_corr()` |
| Evaluation | `compare_synthetic()`, `marginal_metrics()`, `association_matrix()`, `pmse()`, `discriminator_auc()`, `ci_overlap()`, `dcr()` |
| Visualisation | `plot_marginals()`, `plot_categorical()`, `plot_association()`, `plot_pairs()`, `plot_qq()`, `plot_ts()`, `plot_acf()`, `plot_dcr()`, `autoplot()` |
| I/O | `read_data()`, `write_data()`, `supported_formats()`, `match_types()` |

Supported file formats: CSV, TSV (optionally `.gz`/`.bz2`/`.xz`), RDS,
RData, Parquet, Feather/Arrow IPC, HDF5, JSON, NDJSON, Excel, SPSS, Stata,
SAS (`sas7bdat` read, `xpt` read/write) and fst.

## Example scripts

Runnable scripts are installed with the package:

```r
ex <- function(f) system.file("examples", f, package = "SynthesizerPlus")
source(ex("demo.R"), echo = TRUE)          # a tour of the whole package
source(ex("closeness.R"), echo = TRUE)     # closeness, perturbations, calibration
source(ex("privacy.R"), echo = TRUE)       # releasing confidential microdata
source(ex("augmentation.R"), echo = TRUE)  # more training data for ML models
```

## Vignettes

* *Getting started* — the model, column types, missing data, strata, dependence control
* *Evaluating and visualising synthetic data* — fidelity, utility and disclosure risk
* *Releasing confidential microdata as synthetic data*
* *Controlling how close synthetic data are to the real data*
* *Augmenting small training sets for machine learning*
* *Regularising high-variance models with synthetic data*
* *Synthetic time series*
* *Generating data from multivariate distributions*
* *Reading and writing data: CSV, Parquet, HDF5 and more*

## Acknowledgement

The idea of combining empirical quantile functions with a Gaussian copula
follows the [`synthesizer`](https://github.com/markvanderloo/synthesizer)
package by Mark van der Loo. SynthesizerPlus is an independent
implementation.

## License

EUPL (European Union Public Licence)
