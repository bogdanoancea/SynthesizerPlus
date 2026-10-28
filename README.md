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

Install the current CRAN release:

```r
install.packages("SynthesizerPlus")
```

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

## Conditional synthesis of one target (experimental)

The copula represents each variable by one latent coordinate and all
dependence by correlations. When an analysis centres on **one response**
whose relation to the other columns is nonlinear, or on an **unordered
categorical** response, that response can instead be drawn from a
conditional model given the copula-generated predictors:

```r
# continuous response: CART (default), penalised GAM, or automatic choice
set.seed(1)
d <- data.frame(x = runif(600, -2, 2), z = rnorm(600))
d$y <- sin(3 * d$x) + 0.5 * d$z + rnorm(600, sd = 0.2 + 0.2 * abs(d$x))

fit_cart <- fit_synthesizer(d, nonlinear_target = "y", seed = 10)
fit_gam  <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "gam", seed = 10)
fit_auto <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto", seed = 10)
fit_auto$pooled$nonlinear$selection$scores   # cross-validated scores
syn <- generate(fit_auto, n = 1000, seed = 11)

# unordered categorical response: multinomial model with spline terms
fit_nom <- fit_synthesizer(iris, nominal_target = "Species", seed = 42)
syn_nom <- generate(fit_nom, n = 500, seed = 43)
```

* **Continuous response** (`nonlinear_target`, needs `rpart` or `mgcv`):
  a CART or GAM conditional mean plus resampled residuals, calibrated to the
  observed marginal of the response. `nonlinear_engine = "auto"` chooses
  the engine by training-only cross-validation of the conditional
  predictive distribution (CRPS).
* **Unordered categorical response** (`nominal_target`, needs `nnet`): a
  multinomial logit with natural-spline terms for numeric predictors (kept
  when they improve the AIC over linear terms), so category probabilities
  can depend nonlinearly on them.
* Only one target per fit; non-default `closeness`, `dependence` and
  `perturb` are not supported for such fits. Strata too small for the
  conditional model fall back to their copula, with a warning.
* The fitted conditional models contain training-data information and must
  not be released as if privacy-safe. Cross-validated or held-out accuracy
  does not certify synthetic-data quality or disclosure safety: use the
  evaluation and disclosure-risk tools, and
  `inst/examples/hybrid_nonlinear_benchmark.R` for paired comparisons on
  held-out data.

See `vignette("hybrid-nonlinear")` and `vignette("nominal-synthesis")`.

## Features

| Area | Functions |
|---|---|
| Synthesis | `fit_synthesizer()`, `generate()`, `synthesize()`, `simulate()` |
| Conditional synthesis (experimental) | `fit_synthesizer(nonlinear_target = , nonlinear_engine = "cart"/"gam"/"auto")`, `fit_synthesizer(nominal_target = )` |
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
* *Hybrid copula and conditional CART/GAM synthesis for nonlinear relationships* (experimental)
* *Nominal conditional synthesis* (experimental)

## Acknowledgement

The idea of combining empirical quantile functions with a Gaussian copula
follows the [`synthesizer`](https://github.com/markvanderloo/synthesizer)
package by Mark van der Loo. SynthesizerPlus is an independent
implementation.

## License

EUPL (European Union Public Licence)
