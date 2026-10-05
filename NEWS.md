# SynthesizerPlus 0.2.1

Fixes from an external code review.

* `r_copula()` and `r_mvdist()` now require `corr` to be a correlation
  matrix (square, symmetric, unit diagonal, entries in [-1, 1], positive
  semi-definite). Previously a covariance matrix was accepted and produced
  margins that were not uniform. The *t* copula and `r_mvt()` validate
  `df` (positive; `Inf` gives the Gaussian copula); Archimedean `theta`
  must be finite.
* Row keys used for strata, `disclosure_risk()` and exact-copy detection
  are now built from integer codes: a missing value can no longer collide
  with the string `"<NA>"`, and values containing the old separator can no
  longer make different rows look identical. Strata are listed in order of
  first appearance.
* `dcr()` gains `na = c("category", "exclude")`. The default treats
  missingness as a category (as before, now documented); `"exclude"` gives
  the standard Gower coefficient, ignoring variables missing in either
  record.
* `discriminator_auc()` learns all preprocessing (imputation, scaling,
  retained categories) within the training folds, so the cross-validated
  AUC has no leakage from the held-out fold. Its documentation explains
  that values clearly below 0.5 indicate copies of real records.
* `r_mvmixture()` rejects all-zero or non-finite weights and components of
  inconsistent dimension with clear messages; `margin_categorical()` rejects
  duplicated levels and non-finite probabilities; the quantile functions it
  and `margin_empirical()` return check that probabilities lie in [0, 1].
* `match_types()` reads date-time strings without zone information as
  local times in the template's time zone (strings ending in `Z` are UTC,
  explicit offsets such as `+02:00` are honoured), and accepts mixed
  formats within a column.
* Documentation: continuous marginals are described as *interpolated*
  empirical quantile functions; the `noise` perturbation preserves mean and
  variance approximately rather than exactly.

# SynthesizerPlus 0.2.0

Complete rewrite. the API changed and is not backward compatible with 0.1.0.

## New features

* `fit_synthesizer()` / `generate()` / `synthesize()` with Gaussian, Student-t
  and independence copulas; categorical, ordinal and logical variables take
  part in the dependence structure through the distributional transform.
* Joint missing-data model, stratified models (`by`) with pooled fallback,
  identifier columns, compressed quantile functions (`knots`), dependence
  strength control and `simulate()` for multiple synthetic data sets.
* Control of how close synthetic data are to the real data: `closeness`
  argument of `generate()`, `synthesize()` and `augment_data()`;
  `perturbation()` for record-level (noise, swap) and distribution-level
  (shift, scale, temperature, dependence up to 2) settings; and
  `calibrate_closeness()` to reach a target discriminator AUC, DCR ratio,
  pMSE ratio or marginal distance.
* `discriminator_auc()` now uses squares and pairwise interactions of the
  numeric variables by default (`model = "quadratic"`, ridge-penalised), so
  that differences in spread and correlation are detected; `pmse()` gains the
  same option.
* `disclosure_risk()` for statistical disclosure control: replicated uniques
  and correct attribution probability (with no-data baseline) for chosen key
  and sensitive variables. `dcr()` gains `close_share`, the share of
  near-copies of real records. `calibrate_closeness()` accepts any custom
  metric, e.g. a risk threshold.
* New simulated census-like data set `census_sim` and vignette "Releasing
  confidential microdata as synthetic data".
* `augment_data()` and `ts_windows()` for enlarging small training sets of
  tabular data and time series for machine learning.
* Time series: copula-VAR (with optional fixed order), moving-block and
  stationary bootstrap, iid.
* Multivariate distributions: normal, t, log-normal, skew-normal, Dirichlet,
  normal mixtures; Gaussian, t, Clayton, Gumbel and Frank copulas;
  `r_mvdist()` for arbitrary margins.
* Evaluation: marginal distances, mixed-type association matrix, pMSE,
  cross-validated discriminator AUC, confidence-interval overlap, distance to
  closest record; `compare_synthetic()` bundles them.
* ggplot2 visualisations for marginals, associations, pairs, Q-Q, time
  series, autocorrelation and disclosure risk.
* `read_data()` / `write_data()` for CSV, TSV, RDS, RData, Parquet, Feather,
  HDF5, JSON, NDJSON, Excel, SPSS, Stata, SAS and fst; `generate_to_file()`
  for chunked generation of very large data sets.

## Robustness

* Arguments that are misspelled or passed to the wrong function (e.g.
  `closeness` given to `fit_synthesizer()`) now trigger a warning instead of
  being silently ignored.
* `closeness`, `perturb` and `dependence` are not available for time-series
  synthesizers; passing them now warns instead of being silently ignored.
* Small samples of continuous variables are no longer treated as discrete:
  a numeric variable is discrete only if it has at most
  `discrete_threshold` values *and* these are whole numbers or repeat.
* Numeric perturbations (`noise`, `shift`, `scale`) also apply to numeric
  variables modelled as discrete, e.g. within strata; whole-number variables
  are rounded stochastically, so shifts are not biased by rounding.
* Numeric columns with a class, such as `difftime`, keep their class and
  units; 64-bit integers (`integer64`) are rejected with a clear message.
* `disclosure_risk()` gains `ignore` to leave out direct identifiers, which
  otherwise mask exact copies.
* `augment_data()` accepts plain vectors.
* Clearer errors for perturbing stratification or identifier variables and
  for non-positive-definite correlation matrices with `pd_method = "none"`.

## Bugs fixed relative to 0.1.0

* Categorical columns were sampled independently of the copula.
* Factor synthesizers forced every level to appear at least once.
* Missing values were silently dropped.
* `rankcor < 1` changed the marginal distributions.
* `make_synthesizer.POSIXct()` returned a character vector.
* `seed = ` left a global `.Random.seed` behind.
* Time series could only be generated at their original length, and
  `synthesize()` failed for series longer than 1000 observations.
