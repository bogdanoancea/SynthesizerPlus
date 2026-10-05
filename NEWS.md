# SynthesizerPlus 0.2.3

Fixes from a third code review (robustness and transparency; no change to
the synthesis model).

* Covariance and scale matrices passed to `r_mvnorm()`, `r_mvt()`,
  `r_mvlnorm()`, `r_mvskewnorm()` and `r_mvmixture()` are validated by one
  shared check: non-empty, square, numeric, finite, symmetric and positive
  semi-definite, each with a clear error message (previously some invalid
  inputs failed with internal errors such as "infinite or missing values in 'x'").
* Location and shape vectors (`mean`, `meanlog`, `xi`, skew-normal `alpha`,
  mixture means) must be finite; `Inf`/`NA` used to propagate silently into
  the output.
* `r_mvskewnorm()` requires a strictly positive diagonal of `omega`, as its
  parameterisation standardises by `sqrt(diag(omega))`.
* `r_dirichlet()` now draws the gamma variates on the log scale (with the
  `Y * U^(1/a)` boost for shapes below 1) and normalises with log-sum-exp, so
  very small concentration parameters no longer give `NaN` rows. Draws for
  `alpha < 1` therefore differ from 0.2.2 for the same seed.
* `make_corr(d = 1)` returns the 1 x 1 identity without consulting `rho`.
* t copula: when fewer than 10 complete rows are available, `df` is not
  estimated but fixed at 10; this now triggers a warning, is recorded as
  `df_estimated = FALSE` in the fitted copula and is shown by `summary()`.
* `?fit_synthesizer` gains an "Estimation details" section documenting the
  pairwise-complete estimation of latent correlations, the subsequent
  projection to a valid correlation matrix, and the grid/profile-likelihood
  estimate of the t-copula `df`.
* New distributional tests: Dirichlet means, variances and covariances;
  Kendall's tau of the Gaussian, t, Clayton, Gumbel and Frank copulas against
  their closed forms for several parameters and all pairs; uniform margins of
  every family; moments of the multivariate t, skew-normal, log-normal and
  normal-mixture generators.

# SynthesizerPlus 0.2.2

Fixes from a second code review.

* Argument validation: `r_dirichlet()` requires finite positive `alpha`;
  `make_corr()` requires finite numeric `rho`; `disclosure_risk()` requires a
  finite non-negative `tolerance`; `seed` arguments must be finite and in
  the integer range.
* `r_mvmixture()` validates every component (finite means, finite and
  positive semi-definite covariance matrices) before sampling, so an
  invalid component is reported even when its weight is zero.
* Row keys treat `NaN` like `NA` (both are missing values).
* `dcr()`: `exact_match_rate` now counts only records identical to a real
  record on *all* variables. New `zero_distance_rate` (distance zero on the
  variables compared) and `compared_share` (share of variables on which the
  nearest real record was compared) make the meaning explicit under
  `na = "exclude"`. `real_duplicate_rate` is defined in the same way.
* `compare_synthetic()` flags a discriminator AUC well below 0.5, which
  indicates synthetic records that copy real ones. The AUC is deliberately
  not converted to a symmetric measure such as `max(AUC, 1 - AUC)`, which
  would report copies of the real data as easy to distinguish.

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
