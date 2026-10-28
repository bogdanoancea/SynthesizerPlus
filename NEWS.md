# SynthesizerPlus 0.4.0

## Change in results for unordered categorical columns

* Unordered factors and character columns are now placed on their latent
  copula coordinate in a data-driven order (the first principal component of
  each category's mean normal scores on the numeric, ordered and logical
  columns), instead of the order of their labels. Synthetic data, fitted
  correlations and closeness calibrations involving such columns therefore
  differ from those of earlier versions with the same seed. Renaming or
  reordering the categories no longer changes the fitted model, and
  associations between a nominal column and the other columns are better
  preserved. Ordered factors and logical columns keep their natural order;
  results for data without unordered categorical columns are unchanged.

## Conditional synthesis of one continuous response

* `nonlinear_engine = "gam"` adds a penalised GAM (`mgcv`) alternative to the
  CART conditional model, with optional tensor interactions
  (`nonlinear_gam_interactions`) and local residual resampling
  (`nonlinear_local_k`).
* `nonlinear_engine = "auto"` selects CART, additive GAM or GAM with
  interactions by training-only K-fold cross-validation
  (`nonlinear_cv_folds`, `nonlinear_cv_rule`). Candidates are scored by the
  continuous ranked probability score of their conditional predictive
  distribution (`nonlinear_cv_score = "crps"`, default), which accounts for
  the conditional spread as well as the conditional mean;
  `nonlinear_cv_score = "mse"` scores the conditional mean only. Both scores
  are stored in `fit$pooled$nonlinear$selection`.
* `nonlinear_residual` (`"leaf"`, `"global"`, `"none"`) and
  `nonlinear_calibration` (`"quantile"`, `"none"`) control the residual
  generator and the marginal calibration of the response.
* The quantile calibration now maps each generated value through a reference
  distribution of conditional draws stored at fitting time, instead of
  ranking the values generated in the same call. Results no longer depend on
  the number of rows generated per call, chunk or stratum (previously a
  single generated row always received the median response).
* The CART model now encodes predictor categories from all training rows, so
  a category that occurs only where the response is missing no longer stops
  generation with an "unseen categorical predictor" error.
* When a stratum has its own copula but too few observed responses for the
  conditional model, `fit_synthesizer()` now warns and lists such strata in
  `fit$strata$conditional_fallback`; `print()` reports them.

## Nominal categorical response (experimental)

* `nominal_target` generates one unordered factor, character or logical
  response from a multinomial model (`nnet`) of its conditional
  distribution given the copula-generated predictors. Numeric predictors
  enter through natural cubic spline bases (`nominal_spline_df`, default 3),
  so category probabilities can depend nonlinearly on them; the linear model
  is fitted too and the one with the smaller AIC is used (`0` fits the
  linear model only). The response is omitted from the predictor copula; its missingness
  indicator stays in the joint missingness model. Mutually exclusive with
  `nonlinear_target`. Strata without a converged model fall back to the
  copula, with the same warning and record as above.
* `nominal_decay` and `nominal_maxit` (default 1000 iterations) are passed
  to `nnet::multinom()`.
* `splines` (a base R package) is now imported.

## Student-t copula fitting

* The degrees-of-freedom profile log-likelihood is evaluated through a
  Cholesky factor with common normalising terms cancelled. It is
  algebraically equivalent to the previous likelihood and gives the same
  estimates.
* `t_profile = TRUE` records elapsed times of the fitting stages under
  `fit$pooled$copula$profile`.

## Tests

* New tests cover the category order of unordered columns, the GAM and
  automatic engines, the residual and calibration options, batch-size
  independence of the calibration, categories seen only with a missing
  response, conditional-model fallback in strata, CRPS selection, nonlinear
  nominal dependence, and the t-copula likelihood.
* The synthetic missingness-rate test now allows for binomial Monte Carlo
  variation.

# SynthesizerPlus 0.3.0

## Experimental automatic conditional engine selection

* Added optional `nonlinear_engine = "auto"` with training-only K-fold
  selection among CART, additive GAM and GAM with numeric tensor
  interactions. The default engine remains CART.
* Added `nonlinear_cv_folds` (default 5) and `nonlinear_cv_rule`
  (`"one_se"` or `"min"`). Model-selection diagnostics, including per-fold
  MSE and candidate failures, are saved in the fitted nonlinear model.
* Included regression tests. CV predictive loss is not a guarantee of
  synthetic-data utility or privacy.


## Experimental GAM conditional engine

* Added optional `nonlinear_engine = "gam"` for one continuous hybrid target.
  Numeric predictors enter through penalised cubic regression splines and
  categorical predictors through factor effects. Tensor interactions may be
  enabled explicitly with `nonlinear_gam_interactions = TRUE`.
* For GAM, the existing `nonlinear_residual = "leaf"` option uses nearest
  fitted-mean residual neighbours rather than tree leaves. CART, Gaussian
  copula, and Student-t copula defaults and results are unaffected.
* Added GAM tests. GAM fitting requires suggested package `mgcv`. This is experimental and no
  superiority or privacy guarantee is claimed.


## Controlled hybrid ablation (experimental)

* The optional hybrid model now supports `nonlinear_residual = "leaf"`,
  `"global"`, or `"none"`, and `nonlinear_calibration = "quantile"` or
  `"none"`. Defaults reproduce the original hybrid method. The six
  combinations allow controlled study of stochastic residual resampling and
  response marginal calibration without changing the public default.
* Added reproducibility and input-validation tests. No superiority is claimed
  before external replicated comparisons.

## Experimental extension (source candidate)

* Added optional `nonlinear_target` argument to `fit_synthesizer()` to generate
  one continuous response from a CART conditional model while retaining
  copula-based generation of predictors. Leaf-grouped residual resampling
  and empirical quantile calibration are used to approximate the response
  conditional distribution and preserve its marginal. The default copula
  path is unchanged.
* Added regression tests and an exploratory held-out nonlinear prediction
  benchmark script. This feature requires the suggested `rpart` package.
* The hybrid's final joint distribution is not an elliptical copula;
  non-default closeness/dependence controls are not supported. Fitted trees
  and residual pools contain training data information and are not
  disclosure-safe model artifacts.


## Bug fixes

* Fixed a dependence on the declared `levels()` order of unordered factors
  during marginal fitting and latent-copula synthesis. Permuting the levels
  without changing observed labels no longer changes synthetic observations
  when the fitting and generation seeds are held fixed.
* Unordered factors now use a deterministic (radix-sorted) internal order of
  observed category labels. Explicit ordering of `ordered` factors is
  unchanged. Generated factors retain their original class and level
  metadata, including unused levels.

## Tests

* Added regression tests for unordered-factor level permutation invariance,
  preservation of original output factor metadata, and ordered-factor
  behaviour.

## Compatibility and limitations

* No changes to the public function signatures. As the internal category
  encoding has changed, fitted models involving unordered factors may give
  different synthetic results from earlier package versions. Refit models
  when moving from 0.2.x to 0.3.0.
* The fix does not make one-dimensional latent intervals a fully nominal
  statistical model. Changing the *labels* of nominal categories (rather
  than reordering the same labels) can still change the fitted model.


## Testing

* Added regression tests verifying that reordering the levels
  of unordered factors does not change synthetic observations
  when random seeds and all other settings are identical.

* Added checks for preservation of factor classes, original
  level metadata, and ordered-factor behaviour.

## Compatibility

* Existing public function signatures remain unchanged.

* Fitted models involving unordered factors may produce
  different results compared with earlier releases because
  of the corrected internal category ordering.

* Users requiring the corrected behaviour should refit
  models created with previous versions.

## Known limitations

* Nominal categories continue to be represented using
  one-dimensional latent intervals. Although the correction
  removes sensitivity to permutations of existing factor
  levels, invariance to arbitrary relabelling of nominal
  categories is not guaranteed.
  
# SynthesizerPlus 0.2.9

* Excel files: date-times are now written as UTC clock time. 'writexl' 2.0
  writes local wall-clock time instead of UTC, while 'readxl' reads Excel
  date-times as UTC, so round trips through `write_data()`/`read_data()` were
  shifted by the time-zone offset (e.g. 2 or 3 hours for Europe/Bucharest).
  The instant is now preserved with every 'writexl' version and with
  'openxlsx'. `?write_data` documents how each format stores date-times.
* `inst/CITATION` and `CITATION.cff` give the Zenodo DOI.
* The multivariate example in `?fit_synthesizer.ts` selected rows with
  `EuStockMarkets[1:500, ]`, which returns a matrix rather than a `ts`, so the
  stationary bootstrap was never used and a warning was printed; it now uses
  `window()`. Passing time-series arguments (`method`, `order`, ...) for a
  non-`ts` input now gives a warning that explains the cause.

# SynthesizerPlus 0.2.8

* DESCRIPTION: the Description now states how discrete and categorical
  variables enter the model (distributional transform, maximum-likelihood
  polychoric and polyserial correlations, with references) and mentions the
  missing-data model, closeness calibration and the correct attribution
  probability.
* Getting-started vignette: the `synthesize(file = )` example now writes a
  file, shows its first lines and reads it back with `read_data()`.

# SynthesizerPlus 0.2.7

Statistical validation, and the bias it uncovered.

## Changes to fitted models

* **Dependence involving discrete or categorical variables was attenuated.**
  The latent correlation of a pair with a discrete margin was estimated by
  correlating its jittered pseudo-observations (distributional transform),
  which biases it towards zero, most strongly for binary variables. In a
  validation design, Kendall's tau between a count and a 4-category factor
  was 0.29 in synthetic data against 0.35 in the training data; under
  missing-at-random (MAR) missingness, the association between `x` and the
  missingness of `y` was reproduced at about 60% of its strength (0.30
  instead of 0.51). Such pairs are now estimated by maximum likelihood under
  the Gaussian copula: polyserial (discrete-continuous) and polychoric
  (discrete-discrete) correlations, with thresholds from the marginal
  probabilities (Olsson 1979; Olsson, Drasgow and Dorans 1982) and the
  bivariate normal CDF of Genz (2004). The estimates agree with `polycor`'s
  two-step estimators to about 1e-5. The t copula uses them for these pairs.
* **Missing data: range restriction.** With incomplete data, continuous
  pairs are now estimated by the pairwise covariance of the normal scores
  scaled by each variable's standard deviation over all its observed values
  (for the t copula: the tau-based correlation rescaled in the same way),
  instead of the pairwise correlation. When missingness of `y` depends on
  `x`, the rows where both are observed have a restricted range; the new
  estimator reproduces the relation among the observed rows (e.g.
  `cor(x, y)` among observed rows 0.568 synthetic vs 0.572 real, where the
  old estimator, combined with the first fix, would give 0.646). With
  complete data nothing changes.
* The latent correlation between a variable and its own missingness
  indicator is now fixed at 0 by design (it is not identified, and this
  keeps the distribution of the observed values); it used to be ~0 by
  accident of the jitter. A discrete margin with a single category has zero
  correlations without being reported as unidentified.
* Fitted models therefore differ from 0.2.6 whenever discrete/categorical
  variables or missing values are present; with only continuous, complete
  variables they are unchanged. Fitting is slower with many discrete
  variables or missingness indicators (e.g. 1.5 s instead of 0.1 s for
  5000 rows with 55 copula columns); the polyserial likelihood is evaluated
  on 500 quantile bins for large samples (difference to the exact estimate
  < 5e-4) and the polychoric likelihood only at the corners of observed
  cells. Synthetic `iris` now keeps the species separation much better (a
  petal-based rule classifies 82% of synthetic records correctly vs 66%;
  real data 96%).
* The fitted copula records the estimator used (`latent`:
  `"normal_scores"` or `"polyserial_polychoric"`).

## Statistical validation suite

* New `tests/testthat/test-validation.R`: data simulated from known
  distributions, dependence structures, missingness mechanisms (MCAR, MAR),
  strata and VAR processes; synthetic samples are compared with the truth
  or the training data within Monte Carlo tolerances (marginals, latent
  correlations, rank association with discrete margins, 2 x 2 tables, MAR
  missingness probabilities, observed-row relations, stratum means and
  correlations, t-copula df and tail dependence, VAR auto- and
  cross-correlations). Run against 0.2.6, it fails 11 expectations in
  exactly the areas affected by the bias above. The large designs are
  skipped on CRAN; one fast MAR check always runs.

## Other changes (seventh code review)

* `fit_synthesizer.ts()` warns when an argument does not apply to the chosen
  method (`block_length` outside the bootstraps, `order`/`burn_in` outside
  `"copula_var"`, `max_lag` when `order` is given).
* `disclosure_risk()` gains `abs_tolerance` for numeric targets: a
  prediction is correct if within `max(tolerance * |y|, abs_tolerance)`, so
  targets with meaningful zeros can be assessed.
* Documentation: `missing = "drop"` leaves columns without any observed value
  entirely missing; missingness that depends on the unobserved value itself
  (MNAR) is not reproduced; the time-series methods assume approximately
  stationary series (`?fit_synthesizer.ts`, time-series vignette).
* The getting-started vignette claimed that the pooled model reproduces the
  June peak of missing ozone values in `airquality`; its own output showed
  otherwise (also in 0.2.6). A single copula only represents monotone
  relations; the section now shows a monotone MAR example and that
  `by = "Month"` reproduces the non-monotone pattern.
* Tests no longer use `testthat::local_mocked_bindings()` (testthat >=
  3.1.7) or `data.frame()` on an `integer64` column (needs bit64), which
  failed on R 4.2 with an older testthat.

# SynthesizerPlus 0.2.6

Bugs fixed.

* New `seed` argument in `fit_synthesizer()` (data frame, matrix, vector and
  `ts` methods). Fitting is stochastic (random tie-breaking, the randomised
  distributional transform, Kendall subsampling for the t copula), so the
  fitted model can now be made reproducible when fitting and generation are
  separate steps; the global RNG state is restored afterwards.
* The warning about unidentified copula correlations now reports the
  maximum number of affected pairs *per fitted model* ("Up to k pair(s) ...
  per fitted model ... affecting m of M fitted model(s)"). Previously, with
  stratified models, `k` could be read as a total.
* Documentation: unidentified correlations are those with fewer than 3
  jointly observed rows *or* a non-finite estimate, so `n_pair` alone does
  not flag all of them (`?copula_correlation`, `?fit_synthesizer`).
  `summary()` labels its count as referring to the pooled model.

# SynthesizerPlus 0.2.5

Bugs fixed.

* t copula: the pairwise counts (`n_pair`) and the identification of
  pairwise correlations were computed on the random 2000-row subsample used
  for Kendall's tau (bug introduced in 0.2.4). With more than 2000 rows,
  `copula_correlation(fit, what = "n_pair")` under-reported the counts, and
  a pair with sparse overlap (e.g. 10 joint rows out of 10,000) could be
  declared unidentified, or estimated from only a few of its rows, depending
  on the seed. Counts and identification now always use the full data, and
  any pair that the subsample leaves with fewer than `min(n_pair, 500)` joint
  rows is re-estimated from up to 2000 of its own joint rows.
* t copula fitting is about 3-4 times faster with many variables: Kendall's
  tau for all pairs is computed as one blocked cross-product of sign
  matrices (exact for the tie-free pseudo-observations) instead of
  `p(p - 1)/2` separate O(n^2) calls of `cor(method = "kendall")`. Estimates
  for complete data are unchanged.
* The fitted copula records how unidentified correlations were filled
  (`completion`: `"none"`, `"maxdet"` or `"zero_fallback"`) and the largest
  change made by the subsequent positive-definiteness repair
  (`pd_adjustment`); `summary()` reports both. The fit-time warning now
  states the method actually used and that the repair may follow, instead of
  always claiming maximum-determinant completion.

# SynthesizerPlus 0.2.4

Bugs fixed.

* Unidentified copula correlations. A pair of variables observed together
  in fewer than 3 rows (e.g. questions from different survey modules) used
  to get correlation 0, read as estimated independence. Worse, a pair seen
  together in exactly 2 rows got a spurious correlation of +/-1, and its
  projection to a valid matrix distorted the *well-estimated* pairs. Such
  pairs are now filled by the maximum-determinant positive-definite
  completion of the identified entries (Dempster, 1972). This keeps every
  identified correlation and makes the pair conditionally independent given
  the other variables; the unknown entries fall back to 0 only if no
  completion exists. `fit_synthesizer()` warns when this happens, the
  fitted copula stores the pairwise counts (`n_pair`) and the number of
  unidentified pairs, `copula_correlation(fit, what = "n_pair")` returns
  the counts, and `summary()` reports them. Fully observed data are
  unaffected.
* `r_mvskewnorm()` requires a positive-definite `omega`, as in the standard
  definition of the skew-normal density; singular `sigma` remains allowed
  for `r_mvnorm()`, `r_mvt()` and `r_mvlnorm()` (degenerate distributions),
  and this is now documented.
* `.check_prob()` (used for `closeness`, `close_quantile`, ...) rejects
  non-finite values explicitly, consistent with the other validators.

# SynthesizerPlus 0.2.3

Bugs fixed.

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

Bugs fixed.

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

Bugs fixed. 

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
