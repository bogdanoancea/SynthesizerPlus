# Synthesizer: fitting and generation ----------------------------------------

#' Fit a synthetic-data model
#'
#' @description
#' `fit_synthesizer()` learns a generative model from a data set. The model
#' combines
#'
#' * a **marginal model per column**: an *interpolated* empirical quantile
#'   function for continuous variables (linear or monotone-spline
#'   interpolation between order statistics), and an empirical probability
#'   table for discrete, categorical, ordinal and logical variables.
#'   Interpolation means that synthetic values of a continuous variable lie
#'   between, rather than on, the observed values, always within the
#'   observed range. With small samples this slightly reduces the spread in
#'   the tails compared with the data;
#' * a **copula** (Gaussian, Student-*t* or independence) that couples the
#'   columns. Discrete and categorical variables take part in the
#'   dependence structure: their pseudo-observations use the distributional
#'   transform, and their latent correlations are estimated by maximum
#'   likelihood (polyserial / polychoric), see *Estimation details*;
#' * an optional **missing-data model**: for each column with missing values
#'   a missingness indicator is added to the copula, so that the rate of
#'   missingness *and* its (monotone) association with other variables are
#'   reproduced;
#' * optional **stratification** (`by`), fitting separate models within
#'   groups, with a pooled fallback for small groups.
#'
#' Use [generate()] to draw synthetic data from the fitted model, or
#' [synthesize()] to fit and generate in one step.
#'
#' @param data A data frame, matrix, atomic vector (numeric, integer,
#'   logical, character, factor, `Date`, `POSIXct`) or `ts` object. Time
#'   series are handled by [fit_synthesizer.ts()].
#' @param copula Dependence model: `"gaussian"` (default), `"t"` (heavier
#'   joint tails) or `"independence"`. See *Estimation details*.
#' @param by Optional character vector of column names used to stratify the
#'   model. Each stratum gets its own marginals and copula; strata are sampled
#'   in proportion to their training size.
#' @param min_stratum_size Strata with fewer rows than this use a pooled model
#'   fitted on all rows.
#' @param missing How missing values are treated: `"joint"` (default; model
#'   missingness inside the copula), `"independent"` (reproduce the missing
#'   rate of each column independently) or `"drop"` (do not reproduce
#'   missingness: partially observed columns are generated without missing
#'   values; a column with *no* observed training value has no distribution
#'   to draw from and necessarily stays entirely missing). Under `"joint"`,
#'   missingness may depend on the other variables (missing at random); it
#'   cannot depend on the unobserved value itself, which the data do not
#'   identify (see *Estimation details*).
#' @param interpolation Interpolation of the empirical quantile function for
#'   continuous variables: `"linear"` (default) or `"spline"` (monotone Hyman
#'   spline, smoother densities).
#' @param discrete_threshold Numeric columns with at most this many distinct
#'   values are treated as discrete (sampled from their observed support),
#'   provided their values are whole numbers or actually repeat (fewer
#'   distinct values than half the observations). A small sample of a
#'   continuous variable is therefore still modelled as continuous.
#' @param knots Optional maximum number of quantile knots stored per
#'   continuous column. `NULL` (default) keeps every observation; a smaller
#'   value shrinks the model and avoids storing the original values.
#' @param id_cols Optional character vector of identifier columns. They are
#'   not modelled; synthetic data get fresh sequential identifiers.
#' @param pd_method How to repair a correlation matrix that is not positive
#'   definite: `"auto"`/`"higham"` (nearest correlation matrix), `"eigen"`
#'   (eigenvalue flooring) or `"none"` (fail if not positive definite).
#' @param seed Optional random seed for the fit. Fitting is stochastic
#'   (random tie-breaking and the randomised distributional transform used
#'   for the pseudo-observations, row subsampling for the *t* copula), so a
#'   seed makes the fitted model reproducible when fitting and generation are
#'   separate steps. The global random number generator state is restored
#'   afterwards. ([synthesize()] has its own `seed`, covering both steps.)
#' @param nonlinear_target Optional name of ONE continuous numeric response to
#'   generate using a hybrid copula--conditional model. `NULL` (default)
#'   leaves the existing copula synthesizer unchanged. The copula generates
#'   all predictors; CART or a penalised GAM predicts the response, local
#'   residuals can be resampled, and the response is quantile-calibrated to its
#'   empirical marginal. Requires `rpart` for CART or `mgcv` for GAM. Not supported
#'   for a `ts` model or a response with discrete/categorical marginals.
#'   The fitted tree and residual pools contain training-data information;
#'   do not release the fitted model as if it were privacy-safe.
#' @param nominal_target Optional name of one unordered categorical response
#'   (factor, character or logical). When supplied, predictors are generated
#'   using the copula and this response is drawn from a fitted multinomial
#'   conditional model (`nnet::multinom`). Nominal categories are not
#'   thresholded along a one-dimensional latent variable in the final draw.
#'   Mutually exclusive with `nonlinear_target`. Experimental, maximum of 20
#'   observed categories and at least 40 observed responses; strata without a
#'   converged model fall back to copula synthesis.
#' @param nominal_decay Nonnegative weight-decay penalty passed to
#'   `nnet::multinom()` (default `0.0001`).
#' @param nominal_maxit Maximum number of optimisation iterations for the
#'   nominal model (default 1000). A model that has not converged within
#'   this limit is not used.
#' @param nominal_spline_df Degrees of freedom of the natural cubic spline
#'   basis used for each numeric predictor (with at least five distinct
#'   values) in the nominal multinomial model, so category probabilities can
#'   depend nonlinearly on numeric predictors. Default 3. The spline model
#'   and the linear model are both fitted and the one with the smaller AIC
#'   is used (recorded in `nominal$spline_df` and `nominal$selection`); `0`
#'   fits the linear model only.
#' @param t_profile If `TRUE` with `copula = "t"`, add exploratory elapsed-time
#'   stage measurements under `fit$pooled$copula$profile` and per-stratum
#'   copulas. Profile fields are absent by default.
#' @param nonlinear_engine Conditional response model for `nonlinear_target`:
#'   `"cart"` (default, unchanged), `"gam"` (penalised additive model
#'   fitted by **mgcv**), or `"auto"` (experimental training-only K-fold
#'   selection among CART, additive GAM and GAM with tensor interactions).
#'   When `"auto"`, the fitted model records the selection scores in
#'   `fit$pooled$nonlinear$selection` (and separately per fitted stratum).
#'   The option is ignored when `nonlinear_target = NULL`.
#' @param nonlinear_cv_folds Integer number of folds for automatic engine
#'   selection, default 5 (2--10). Only non-missing target observations are
#'   used for validation; predictors and target are held out together.
#' @param nonlinear_cv_rule Either `"one_se"` (default: prefer CART, then
#'   additive GAM, then interaction GAM when their average fold score falls
#'   within one standard error of the best) or `"min"` (select the minimum
#'   average fold score). All choices are made using original training
#'   data only; no benchmark holdout observations are consulted.
#' @param nonlinear_cv_score Criterion for automatic engine selection:
#'   `"crps"` (default) scores the whole conditional predictive distribution
#'   (fitted mean plus residual resampling, 50 draws per held-out row) by the
#'   continuous ranked probability score, so a candidate that misrepresents
#'   the conditional spread is penalised; `"mse"` scores only the conditional
#'   mean by its squared prediction error. Both are reported in the
#'   selection table.
#' @param nonlinear_gam_k Maximum basis dimension per GAM smooth (default 10,
#'   permitted 4--50), truncated when too few distinct covariate values occur.
#' @param nonlinear_gam_interactions Whether to include pairwise tensor
#'   interactions among up to three eligible numeric smooth predictors in the
#'   GAM (default `FALSE`). An additive GAM may miss interactions.
#' @param nonlinear_local_k Number of nearest training fitted values for
#'   local GAM residual resampling (default 30). For `nonlinear_engine = "gam"`,
#'   the existing residual mode `"leaf"` means this local fitted-mean
#'   neighbourhood, since a GAM has no CART terminal leaves. This does not
#'   necessarily reproduce variance changes orthogonal to the fitted mean.
#' @param nonlinear_minbucket Minimum observations per CART leaf (default 20).
#' @param nonlinear_cp Complexity penalty for CART splitting (default 0.001).
#' @param nonlinear_maxdepth Maximum CART depth (default 10).
#' @param nonlinear_residual Experimental residual resampling mode for a
#'   hybrid `nonlinear_target`: `"leaf"` (default: CART leaf groups, or
#'   GAM local neighbourhoods defined by fitted conditional means), `"global"` (resample from all fitted residuals),
#'   or `"none"` (use only fitted conditional means). With `"none"`,
#'   identical predictions are calibrated to the same value, so no artificial
#'   within-leaf randomness is added.
#' @param nonlinear_calibration Experimental marginal recalibration of the
#'   hybrid response: `"quantile"` (default) maps each simulated value
#'   through the empirical CDF of a reference sample of conditional draws,
#'   stored at fitting time, and then through the observed response's
#'   empirical quantile function. The result does not depend on how many
#'   rows are generated per call, chunk or stratum. `"none"` uses raw
#'   conditional draws, so its synthetic marginal need not equal the observed
#'   response marginal and values can exceed its observed range.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `sp_synthesizer`: a list whose main
#'   components are `pooled` (the model fitted to all rows: `marginals`,
#'   `copula` and, when requested, `nonlinear` or `nominal`), `strata`
#'   (with `by`: stratum `keys`, `values`, `counts`, `probs`, per-stratum
#'   `models`, `NULL` for strata smaller than `min_stratum_size`, and
#'   `conditional_fallback`, a data frame of the strata whose own model has
#'   no conditional model for the target), and `settings`.
#'
#'   With `nonlinear_engine = "auto"`, `pooled$nonlinear$selection` (and the
#'   same element of each stratum model) holds the selected candidate
#'   (`selected`, `engine`, `interactions`), the `reason` (`"cv"`, or a
#'   fallback reason), `rule`, `score` (`"crps"` or `"mse"`), `folds`, a
#'   `scores` table with per-candidate `mean_score`, `se_score`, `mean_mse`,
#'   `se_mse` and any `issue`, and the per-fold matrices `fold_scores` and
#'   `fold_mse`.
#'
#' @section Conditional models for one target:
#' *Continuous response (`nonlinear_target`).* The copula generates all
#' columns; the response is then replaced by a draw from a conditional model
#' given the generated predictors. Predictors are encoded from all training
#' rows (numeric columns median-imputed plus a missingness indicator;
#' categories, at most 30, plus a `MISSING` level); a category never observed
#' with the response is passed to the tree as unknown and routed by surrogate
#' splits. CART (`rpart`) gives a piecewise-constant conditional mean
#' \eqn{\hat m(x)}{m(x)} and, with residual mode `"leaf"`, residuals are resampled
#' from the centred residuals of the same leaf. The GAM (`mgcv::gam`, REML,
#' cubic regression splines for numeric predictors with at least four
#' distinct values, `select = TRUE` so terms can be penalised to zero,
#' optional pairwise `ti()` interactions among at most three of them)
#' resamples centred residuals from the `nonlinear_local_k` training
#' observations whose fitted means are nearest. Calibration `"quantile"`
#' maps a raw draw \eqn{y^*} to \eqn{F_Y^{-1}(\tilde G(y^*))}{F_Y^{-1}(G(y*))}, where
#' \eqn{F_Y^{-1}} is the empirical quantile function of the observed response
#' and \eqn{\tilde G}{G} the mid-point empirical CDF of 20000 raw conditional
#' draws for resampled training predictor rows, drawn once at fitting time
#' from a fixed random stream. The synthetic response marginal therefore
#' matches the observed one up to Monte Carlo error, independently of batch
#' size.
#'
#' *Automatic engine selection.* CART, an additive GAM and a GAM with tensor
#' interactions are compared by K-fold cross-validation on the training rows
#' with an observed response. With `nonlinear_cv_score = "crps"` each
#' held-out response \eqn{y} is scored by the continuous ranked probability
#' score of \eqn{M = 50} conditional draws \eqn{X_1,\dots,X_M} from the
#' fold model, \eqn{\frac1M\sum_i|X_i-y| - \frac1{2M^2}\sum_{i,j}|X_i-X_j|}{mean|X_i - y| - mean|X_i - X_j| / 2};
#' with `"mse"` by \eqn{(y-\hat m(x))^2}{(y - m(x))^2}. A candidate failing in any fold is
#' discarded. The `"one_se"` rule picks the simplest candidate (CART, then
#' additive GAM, then interaction GAM) whose mean score is within one
#' standard error of the best; the chosen engine is then refitted on all
#' rows.
#'
#' *Unordered categorical response (`nominal_target`).* The response is
#' removed from the copula (its missingness indicator stays) and drawn from
#' a multinomial logit model (`nnet::multinom`, weight decay
#' `nominal_decay`), \eqn{P(Y=k\mid x) \propto \exp(\beta_{k0} + \sum_j
#' f_{kj}(x_j))}{P(Y = k | x) proportional to exp(b_k0 + sum_j f_kj(x_j))}, with \eqn{f_{kj}} a natural cubic spline with
#' `nominal_spline_df` degrees of freedom (knots at quantiles) for each
#' numeric predictor with at least five distinct values, and linear or
#' dummy terms otherwise. The model with linear terms only is also fitted,
#' and of the two converged models the one with the smaller AIC is used:
#' splines are kept when the data support the extra parameters, while for
#' (nearly) separated classes the better-determined linear logit is used.
#' Category probabilities have no interactions between predictors.
#'
#' In a stratified fit, a stratum whose own model cannot fit the requested
#' conditional model (too few observed responses, fewer than two categories,
#' or no convergence) generates the target from its stratum copula; this is
#' reported by a warning and in `strata$conditional_fallback`. The fitted
#' conditional models contain training-data information and are not
#' disclosure-safe artefacts.
#'
#' @section Estimation details:
#' *Latent correlations.* For a pair of continuous variables the Gaussian
#' copula uses the normal scores \eqn{\Phi^{-1}(U)} and the t copula
#' inverts Kendall's tau, \eqn{\rho = \sin(\pi\tau/2)}. As Kendall's tau is
#' quadratic in the sample size, it is computed on a random subsample of 2000
#' rows; a pair that the subsample leaves with fewer than
#' \eqn{\min(n_{ij}, 500)} jointly observed rows is re-estimated from up to
#' 2000 of its own \eqn{n_{ij}} joint rows, so every pair with
#' \eqn{n_{ij} \ge 3} is estimated whatever the subsample.
#'
#' *Discrete and categorical variables.* The pseudo-observations of a
#' discrete margin fill each category's probability interval with
#' independent uniform noise (the distributional transform). Correlating
#' these jittered scores would bias every dependence involving a discrete
#' variable towards zero, most strongly for binary variables. Pairs with at
#' least one discrete margin are therefore estimated by maximum likelihood
#' under the Gaussian copula, with the thresholds fixed at the normal
#' quantiles of the marginal cumulative probabilities: the polyserial
#' correlation for a discrete-continuous pair (Olsson, Drasgow and Dorans,
#' 1982) and the polychoric correlation for a discrete-discrete pair
#' (Olsson, 1979), using the bivariate normal distribution function of Genz
#' (2004). The t copula uses the same estimates for these pairs.
#'
#' *Order of unordered categories.* A categorical variable occupies one
#' latent coordinate, so its categories need an order on it. Ordered factors
#' keep their level order and logical variables `FALSE < TRUE`. For unordered
#' factors and character columns the order is estimated from the data (a
#' form of optimal scaling): each numeric, ordered or logical column is
#' converted to normal scores \eqn{\Phi^{-1}(r/(n+1))}{qnorm(r/(n+1))} from its ranks
#' \eqn{r}; the matrix \eqn{M} of category means of these scores (one row
#' per category) is centred with the category frequencies \eqn{w} as
#' weights, and the categories are sorted by their scores on the first
#' principal component of \eqn{\mathrm{diag}(\sqrt{w})M}{diag(sqrt(w)) M}, the sign fixed so
#' that the first non-zero loading is positive, with ties broken by
#' decreasing frequency and then by label. Categories that are adjacent in
#' their association with the other columns are thus adjacent on the latent
#' coordinate, and renaming or reordering the categories does not change the
#' fitted model. Without informative reference columns, the declared level
#' order (factors) or sorted labels (character) are used. The same order is
#' used when generating, so the latent correlations are the maximum
#' likelihood estimates of the model actually fitted. One coordinate cannot
#' represent every association pattern of a nominal variable; see
#' `nominal_target` and stratification (`by`).
#'
#' *Missing data.* All estimates use **pairwise complete observations**. The
#' rows in which two variables are both observed may be a selected subset
#' (e.g. `y` is observed mainly for low `x`), with a restricted range. For
#' continuous pairs the estimate is therefore the pairwise *covariance* of
#' the normal scores, scaled by each variable's standard deviation over all
#' of its observed values (for the t copula, the tau-based correlation
#' rescaled in the same way). This reproduces the relation among the
#' observed rows under such selection; with complete data it equals the
#' ordinary correlation. With `missing = "joint"`, the missingness
#' indicators are binary variables of the copula, so missingness can depend
#' on the other variables (missing at random). The latent correlation
#' between a variable and its *own* indicator is fixed at 0: the data do not
#' identify it (the value is never seen when missing), and this choice keeps
#' the distribution of the observed values of each variable as in the
#' training data. Pairwise estimates based on different subsets of rows
#' need not form a coherent joint estimate, and the resulting matrix need not
#' be positive definite; it is then projected onto the set of valid
#' correlation matrices according to `pd_method` (Higham's nearest
#' correlation matrix by default). The projection guarantees a valid model,
#' not an optimal joint estimator.
#'
#' *Unidentified pairs.* If two variables are observed together in fewer than
#' 3 rows (e.g. questions from different survey modules), their correlation
#' cannot be estimated (with 2 rows it is trivially \eqn{\pm 1}); the same
#' applies when the pairwise estimate is not finite. Such
#' entries are not set to zero, which would assert independence and, through
#' the projection, distort the identified entries. They are filled by the
#' maximum-determinant positive-definite completion of the identified
#' entries (Dempster, 1972), which keeps every identified correlation and
#' makes the pair conditionally independent given the other variables. If
#' the identified entries cannot be completed, the unknown entries fall back
#' to 0. In either case the matrix may then still be modified by the
#' positive-definiteness repair. A warning reports how many pairs were
#' affected and how they were filled; `copula_correlation(fit, what =
#' "n_pair")` gives the pairwise counts, and the fitted copula records the
#' method (`completion`: `"none"`, `"maxdet"` or `"zero_fallback"`) and the
#' largest change made by the repair (`pd_adjustment`), both shown by
#' `summary()`. Pair counts and identification always refer to the full
#' data.
#'
#' *t-copula degrees of freedom.* With the correlation matrix held fixed,
#' `df` is chosen by profile likelihood over the grid
#' 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100 using the complete rows (at most
#' 5000), so it is a discrete grid estimate rather than a continuous maximum
#' likelihood estimate. If fewer than 10 complete rows are available, `df`
#' is **not estimated** and is fixed at 10; a warning is issued and the
#' fitted copula stores `df_estimated = FALSE` (see [summary()]).
#'
#' @seealso [generate()], [synthesize()], [compare_synthetic()]
#'
#' @references
#' Sklar, A. (1959). Fonctions de répartition à n dimensions et leurs marges.
#' *Publications de l'Institut de Statistique de l'Université de Paris*, 8,
#' 229--231.
#'
#' Dempster, A. P. (1972). Covariance selection. *Biometrics*, 28(1),
#' 157--175. \doi{10.2307/2528966}
#'
#' Genz, A. (2004). Numerical computation of rectangular bivariate and
#' trivariate normal and t probabilities. *Statistics and Computing*, 14(3),
#' 251--260. \doi{10.1023/B:STCO.0000035304.20635.31}
#'
#' Olsson, U. (1979). Maximum likelihood estimation of the polychoric
#' correlation coefficient. *Psychometrika*, 44(4), 443--460.
#' \doi{10.1007/BF02296207}
#'
#' Olsson, U., Drasgow, F. and Dorans, N. J. (1982). The polyserial
#' correlation coefficient. *Psychometrika*, 47(3), 337--347.
#' \doi{10.1007/BF02294164}
#'
#' Rüschendorf, L. (2009). On the distributional transform, Sklar's theorem,
#' and the empirical copula process. *Journal of Statistical Planning and
#' Inference*, 139(11), 3921--3927. \doi{10.1016/j.jspi.2009.05.030}
#'
#' Higham, N. J. (2002). Computing the nearest correlation matrix -- a
#' problem from finance. *IMA Journal of Numerical Analysis*, 22(3),
#' 329--343. \doi{10.1093/imanum/22.3.329}
#'
#' @examples
#' fit <- fit_synthesizer(iris)
#' fit
#' syn <- generate(fit, n = 300, seed = 1)
#' head(syn)
#'
#' # Stratified model and missing-data model
#' fit_aq <- fit_synthesizer(airquality, by = "Month")
#' colSums(is.na(generate(fit_aq, seed = 2)))
#' @export
fit_synthesizer <- function(data, ...) {
  UseMethod("fit_synthesizer")
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.data.frame <- function(data,
                                       copula = c("gaussian", "t", "independence"),
                                       by = NULL,
                                       min_stratum_size = 20L,
                                       missing = c("joint", "independent", "drop"),
                                       interpolation = c("linear", "spline"),
                                       discrete_threshold = 20L,
                                       knots = NULL,
                                       id_cols = NULL,
                                       pd_method = c("auto", "higham", "eigen", "none"),
                                       nonlinear_target = NULL,
                                       nominal_target = NULL,
                                       nominal_decay = 0.0001,
                                       nominal_maxit = 1000L,
                                       nominal_spline_df = 3L,
                                       t_profile = FALSE,
                                       nonlinear_engine = c("cart", "gam", "auto"),
                                       nonlinear_cv_folds = 5L,
                                       nonlinear_cv_rule = c("one_se", "min"),
                                       nonlinear_cv_score = c("crps", "mse"),
                                       nonlinear_gam_k = 10L,
                                       nonlinear_gam_interactions = FALSE,
                                       nonlinear_local_k = 30L,
                                       nonlinear_minbucket = 20L,
                                       nonlinear_cp = 0.001,
                                       nonlinear_maxdepth = 10L,
                                       nonlinear_residual = c("leaf", "global", "none"),
                                       nonlinear_calibration = c("quantile", "none"),
                                       seed = NULL, ...) {
  if (!is.null(seed)) {
    # fitting is stochastic (random tie-breaking and the randomised
    # distributional transform of the pseudo-observations, Kendall
    # subsampling); a seed makes the fitted model reproducible
    return(.with_seed(seed, fit_synthesizer.data.frame(
      data, copula = copula, by = by, min_stratum_size = min_stratum_size,
      missing = missing, interpolation = interpolation,
      discrete_threshold = discrete_threshold, knots = knots, id_cols = id_cols,
      pd_method = pd_method, nonlinear_target = nonlinear_target,
      nominal_target = nominal_target, nominal_decay = nominal_decay,
      nominal_maxit = nominal_maxit, nominal_spline_df = nominal_spline_df,
      t_profile = t_profile,
      nonlinear_engine = nonlinear_engine,
      nonlinear_cv_folds = nonlinear_cv_folds,
      nonlinear_cv_rule = nonlinear_cv_rule,
      nonlinear_cv_score = nonlinear_cv_score,
      nonlinear_gam_k = nonlinear_gam_k,
      nonlinear_gam_interactions = nonlinear_gam_interactions,
      nonlinear_local_k = nonlinear_local_k,
      nonlinear_minbucket = nonlinear_minbucket, nonlinear_cp = nonlinear_cp,
      nonlinear_maxdepth = nonlinear_maxdepth,
      nonlinear_residual = nonlinear_residual,
      nonlinear_calibration = nonlinear_calibration, seed = NULL, ...)))
  }
  copula <- match.arg(copula)
  missing <- match.arg(missing)
  interpolation <- match.arg(interpolation)
  pd_method <- match.arg(pd_method)
  dot_names <- names(list(...))
  .warn_unused_dots(list(...), "fit_synthesizer",
                    if (any(c("closeness", "perturb", "dependence") %in% dot_names)) {
                      "'closeness', 'perturb' and 'dependence' are arguments of generate()."
                    } else if (any(c("method", "max_lag", "order", "block_length", "burn_in") %in% dot_names)) {
                      paste("Time-series arguments apply only to 'ts' objects; note that row",
                            "subsetting such as x[1:100, ] drops the ts class (use window()).")
                    })
  data <- .as_plain_df(data)
  if (nrow(data) < 2L) stop("'data' must have at least 2 rows.", call. = FALSE)
  if (ncol(data) < 1L) stop("'data' must have at least 1 column.", call. = FALSE)
  if (anyDuplicated(names(data)) || any(is.na(names(data)) | names(data) == "")) {
    stop("'data' must have unique, non-empty column names.", call. = FALSE)
  }
  discrete_threshold <- .check_count(discrete_threshold, "discrete_threshold")
  min_stratum_size <- .check_count(min_stratum_size, "min_stratum_size")
  if (!is.null(knots)) knots <- .check_count(knots, "knots", allow_zero = FALSE)
  if (!is.null(knots) && knots < 2L) stop("'knots' must be at least 2.", call. = FALSE)

  if (!is.logical(t_profile) || length(t_profile) != 1L || is.na(t_profile)) {
    stop("'t_profile' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(nominal_target)) {
    if (!is.null(nonlinear_target)) {
      stop("'nominal_target' and 'nonlinear_target' are currently mutually exclusive.",
           call. = FALSE)
    }
    if (!is.character(nominal_target) || length(nominal_target) != 1L ||
        is.na(nominal_target) || !nzchar(nominal_target)) {
      stop("'nominal_target' must identify one categorical column.", call. = FALSE)
    }
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("'nominal_target' requires the suggested package 'nnet'.", call. = FALSE)
    }
    nominal_maxit <- .check_count(nominal_maxit, "nominal_maxit", allow_zero = FALSE)
    nominal_spline_df <- .check_count(nominal_spline_df, "nominal_spline_df", allow_zero = TRUE)
    if (nominal_spline_df == 1L || nominal_spline_df > 10L) {
      stop("'nominal_spline_df' must be 0 (linear terms) or between 2 and 10.", call. = FALSE)
    }
    if (!is.numeric(nominal_decay) || length(nominal_decay) != 1L ||
        !is.finite(nominal_decay) || nominal_decay < 0) {
      stop("'nominal_decay' must be a nonnegative finite number.", call. = FALSE)
    }
  }
  if (!is.null(nonlinear_target)) {
    nonlinear_engine <- match.arg(nonlinear_engine)
    nonlinear_cv_rule <- match.arg(nonlinear_cv_rule)
    nonlinear_cv_score <- match.arg(nonlinear_cv_score)
    nonlinear_cv_folds <- .check_count(nonlinear_cv_folds, "nonlinear_cv_folds",
                                       allow_zero = FALSE)
    if (nonlinear_cv_folds < 2L || nonlinear_cv_folds > 10L) {
      stop("'nonlinear_cv_folds' must be between 2 and 10.", call. = FALSE)
    }
    nonlinear_residual <- match.arg(nonlinear_residual)
    nonlinear_calibration <- match.arg(nonlinear_calibration)
    if (!is.character(nonlinear_target) || length(nonlinear_target) != 1L ||
        is.na(nonlinear_target) || !nzchar(nonlinear_target)) {
      stop("'nonlinear_target' must be NULL or the name of one numeric response.", call. = FALSE)
    }
    conditional_packages <- if (nonlinear_engine == "gam") "mgcv" else "rpart"
    if (!all(vapply(conditional_packages, requireNamespace, logical(1), quietly = TRUE))) {
      stop("The hybrid nonlinear engine '", nonlinear_engine,
           "' requires the suggested package '", conditional_packages[1L], "'.",
           call. = FALSE)
    }
    if (nonlinear_engine == "auto" && !requireNamespace("mgcv", quietly = TRUE)) {
      warning("Automatic engine selection: 'mgcv' is unavailable; evaluating CART only.",
              call. = FALSE)
    }
    nonlinear_gam_k <- .check_count(nonlinear_gam_k, "nonlinear_gam_k",
                                    allow_zero = FALSE)
    if (nonlinear_gam_k < 4L || nonlinear_gam_k > 50L) {
      stop("'nonlinear_gam_k' must be between 4 and 50.", call. = FALSE)
    }
    nonlinear_local_k <- .check_count(nonlinear_local_k, "nonlinear_local_k",
                                      allow_zero = FALSE)
    if (!is.logical(nonlinear_gam_interactions) ||
        length(nonlinear_gam_interactions) != 1L ||
        is.na(nonlinear_gam_interactions)) {
      stop("'nonlinear_gam_interactions' must be TRUE or FALSE.", call. = FALSE)
    }
    nonlinear_minbucket <- .check_count(nonlinear_minbucket, "nonlinear_minbucket",
                                         allow_zero = FALSE)
    nonlinear_maxdepth <- .check_count(nonlinear_maxdepth, "nonlinear_maxdepth",
                                       allow_zero = FALSE)
    if (nonlinear_maxdepth > 30L) {
      stop("'nonlinear_maxdepth' must be at most 30.", call. = FALSE)
    }
    if (!is.numeric(nonlinear_cp) || length(nonlinear_cp) != 1L ||
        is.na(nonlinear_cp) || !is.finite(nonlinear_cp) ||
        nonlinear_cp < 0 || nonlinear_cp > 1) {
      stop("'nonlinear_cp' must be a finite number in [0, 1].", call. = FALSE)
    }
  }

  .check_cols(by, names(data), "by")
  .check_cols(id_cols, names(data), "id_cols")
  if (length(intersect(by, id_cols))) {
    stop("A column cannot be both in 'by' and 'id_cols'.", call. = FALSE)
  }
  model_cols <- setdiff(names(data), c(by, id_cols))
  if (!length(model_cols)) stop("No columns left to model.", call. = FALSE)
  if (!is.null(nominal_target)) {
    if (!nominal_target %in% model_cols) {
      stop("'nominal_target' must be a modelled column.", call. = FALSE)
    }
    if (!.value_class(data[[nominal_target]]) %in% c("factor", "character", "logical")) {
      stop("'nominal_target' must be an unordered factor, character or logical column.",
           call. = FALSE)
    }
    if (length(model_cols) < 2L) {
      stop("Nominal conditional synthesis requires predictors.", call. = FALSE)
    }
  }
  if (!is.null(nonlinear_target)) {
    if (!nonlinear_target %in% model_cols) {
      stop("'nonlinear_target' must name a modelled column (not 'by' or 'id_cols').",
           call. = FALSE)
    }
    y <- data[[nonlinear_target]]
    if (!is.double(y) || !identical(.value_class(y), "numeric") ||
        any(!is.na(y) & !is.finite(y))) {
      stop("'nonlinear_target' must be a finite, double-precision numeric column.",
           call. = FALSE)
    }
    if (length(model_cols) < 2L) {
      stop("Hybrid synthesis requires at least one predictor column.", call. = FALSE)
    }
  }
  bad <- model_cols[is.na(vapply(data[model_cols], .value_class, character(1)))]
  if (length(bad)) {
    stop("Unsupported column class in: ", paste(bad, collapse = ", "),
         ". Supported: numeric, integer, logical, character, factor, ordered, Date, POSIXct ",
         "and numeric classes such as difftime.",
         call. = FALSE)
  }

  settings <- list(copula = copula, missing = missing,
                   interpolation = interpolation,
                   discrete_threshold = discrete_threshold, knots = knots,
                   pd_method = pd_method)
  if (isTRUE(t_profile) && copula == "t") settings$t_profile <- TRUE
  if (!is.null(nominal_target)) {
    settings$nominal_target <- nominal_target
    settings$nominal_decay <- nominal_decay
    settings$nominal_maxit <- nominal_maxit
    settings$nominal_spline_df <- nominal_spline_df
  }
  # Preserve the original fitted-object schema when the option is not used.
  if (!is.null(nonlinear_target)) {
    settings$nonlinear_target <- nonlinear_target
    settings$nonlinear_engine <- nonlinear_engine
    if (identical(nonlinear_engine, "auto")) {
      settings$nonlinear_cv_folds <- nonlinear_cv_folds
      settings$nonlinear_cv_rule <- nonlinear_cv_rule
      settings$nonlinear_cv_score <- nonlinear_cv_score
    }
    settings$nonlinear_gam_k <- nonlinear_gam_k
    settings$nonlinear_gam_interactions <- nonlinear_gam_interactions
    settings$nonlinear_local_k <- nonlinear_local_k
    settings$nonlinear_minbucket <- nonlinear_minbucket
    settings$nonlinear_cp <- nonlinear_cp
    settings$nonlinear_maxdepth <- nonlinear_maxdepth
    settings$nonlinear_residual <- nonlinear_residual
    settings$nonlinear_calibration <- nonlinear_calibration
  }

  pooled <- .fit_stratum(data[model_cols], settings)
  if (!is.null(nonlinear_target) &&
      pooled$marginals[[nonlinear_target]]$type != "continuous") {
    stop("'nonlinear_target' must have a continuous empirical marginal.",
         call. = FALSE)
  }
  if (!is.null(nominal_target) && is.null(pooled$nominal)) {
    stop("Could not fit a converged nominal conditional model: require at least 40 ",
         "observed responses, 2-20 categories, and a multinomial fit that converges ",
         "within 'nominal_maxit' iterations (consider increasing it, or a smaller ",
         "'nominal_spline_df').",
         call. = FALSE)
  }
  if (!is.null(nonlinear_target) && is.null(pooled$nonlinear)) {
    stop("Not enough observed target values to fit the nonlinear conditional model: ",
         "at least max(10, 2 * nonlinear_minbucket) are required.", call. = FALSE)
  }

  strata <- NULL
  if (length(by)) {
    key <- .stratum_key(data[by])
    keys <- unique(key)                      # strata in order of first appearance
    counts <- tabulate(match(key, keys), nbins = length(keys))
    first <- match(keys, key)
    values <- data[first, by, drop = FALSE]
    attr(values, "row.names") <- .set_row_names(length(keys))
    models <- lapply(keys, function(k) {
      rows <- which(key == k)
      if (length(rows) < min_stratum_size) return(NULL)
      .fit_stratum(data[rows, model_cols, drop = FALSE], settings)
    })
    names(models) <- keys
    strata <- list(keys = keys, values = values,
                   probs = counts / sum(counts),
                   counts = as.integer(counts), models = models)
    # Strata with their own copula but without the requested conditional
    # model fall back to the stratum copula for the target: record and say so.
    cond <- if (!is.null(nonlinear_target)) "nonlinear" else if (!is.null(nominal_target)) "nominal"
    if (!is.null(cond)) {
      lacking <- vapply(models, function(m) !is.null(m) && is.null(m[[cond]]), logical(1))
      fallback <- values[lacking, , drop = FALSE]
      attr(fallback, "row.names") <- .set_row_names(nrow(fallback))
      strata$conditional_fallback <- fallback
      if (any(lacking)) {
        warning(sprintf(paste0("The %s conditional model could not be fitted in %d of %d ",
                               "stratum model(s) (too few observed target values or no ",
                               "converged fit); the target is generated by the stratum copula ",
                               "there. They are listed in 'fit$strata$conditional_fallback'."),
                        cond, sum(lacking), sum(!vapply(models, is.null, logical(1)))),
                call. = FALSE)
      }
    }
  }

  mods <- c(list(pooled), if (!is.null(strata)) strata$models)
  mods <- mods[!vapply(mods, is.null, logical(1))]
  # correlations of pairs never (or almost never) observed together
  msg <- .unidentified_message(mods)
  if (!is.null(msg)) warning(msg, call. = FALSE)
  # t copula: say so when df could not be estimated (too few complete rows)
  if (copula == "t") {
    fixed <- vapply(mods, function(m) identical(m$copula$df_estimated, FALSE), logical(1))
    if (any(fixed)) {
      warning(sprintf(paste0("t copula: fewer than 10 complete rows in %d of %d fitted model(s); ",
                             "df was not estimated there and is fixed at 10."),
                      sum(fixed), length(fixed)), call. = FALSE)
    }
  }

  structure(list(
    prototype = .prototype(data),
    columns = names(data),
    model_cols = model_cols,
    by = by,
    id_cols = id_cols,
    pooled = pooled,
    strata = strata,
    settings = settings,
    n_train = nrow(data),
    input = "data.frame"
  ), class = "sp_synthesizer")
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.matrix <- function(data, ...) {
  if (is.null(colnames(data))) colnames(data) <- paste0("V", seq_len(ncol(data)))
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  obj <- fit_synthesizer.data.frame(df, ...)
  obj$input <- "matrix"
  obj
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.default <- function(data, ...) {
  if (!is.atomic(data) || !is.null(dim(data)) || is.na(.value_class(data))) {
    stop("Cannot fit a synthesizer to an object of class '",
         paste(class(data), collapse = "/"), "'.", call. = FALSE)
  }
  obj <- fit_synthesizer.data.frame(data.frame(x = data, stringsAsFactors = FALSE), ...)
  obj$input <- "vector"
  obj
}

.check_cols <- function(cols, available, arg) {
  if (is.null(cols)) return(invisible(NULL))
  if (!is.character(cols) || anyNA(cols)) {
    stop(sprintf("'%s' must be a character vector of column names.", arg), call. = FALSE)
  }
  miss <- setdiff(cols, available)
  if (length(miss)) {
    stop(sprintf("'%s' columns not found in data: %s", arg,
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

# Collision-free row keys. Each column is coded by matching its values
# against the union of values in all data frames given (NA gets its own code,
# distinct from any string such as "<NA>"); keys combine the integer codes,
# so no user value can create a collision. Returns one key vector per data
# frame, with codes consistent across them.
.row_keys <- function(...) {
  dfs <- list(...)
  n <- vapply(dfs, nrow, integer(1))
  cols <- names(dfs[[1L]])
  codes <- lapply(cols, function(nm) {
    vals <- lapply(dfs, function(d) {
      v <- d[[nm]]
      if (is.factor(v)) v <- as.character(v)
      if (is.double(v)) v[is.nan(v)] <- NA     # NaN and NA are both missing
      v
    })
    all <- do.call(c, unname(lapply(vals, function(v) if (inherits(v, c("Date", "POSIXct"))) as.numeric(v) else v)))
    match(all, unique(all))
  })
  key <- if (length(codes)) do.call(paste, c(codes, sep = ".")) else rep("", sum(n))
  split(key, rep(seq_along(dfs), n))
}

.stratum_key <- function(df) .row_keys(df)[[1L]]

# Human-readable labels for strata, built from the stratum values
.stratum_labels <- function(values) {
  parts <- lapply(values, function(v) ifelse(is.na(v), "NA", as.character(v)))
  do.call(paste, c(parts, sep = " / "))
}

# Fit marginals + copula on one block of data
.fit_stratum <- function(df, settings) {
  marg <- lapply(names(df), function(nm) {
    .fit_marginal(df[[nm]], nm,
                  interpolation = settings$interpolation,
                  discrete_threshold = settings$discrete_threshold,
                  knots = settings$knots)
  })
  names(marg) <- names(df)
  # categories of unordered categorical columns: data-driven latent order
  marg <- .order_categories(marg, df)

  # Fit before assembling copula coordinates: if multinomial optimisation is
  # unavailable in a sparse stratum, retain the COMPLETE baseline copula,
  # including the target. This is a true copula fallback, not independence.
  nominal <- if (!is.null(settings$nominal_target)) .fit_nominal(df, settings) else NULL

  # A conditionally generated nominal target is deliberately omitted from the
  # latent copula: otherwise arbitrary category labels can still contaminate
  # correlations among predictor columns through positive-definiteness repair.
  # Its missingness indicator (if requested) remains in the joint copula.
  copula_marg <- marg[vapply(marg, function(m) m$type != "empty", logical(1))]
  if (!is.null(nominal)) {
    copula_marg[[settings$nominal_target]] <- NULL
  }
  U <- lapply(names(copula_marg), function(nm) .pseudo_obs(copula_marg[[nm]], df[[nm]]))
  names(U) <- names(copula_marg)

  miss_marg <- list()
  if (settings$missing == "joint") {
    has_na <- vapply(marg, function(m) m$na_rate > 0 && m$na_rate < 1, logical(1))
    for (nm in names(marg)[has_na]) {
      mm <- .missing_marginal(marg[[nm]])
      miss_marg[[nm]] <- mm
      ind <- ifelse(is.na(df[[nm]]), "missing", "observed")
      U[[mm$name]] <- .pseudo_obs(mm, ind)
    }
  }

  Umat <- if (length(U)) do.call(cbind, U) else matrix(numeric(0), nrow(df), 0L)
  if (length(U)) colnames(Umat) <- names(U)
  # cut points of discrete margins (for maximum-likelihood latent
  # correlations) and the pairs fixed at 0 (variable, own missingness)
  cut_of <- function(m) {
    if (m$type %in% c("discrete", "categorical")) {
      cc <- c(0, cumsum(m$probs))
      cc[length(cc)] <- 1
      cc
    } else {
      NULL
    }
  }
  all_marg <- c(copula_marg, stats::setNames(miss_marg, vapply(miss_marg, `[[`, "", "name")))
  cuts <- lapply(colnames(Umat), function(nm) cut_of(all_marg[[nm]]))
  fixed_zero <- matrix(FALSE, ncol(Umat), ncol(Umat),
                       dimnames = list(colnames(Umat), colnames(Umat)))
  for (nm in names(miss_marg)) {
    if (nm %in% colnames(Umat)) {
      fixed_zero[nm, miss_marg[[nm]]$name] <- TRUE
      fixed_zero[miss_marg[[nm]]$name, nm] <- TRUE
    }
  }
  cop <- .fit_copula(Umat, family = settings$copula,
                     pd_method = settings$pd_method, cuts = cuts,
                     fixed_zero = fixed_zero,
                     profile = isTRUE(settings$t_profile))
  hybrid <- if (!is.null(settings$nonlinear_target) &&
                marg[[settings$nonlinear_target]]$type == "continuous") {
    .fit_nonlinear(df, settings)
  } else {
    NULL
  }
  ans <- list(marginals = marg, missing_marginals = miss_marg,
              copula = cop, n = nrow(df))
  if (!is.null(settings$nominal_target)) ans$nominal <- nominal
  if (!is.null(settings$nonlinear_target)) ans$nonlinear <- hybrid
  ans
}

# Warning text about unidentified copula correlations across the fitted
# models (pooled and strata), or NULL when there are none
.unidentified_message <- function(mods) {
  n_unid <- vapply(mods, function(m) m$copula$unidentified %||% 0L, integer(1))
  if (!any(n_unid > 0L)) return(NULL)
  meth <- vapply(mods[n_unid > 0L], function(m) m$copula$completion %||% "maxdet",
                 character(1))
  how <- if (all(meth == "maxdet")) {
    "were filled by maximum-determinant completion"
  } else if (all(meth == "zero_fallback")) {
    "could not be completed consistently and were set to 0"
  } else {
    sprintf("were filled by maximum-determinant completion (set to 0 in %d model(s) where no completion exists)",
            sum(meth == "zero_fallback"))
  }
  sprintf(paste0("Up to %d pair(s) of variables per fitted model are observed together in ",
                 "fewer than 3 rows or have a non-finite estimate (affecting %d of %d ",
                 "fitted model(s)); their copula correlations are not identified and %s, ",
                 "before any positive-definiteness repair. ",
                 "See copula_correlation(fit, what = \"n_pair\") and summary(fit)."),
          max(n_unid), sum(n_unid > 0L), length(mods), how)
}

# Generate n rows from one fitted block
.generate_stratum <- function(model, n, pert, prototype, missing) {
  marg <- model$marginals
  cop <- model$copula
  vars <- colnames(cop$R)
  U <- .rcopula_fitted(n, cop, dependence = pert$dependence)
  colnames(U) <- vars
  out <- vector("list", length(marg))
  names(out) <- names(marg)
  for (nm in names(marg)) {
    m <- marg[[nm]]
    u <- if (nm %in% vars) U[, nm] else stats::runif(n)
    v <- .restore_class(.perturbed_values(m, u, pert), m, prototype[[nm]])
    if (m$type != "empty" && m$na_rate > 0) {
      if (missing == "joint" && !is.null(model$missing_marginals[[nm]])) {
        mm <- model$missing_marginals[[nm]]
        is_na <- .marginal_quantile(mm, U[, mm$name]) == "missing"
      } else if (missing == "independent") {
        is_na <- stats::runif(n) < m$na_rate
      } else {
        is_na <- rep(FALSE, n)
      }
      v[is_na] <- NA
    }
    out[[nm]] <- v
  }
  attr(out, "row.names") <- .set_row_names(n)
  class(out) <- "data.frame"
  if (!is.null(model$nonlinear) && n > 0L) {
    target <- model$nonlinear$target
    target_missing <- is.na(out[[target]])
    transformed <- .generate_nonlinear(model$nonlinear, out, marg[[target]])
    out[[target]] <- .restore_class(transformed, marg[[target]], prototype[[target]])
    if (marg[[target]]$na_rate > 0) {
      # Keep the missingness indicator already drawn from the fitted copula.
      # Save the target mask before replacing its simulated values.
      out[[target]][target_missing] <- NA
    }
  }
  if (!is.null(model$nominal) && n > 0L) {
    target <- model$nominal$target
    target_missing <- is.na(out[[target]])
    sampled <- .generate_nominal(model$nominal, out)
    out[[target]] <- .restore_class(sampled, marg[[target]], prototype[[target]])
    if (marg[[target]]$na_rate > 0) out[[target]][target_missing] <- NA
  }
  out
}

#' Generate synthetic data from a fitted model
#'
#' @description
#' `generate()` draws synthetic records from a fitted model such as an
#' `sp_synthesizer` (see [fit_synthesizer()]) or an `sp_ts_synthesizer` (see
#' [fit_synthesizer.ts()]).
#'
#' @param object A fitted model.
#' @param n Number of records (rows) to generate. Defaults to the number of
#'   rows in the training data.
#' @param seed Optional random seed. The global random number generator state
#'   is restored afterwards.
#' @param dependence Number in \[0, 1\] scaling the strength of the dependence
#'   between variables: 1 (default) reproduces the fitted copula, 0 generates
#'   independent columns. Intermediate values shrink the copula correlation
#'   matrix towards the identity; the marginal distributions are unaffected.
#' @param closeness Number in \[0, 1\] controlling how close the synthetic
#'   data are to the real data: 1 (default) reproduces the fitted model as
#'   closely as possible; smaller values move records away from the real
#'   ones (noise on continuous variables, random re-drawing of categories) and
#'   loosen the distributions (flatter category probabilities, weaker
#'   dependence). With `closeness = 0` the variables are independent and
#'   strongly smoothed. See [perturbation()] for the exact mapping and
#'   [calibrate_closeness()] to choose a value that reaches a target distance.
#' @param perturb Optional [perturbation()] for fine control: record-level
#'   noise and swapping, and distribution-level shifts, spread, category
#'   temperature and dependence strength. Combined with `closeness` and
#'   `dependence`.
#' @param stratum_sizes For stratified models: `"multinomial"` (default;
#'   random stratum sizes) or `"proportional"` (deterministic sizes
#'   proportional to the training data).
#' @param ... Further arguments passed to methods.
#'
#' @details
#' For a fit with `nonlinear_target` or `nominal_target`, that one response is
#' generated from its conditional model given the copula-generated
#' predictors, rather than directly from the fitted copula. The fitted copula
#' correlation is therefore a *baseline* model, not necessarily the final
#' synthetic joint dependence. The missingness of the response still follows
#' the copula. Non-default `dependence`, `closeness` and `perturb` are not
#' supported for such fits.
#'
#' @return Synthetic data of the same type as the training data: a data
#'   frame, matrix, vector or `ts` object.
#'
#' @seealso [fit_synthesizer()], [synthesize()], [generate_to_file()]
#' @examples
#' fit <- fit_synthesizer(mtcars)
#' generate(fit, n = 5, seed = 42)
#'
#' # weaker dependence
#' syn <- generate(fit, n = 500, dependence = 0.5, seed = 1)
#' @export
generate <- function(object, ...) {
  UseMethod("generate")
}

#' @rdname generate
#' @export
generate.sp_synthesizer <- function(object, n = object$n_train, seed = NULL,
                                    dependence = 1, closeness = 1, perturb = NULL,
                                    stratum_sizes = c("multinomial", "proportional"),
                                    ...) {
  .warn_unused_dots(list(...), "generate")
  n <- .check_count(n)
  dependence <- .check_prob(dependence, "dependence")
  closeness <- .check_prob(closeness, "closeness")
  stratum_sizes <- match.arg(stratum_sizes)
  if ((!is.null(object$settings$nonlinear_target) ||
       !is.null(object$settings$nominal_target)) &&
      (dependence != 1 || closeness != 1 || !is.null(perturb))) {
    stop("Conditional synthesis currently requires dependence = 1, ",
         "closeness = 1 and perturb = NULL; these controls do not have ",
         "the same interpretation after conditional synthesis.", call. = FALSE)
  }
  pert <- .combine_perturbation(.closeness_perturbation(closeness), perturb)
  pert$dependence <- pert$dependence * dependence
  .check_pert_names(pert, object$model_cols, object$by, object$id_cols)
  .with_seed(seed, .generate_impl(object, n, pert, stratum_sizes))
}

.generate_impl <- function(object, n, pert, stratum_sizes) {
  proto <- object$prototype
  missing <- object$settings$missing

  if (is.null(object$strata)) {
    body <- .generate_stratum(object$pooled, n, pert, proto, missing)
  } else {
    st <- object$strata
    sizes <- if (stratum_sizes == "multinomial") {
      as.integer(stats::rmultinom(1L, n, st$probs))
    } else {
      .largest_remainder(n, st$probs)
    }
    parts <- vector("list", length(st$keys))
    for (i in seq_along(st$keys)) {
      if (sizes[i] == 0L) next
      mod <- st$models[[i]] %||% object$pooled
      blk <- .generate_stratum(mod, sizes[i], pert, proto, missing)
      for (b in object$by) blk[[b]] <- rep(st$values[[b]][i], sizes[i])
      parts[[i]] <- blk
    }
    parts <- parts[!vapply(parts, is.null, logical(1))]
    body <- if (length(parts)) do.call(rbind, parts) else NULL
    if (!is.null(body) && nrow(body) > 1L) {
      body <- body[sample.int(nrow(body)), , drop = FALSE]
    }
  }

  out <- proto[rep(NA_integer_, n), , drop = FALSE]
  if (!is.null(body)) {
    for (nm in setdiff(object$columns, object$id_cols)) out[[nm]] <- body[[nm]]
  }
  for (nm in object$id_cols) out[[nm]] <- .make_ids(n, proto[[nm]])
  attr(out, "row.names") <- .set_row_names(n)

  switch(object$input,
    vector = out[[1L]],
    matrix = as.matrix(out),
    out
  )
}

.largest_remainder <- function(n, p) {
  raw <- n * p
  base <- floor(raw)
  rem <- n - sum(base)
  if (rem > 0) {
    idx <- order(raw - base, decreasing = TRUE)[seq_len(rem)]
    base[idx] <- base[idx] + 1
  }
  as.integer(base)
}

.make_ids <- function(n, template) {
  ids <- seq_len(n)
  if (is.factor(template)) return(factor(sprintf("syn%0*d", nchar(n), ids)))
  if (is.character(template)) return(sprintf("syn%0*d", nchar(n), ids))
  if (is.integer(template)) return(ids)
  if (is.numeric(template)) return(as.numeric(ids))
  ids
}

#' Fit a synthesizer and generate synthetic data in one step
#'
#' @description
#' `synthesize()` is a convenience wrapper around [fit_synthesizer()] and
#' [generate()]. Optionally the result is written to a file in any format
#' supported by [write_data()].
#'
#' @param data Training data (data frame, matrix, vector or `ts`).
#' @param n Number of synthetic records. Defaults to the size of `data`.
#' @param seed Optional random seed (applies to fitting and generation).
#' @param dependence Dependence strength in \[0, 1\]; see [generate()].
#' @param closeness,perturb Closeness to the real data and optional
#'   perturbation; see [generate()].
#' @param file Optional output file path. The format is inferred from the
#'   extension (see [supported_formats()]).
#' @param format Optional explicit file format, overriding the extension.
#' @param ... Arguments passed to [fit_synthesizer()].
#'
#' @return The synthetic data (invisibly if `file` is given).
#'
#' @examples
#' syn <- synthesize(iris, n = 200, seed = 1)
#' summary(syn)
#'
#' out <- tempfile(fileext = ".csv")
#' synthesize(mtcars, file = out, seed = 1)
#' head(read_data(out))
#' @export
synthesize <- function(data, n = NULL, seed = NULL, dependence = 1,
                       closeness = 1, perturb = NULL,
                       file = NULL, format = NULL, ...) {
  res <- .with_seed(seed, {
    fit <- fit_synthesizer(data, ...)
    if (inherits(fit, "sp_ts_synthesizer")) {
      if (dependence != 1 || closeness != 1 || !is.null(perturb)) {
        warning(.ts_closeness_msg, call. = FALSE)
      }
      generate(fit, n = n %||% fit$n_train)
    } else {
      generate(fit, n = n %||% fit$n_train, dependence = dependence,
               closeness = closeness, perturb = perturb)
    }
  })
  if (!is.null(file)) {
    write_data(res, file, format = format)
    return(invisible(res))
  }
  res
}

#' Simulate several synthetic data sets
#'
#' Draws `nsim` independent synthetic data sets from a fitted model, for
#' example to account for synthesis uncertainty with multiple synthetic data
#' sets (Raghunathan, Reiter and Rubin, 2003).
#'
#' @param object A fitted `sp_synthesizer` or `sp_ts_synthesizer`.
#' @param nsim Number of synthetic data sets.
#' @param seed Optional random seed.
#' @param ... Passed to [generate()], e.g. `n`.
#' @return A list of `nsim` synthetic data sets.
#' @references Raghunathan, T. E., Reiter, J. P. and Rubin, D. B. (2003).
#'   Multiple imputation for statistical disclosure limitation. *Journal of
#'   Official Statistics*, 19(1), 1--16.
#' @examples
#' fit <- fit_synthesizer(mtcars)
#' sims <- simulate(fit, nsim = 3, seed = 1)
#' vapply(sims, function(d) mean(d$mpg), numeric(1))
#' @importFrom stats simulate
#' @export
simulate.sp_synthesizer <- function(object, nsim = 1, seed = NULL, ...) {
  nsim <- .check_count(nsim, "nsim", allow_zero = FALSE)
  .with_seed(seed, lapply(seq_len(nsim), function(i) generate(object, ...)))
}

#' @rdname simulate.sp_synthesizer
#' @export
simulate.sp_ts_synthesizer <- simulate.sp_synthesizer

#' Inspect a fitted synthesizer
#'
#' `copula_correlation()` returns the (latent) copula correlation matrix of a
#' fitted model, or with `what = "n_pair"` the number of rows on which each
#' pairwise estimate is based; `summary()` returns a table describing how
#' each column is modelled.
#'
#' A pairwise correlation is treated as *unidentified* when the two variables
#' are observed together in fewer than 3 rows, or when its estimate is not
#' finite (a numerical safeguard). Such
#' entries are not estimated from data but completed as described in the
#' "Estimation details" section of [fit_synthesizer()]; `n_pair` alone
#' therefore does not identify all of them.
#'
#' @param object A fitted `sp_synthesizer`.
#' @param stratum Optional stratum label for stratified models (as shown by
#'   `summary()`); defaults to the pooled model.
#' @param include_missing Include the missingness indicators in the matrix?
#' @param what `"correlation"` (default) or `"n_pair"` (an integer matrix of
#'   jointly observed rows; the diagonal is the number of observed values).
#' @param ... Unused.
#' @return `copula_correlation()`: a correlation matrix (or a count matrix
#'   for `what = "n_pair"`). `summary()`: an
#'   object of class `summary.sp_synthesizer` (a list with a `columns` data
#'   frame).
#' @examples
#' fit <- fit_synthesizer(airquality)
#' round(copula_correlation(fit), 2)
#' summary(fit)
#' @export
copula_correlation <- function(object, stratum = NULL, include_missing = FALSE,
                               what = c("correlation", "n_pair")) {
  what <- match.arg(what)
  if (!inherits(object, "sp_synthesizer")) {
    stop("'object' must be an sp_synthesizer.", call. = FALSE)
  }
  mod <- object$pooled
  if (!is.null(stratum)) {
    if (is.null(object$strata)) stop("The model is not stratified.", call. = FALSE)
    i <- match(stratum, .stratum_labels(object$strata$values))
    if (is.na(i)) stop("Unknown stratum: ", stratum, call. = FALSE)
    mod <- object$strata$models[[i]] %||% object$pooled
  }
  R <- if (what == "n_pair") mod$copula$n_pair else mod$copula$R
  if (!include_missing) {
    keep <- !startsWith(colnames(R), ".missing_")
    R <- R[keep, keep, drop = FALSE]
  }
  R
}


#' @rdname copula_correlation
#' @export
summary.sp_synthesizer <- function(object, ...) {
  m <- object$pooled$marginals
  cols <- data.frame(
    column = names(m),
    class = vapply(m, `[[`, character(1), "class"),
    model = vapply(m, `[[`, character(1), "type"),
    support = vapply(m, function(z) {
      switch(z$type,
        continuous = sprintf("[%s, %s]", format(signif(min(z$quantiles), 4)),
                             format(signif(max(z$quantiles), 4))),
        discrete = , categorical = sprintf("%d values", length(z$values)),
        empty = "-")
    }, character(1)),
    missing = vapply(m, function(z) z$na_rate, numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  strata <- NULL
  if (!is.null(object$strata)) {
    st <- object$strata
    strata <- data.frame(stratum = .stratum_labels(st$values), n = st$counts,
                         model = ifelse(vapply(st$models, is.null, logical(1)),
                                        "pooled", "own"),
                         stringsAsFactors = FALSE)
  }
  structure(list(columns = cols, strata = strata, settings = object$settings,
                 n_train = object$n_train,
                 df = object$pooled$copula$df,
                 df_estimated = object$pooled$copula$df_estimated,
                 unidentified = object$pooled$copula$unidentified %||% 0L,
                 completion = object$pooled$copula$completion %||% "none",
                 pd_adjustment = object$pooled$copula$pd_adjustment %||% 0),
            class = "summary.sp_synthesizer")
}

#' @export
print.summary.sp_synthesizer <- function(x, ...) {
  cat("Synthesizer trained on", x$n_train, "rows\n")
  cat("Copula:", x$settings$copula)
  if (is.finite(x$df)) {
    cat(sprintf(" (df = %g%s)", x$df,
                if (identical(x$df_estimated, FALSE)) ", fixed: too few complete rows" else ""))
  }
  cat(" | missing data:", x$settings$missing, "\n")
  if (isTRUE(x$unidentified > 0L)) {
    cat(sprintf("Unidentified pairwise correlations (pooled model): %d (%s)\n", x$unidentified,
                if (identical(x$completion, "zero_fallback")) {
                  "no consistent completion; set to 0"
                } else {
                  "filled by maximum-determinant completion"
                }))
  }
  if (isTRUE(x$pd_adjustment > 1e-8)) {
    cat(sprintf("Positive-definiteness repair changed correlations by up to %.3g\n",
                x$pd_adjustment))
  }
  cat("\n")
  cols <- x$columns
  cols$missing <- .fmt_pct(cols$missing)
  print(cols, row.names = FALSE)
  if (!is.null(x$strata)) {
    cat("\nStrata:\n")
    print(x$strata, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.sp_synthesizer <- function(x, ...) {
  m <- x$pooled$marginals
  types <- table(vapply(m, `[[`, character(1), "type"))
  cat("<sp_synthesizer>\n")
  cat(sprintf("  trained on : %d rows, %d modelled column(s)\n",
              x$n_train, length(x$model_cols)))
  cat("  marginals  :", paste(sprintf("%d %s", as.integer(types), names(types)),
                              collapse = ", "), "\n")
  cop <- x$pooled$copula
  cat("  copula     :", cop$family,
      if (is.finite(cop$df)) sprintf("(df = %g)", cop$df) else "", "\n")
  cat("  missing    :", x$settings$missing, "\n")
  if (!is.null(x$settings$nonlinear_target)) {
    cat("  nonlinear  :", toupper(x$settings$nonlinear_engine %||% "cart"),
        "conditional target", x$settings$nonlinear_target, "\n")
    cat("  options    : residual=", x$settings$nonlinear_residual %||% "leaf",
        ", calibration=", x$settings$nonlinear_calibration %||% "quantile", "\n", sep = "")
  }
  if (!is.null(x$settings$nominal_target)) {
    cat("  nominal    : multinomial conditional target", x$settings$nominal_target, "\n")
  }
  if (!is.null(x$strata)) {
    own <- sum(!vapply(x$strata$models, is.null, logical(1)))
    cat(sprintf("  strata     : %d by %s (%d with own model)\n",
                length(x$strata$keys), paste(x$by, collapse = " x "), own))
    fb <- x$strata$conditional_fallback
    if (!is.null(fb) && nrow(fb)) {
      cat(sprintf("  fallback   : %d stratum model(s) generate the target by copula\n",
                  nrow(fb)))
    }
  }
  if (length(x$id_cols)) cat("  id columns :", paste(x$id_cols, collapse = ", "), "\n")
  invisible(x)
}
