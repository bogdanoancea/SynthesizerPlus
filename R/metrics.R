# Utility and disclosure-risk metrics ----------------------------------------

.check_pair <- function(real, synthetic) {
  if (is.matrix(real)) real <- as.data.frame(real)
  if (is.matrix(synthetic)) synthetic <- as.data.frame(synthetic)
  if (is.atomic(real) && is.null(dim(real))) real <- data.frame(x = real)
  if (is.atomic(synthetic) && is.null(dim(synthetic))) synthetic <- data.frame(x = synthetic)
  if (!is.data.frame(real) || !is.data.frame(synthetic)) {
    stop("'real' and 'synthetic' must be data frames.", call. = FALSE)
  }
  common <- intersect(names(real), names(synthetic))
  if (!length(common)) stop("'real' and 'synthetic' have no columns in common.", call. = FALSE)
  if (nrow(real) < 2L || nrow(synthetic) < 2L) {
    stop("Both data sets need at least 2 rows.", call. = FALSE)
  }
  list(real = .as_plain_df(real[common]), synthetic = .as_plain_df(synthetic[common]))
}

.select_vars <- function(df, vars) {
  if (is.null(vars)) return(names(df))
  .check_cols(vars, names(df), "vars")
  vars
}

# Marginal metrics -----------------------------------------------------------

.ks_stat <- function(a, b) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  t <- sort(unique(c(a, b)))
  max(abs(stats::ecdf(a)(t) - stats::ecdf(b)(t)))
}

.wasserstein <- function(a, b, grid = 999L) {
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  p <- seq_len(grid) / (grid + 1)
  w <- mean(abs(stats::quantile(a, p, names = FALSE) - stats::quantile(b, p, names = FALSE)))
  s <- stats::sd(a)
  if (!is.finite(s) || s == 0) s <- 1
  w / s
}

.tvd <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  lev <- unique(c(a, b))
  lev <- lev[!is.na(lev)]
  if (!length(lev)) return(NA_real_)
  pa <- tabulate(match(a, lev), length(lev))
  pb <- tabulate(match(b, lev), length(lev))
  if (sum(pa) == 0 || sum(pb) == 0) return(NA_real_)
  0.5 * sum(abs(pa / sum(pa) - pb / sum(pb)))
}

#' Compare marginal distributions of real and synthetic data
#'
#' For each variable, `marginal_metrics()` reports
#' * for numeric variables (incl. dates): the Kolmogorov--Smirnov distance
#'   between the empirical CDFs and the Wasserstein-1 distance scaled by the
#'   standard deviation of the real data;
#' * for categorical variables: the total variation distance between the
#'   category proportions;
#' * means, standard deviations and missing-value rates of both data sets.
#'
#' All distances are 0 for identical distributions; KS and TVD are bounded
#' by 1.
#'
#' @param real,synthetic Data frames with common columns.
#' @param vars Optional subset of variables.
#' @return A data frame with one row per variable.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' marginal_metrics(iris, syn)
#' @export
marginal_metrics <- function(real, synthetic, vars = NULL) {
  p <- .check_pair(real, synthetic)
  vars <- .select_vars(p$real, vars)
  rows <- lapply(vars, function(v) {
    a <- p$real[[v]]
    b <- p$synthetic[[v]]
    num <- .is_numeric_like(a)
    an <- if (num) as.numeric(a) else NULL
    bn <- if (num) as.numeric(b) else NULL
    data.frame(
      variable = v,
      type = if (num) "numeric" else "categorical",
      statistic = if (num) "KS" else "TVD",
      distance = if (num) .ks_stat(an, bn) else .tvd(a, b),
      wasserstein = if (num) .wasserstein(an, bn) else NA_real_,
      mean_real = if (num) mean(an, na.rm = TRUE) else NA_real_,
      mean_synthetic = if (num) mean(bn, na.rm = TRUE) else NA_real_,
      sd_real = if (num) stats::sd(an, na.rm = TRUE) else NA_real_,
      sd_synthetic = if (num) stats::sd(bn, na.rm = TRUE) else NA_real_,
      missing_real = mean(is.na(a)),
      missing_synthetic = mean(is.na(b)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Association ------------------------------------------------------------------

#' Mixed-type association matrix
#'
#' Computes pairwise association for data with numeric and categorical
#' variables: Spearman correlation for numeric--numeric pairs, bias-corrected
#' Cramér's V for categorical--categorical pairs, and the correlation ratio
#' \eqn{\eta} (on ranks) for numeric--categorical pairs. Pairwise complete
#' observations are used.
#'
#' @param data A data frame.
#' @param vars Optional subset of variables.
#' @return A symmetric matrix with ones on the diagonal.
#' @references Bergsma, W. (2013). A bias-correction for Cramér's V and
#'   Tschuprow's T. *Journal of the Korean Statistical Society*, 42(3),
#'   323--328. \doi{10.1016/j.jkss.2012.10.002}
#' @examples
#' round(association_matrix(iris), 2)
#' @export
association_matrix <- function(data, vars = NULL) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  vars <- .select_vars(data, vars)
  d <- length(vars)
  num <- vapply(data[vars], .is_numeric_like, logical(1))
  M <- diag(d)
  dimnames(M) <- list(vars, vars)
  if (d < 2L) return(M)
  for (i in seq_len(d - 1L)) {
    for (j in (i + 1L):d) {
      a <- data[[vars[i]]]
      b <- data[[vars[j]]]
      ok <- !is.na(a) & !is.na(b)
      val <- if (sum(ok) < 3L) {
        NA_real_
      } else if (num[i] && num[j]) {
        suppressWarnings(stats::cor(as.numeric(a[ok]), as.numeric(b[ok]), method = "spearman"))
      } else if (!num[i] && !num[j]) {
        .cramers_v(a[ok], b[ok])
      } else if (num[i]) {
        .corr_ratio(as.numeric(a[ok]), b[ok])
      } else {
        .corr_ratio(as.numeric(b[ok]), a[ok])
      }
      M[i, j] <- M[j, i] <- if (is.finite(val)) val else NA_real_
    }
  }
  M
}

.cramers_v <- function(a, b) {
  tab <- table(as.character(a), as.character(b))
  r <- nrow(tab)
  k <- ncol(tab)
  n <- sum(tab)
  if (r < 2L || k < 2L) return(0)
  e <- outer(rowSums(tab), colSums(tab)) / n
  chi2 <- sum((tab - e)^2 / e)
  phi2 <- max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
  rc <- r - (r - 1)^2 / (n - 1)
  kc <- k - (k - 1)^2 / (n - 1)
  den <- min(kc - 1, rc - 1)
  if (den <= 0) return(0)
  sqrt(phi2 / den)
}

.corr_ratio <- function(x, g) {
  x <- rank(x)
  g <- as.character(g)
  if (length(unique(g)) < 2L) return(0)
  m <- tapply(x, g, mean)
  nk <- tapply(x, g, length)
  ss_b <- sum(nk * (m - mean(x))^2)
  ss_t <- sum((x - mean(x))^2)
  if (ss_t == 0) return(0)
  sqrt(ss_b / ss_t)
}

# Feature preparation for classifiers ---------------------------------------

# Feature preparation for the classifiers, in two steps so that it can be
# learned on training data only (no leakage into cross-validation folds):
# .prep_fit() learns medians, means, standard deviations and the retained
# categories; .prep_apply() turns any data frame into the design matrix.
.prep_fit <- function(df, max_levels = 20L) {
  spec <- list()
  for (nm in names(df)) {
    v <- df[[nm]]
    if (.is_numeric_like(v)) {
      v <- as.numeric(v)
      miss <- is.na(v)
      if (all(miss)) next
      med <- stats::median(v, na.rm = TRUE)
      v[miss] <- med
      s <- stats::sd(v)
      spec[[nm]] <- list(type = "numeric", median = med, na_indicator = any(miss),
                         mean = mean(v), sd = if (is.finite(s)) s else 0)
    } else {
      v <- as.character(v)
      v[is.na(v)] <- "(missing)"
      tab <- sort(table(v), decreasing = TRUE)
      lumped <- length(tab) > max_levels
      if (lumped) {
        keep <- names(tab)[seq_len(max_levels - 1L)]
        v[!v %in% keep] <- "(other)"
        tab <- sort(table(v), decreasing = TRUE)
      }
      if (length(tab) < 2L) next
      spec[[nm]] <- list(type = "categorical", levels = names(tab), lumped = lumped)
    }
  }
  spec
}

.prep_apply <- function(spec, df) {
  cols <- list()
  numeric_cols <- character(0)
  for (nm in names(spec)) {
    sp <- spec[[nm]]
    v <- df[[nm]]
    if (sp$type == "numeric") {
      v <- as.numeric(v)
      miss <- is.na(v)
      v[miss] <- sp$median
      if (sp$na_indicator) cols[[paste0(nm, "..NA")]] <- as.numeric(miss)
      if (sp$sd > 0) {
        cols[[nm]] <- (v - sp$mean) / sp$sd
        numeric_cols <- c(numeric_cols, nm)
      }
    } else {
      v <- as.character(v)
      v[is.na(v)] <- "(missing)"
      if (sp$lumped) v[!v %in% sp$levels] <- "(other)"
      for (l in sp$levels[-1L]) cols[[paste0(nm, "==", l)]] <- as.numeric(v == l)
    }
  }
  if (!length(cols)) {
    X <- matrix(numeric(0), NROW(df), 0L)
  } else {
    X <- do.call(cbind, cols)
    colnames(X) <- names(cols)
  }
  attr(X, "numeric_cols") <- numeric_cols
  X
}

.prep_features <- function(df, max_levels = 20L) {
  .prep_apply(.prep_fit(df, max_levels), df)
}

# Add squares and pairwise products of the numeric features (detects
# differences in spread and in correlation, not only in location)
.add_quadratic <- function(X, max_numeric = 20L) {
  nc <- attr(X, "numeric_cols")
  if (length(nc) > max_numeric) nc <- nc[seq_len(max_numeric)]
  if (!length(nc)) return(X)
  extra <- list()
  for (i in seq_along(nc)) {
    extra[[paste0(nc[i], "^2")]] <- X[, nc[i]]^2
    if (i < length(nc)) {
      for (j in (i + 1L):length(nc)) {
        extra[[paste0(nc[i], ":", nc[j])]] <- X[, nc[i]] * X[, nc[j]]
      }
    }
  }
  cbind(X, do.call(cbind, extra))
}

.fit_logit <- function(X, y, ridge = 0) {
  Xi <- cbind(`(Intercept)` = 1, X)
  if (ridge > 0) return(.fit_logit_ridge(Xi, y, ridge))
  fit <- suppressWarnings(stats::glm.fit(Xi, y, family = stats::binomial()))
  b <- fit$coefficients
  b[is.na(b)] <- 0
  list(coef = b, fitted = fit$fitted.values, rank = fit$rank)
}

# Ridge-penalised logistic regression by iteratively reweighted least squares
# (intercept not penalised)
.fit_logit_ridge <- function(Xi, y, lambda, maxit = 50L, tol = 1e-8) {
  p <- ncol(Xi)
  pen <- c(0, rep(lambda, p - 1L))
  b <- c(stats::qlogis(min(max(mean(y), 1e-6), 1 - 1e-6)), rep(0, p - 1L))
  for (it in seq_len(maxit)) {
    eta <- drop(Xi %*% b)
    mu <- stats::plogis(eta)
    w <- pmax(mu * (1 - mu), 1e-10)
    z <- eta + (y - mu) / w
    b_new <- solve(crossprod(Xi, w * Xi) + diag(pen, p), crossprod(Xi, w * z))
    b_new <- drop(b_new)
    if (max(abs(b_new - b)) < tol) {
      b <- b_new
      break
    }
    b <- b_new
  }
  list(coef = b, fitted = stats::plogis(drop(Xi %*% b)), rank = p)
}

.auc <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.stack <- function(real, synthetic, vars, max_rows, model = "linear") {
  r <- .subsample(real[vars], max_rows)
  s <- .subsample(synthetic[vars], max_rows)
  X <- .prep_features(rbind(r, s))
  if (model == "quadratic") X <- .add_quadratic(X)
  list(X = X, y = c(rep(0, nrow(r)), rep(1, nrow(s))))
}

# Stacked raw data (no preprocessing), for cross-validation
.stack_raw <- function(real, synthetic, vars, max_rows) {
  r <- .subsample(real[vars], max_rows)
  s <- .subsample(synthetic[vars], max_rows)
  list(data = rbind(r, s), y = c(rep(0, nrow(r)), rep(1, nrow(s))))
}

#' Propensity-score utility (pMSE)
#'
#' The real and synthetic records are stacked and a logistic regression
#' (main effects; categorical variables as dummies, missing values as
#' indicators) predicts membership of the synthetic set. The propensity
#' score mean-squared error
#' \deqn{pMSE = \frac{1}{N}\sum_i (\hat p_i - c)^2,\quad c = n_{syn}/N}
#' is 0 when the two sets are indistinguishable. `pmse_ratio` divides it by
#' its expectation under the null hypothesis that both sets come from the
#' same distribution, \eqn{(k-1)(1-c)^2 c / N} with \eqn{k} model
#' parameters; values near 1 indicate high utility.
#'
#' With `model = "quadratic"`, squares and pairwise products of the numeric
#' variables are added to the logistic model (as recommended by Snoke et
#' al. for detecting differences in the dependence structure).
#'
#' @param real,synthetic Data frames with common columns.
#' @param vars Optional subset of variables.
#' @param max_rows Each data set is subsampled to at most this many rows.
#' @param seed Optional random seed (subsampling).
#' @param model Propensity model: `"linear"` (main effects, default) or
#'   `"quadratic"`.
#' @return A list with `pmse`, `pmse_null`, `pmse_ratio` and `n_params`.
#' @references Snoke, J., Raab, G. M., Nowok, B., Dibben, C. and
#'   Slavkovic, A. (2018). General and specific utility measures for
#'   synthetic data. *Journal of the Royal Statistical Society A*, 181(3),
#'   663--688. \doi{10.1111/rssa.12358}
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' pmse(iris, syn)
#' @export
pmse <- function(real, synthetic, vars = NULL, max_rows = 10000L, seed = NULL,
                 model = c("linear", "quadratic")) {
  model <- match.arg(model)
  p <- .check_pair(real, synthetic)
  vars <- .select_vars(p$real, vars)
  .with_seed(seed, {
    st <- .stack(p$real, p$synthetic, vars, max_rows, model)
    fit <- .fit_logit(st$X, st$y)
    N <- length(st$y)
    c0 <- mean(st$y)
    val <- mean((fit$fitted - c0)^2)
    null <- (fit$rank - 1) * (1 - c0)^2 * c0 / N
    list(pmse = val, pmse_null = null,
         pmse_ratio = if (null > 0) val / null else NA_real_,
         n_params = fit$rank)
  })
}

#' Discriminator AUC
#'
#' Cross-validated area under the ROC curve of a logistic-regression
#' classifier trained to tell real from synthetic records. Out-of-fold
#' predictions are used, and all preprocessing (imputation, scaling, choice
#' of categories) is learned within the training folds, so the AUC is not
#' optimistically biased. Values
#' near 0.5 mean the classifier cannot distinguish the data sets; values
#' near 1 mean the synthetic data are easy to spot.
#'
#' Values clearly *below* 0.5 are a warning sign of a different kind: they
#' arise when synthetic records are copies of real ones. A held-out record
#' then has an identical twin with the opposite label in the training folds,
#' and the classifier systematically predicts the wrong class. Check such
#' cases with [dcr()]. The AUC is therefore deliberately not folded into a
#' symmetric measure such as \eqn{\max(AUC, 1 - AUC)}: that would report
#' copies of the real data as easily distinguishable.
#'
#' The default `model = "quadratic"` adds squares and pairwise products of
#' the numeric variables (ridge-penalised), so that differences in spread and
#' in correlation are detected, not only differences in location.
#' `model = "linear"` uses main effects only.
#'
#' @inheritParams pmse
#' @param folds Number of cross-validation folds.
#' @param model Classifier features: `"quadratic"` (default) or `"linear"`.
#' @return The cross-validated AUC (a number).
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' discriminator_auc(iris, syn, seed = 1)
#' @export
discriminator_auc <- function(real, synthetic, vars = NULL, folds = 5L,
                              max_rows = 10000L, seed = NULL,
                              model = c("quadratic", "linear")) {
  model <- match.arg(model)
  p <- .check_pair(real, synthetic)
  vars <- .select_vars(p$real, vars)
  folds <- .check_count(folds, "folds", allow_zero = FALSE)
  if (folds < 2L) stop("'folds' must be at least 2.", call. = FALSE)
  .with_seed(seed, {
    st <- .stack_raw(p$real, p$synthetic, vars, max_rows)
    N <- length(st$y)
    ridge <- if (model == "quadratic") 1 else 0
    fold <- sample(rep_len(seq_len(folds), N))
    score <- numeric(N)
    design <- function(spec, d) {
      X <- .prep_apply(spec, d)
      if (model == "quadratic") X <- .add_quadratic(X)
      X
    }
    for (k in seq_len(folds)) {
      te <- fold == k
      # preprocessing is learned on the training folds only
      spec <- .prep_fit(st$data[!te, , drop = FALSE])
      Xtr <- design(spec, st$data[!te, , drop = FALSE])
      Xte <- design(spec, st$data[te, , drop = FALSE])
      fit <- .fit_logit(Xtr, st$y[!te], ridge = ridge)
      score[te] <- drop(cbind(1, Xte) %*% fit$coef)
    }
    .auc(st$y, score)
  })
}

#' Confidence-interval overlap of regression coefficients
#'
#' A specific utility measure: the same regression model is fitted to the
#' real and the synthetic data and, for each coefficient, the overlap of the
#' two confidence intervals is computed (Karr et al., 2006):
#' \deqn{IO = \frac{1}{2}\left(\frac{U-L}{u_r-l_r} + \frac{U-L}{u_s-l_s}\right)}
#' where \eqn{[L, U]} is the intersection of the intervals. Overlap is 1 for
#' identical intervals and can become negative when they are disjoint.
#'
#' @param real,synthetic Data frames.
#' @param formula Model formula.
#' @param family A [stats::family()] object (default Gaussian, i.e. linear
#'   regression).
#' @param level Confidence level.
#' @return A data frame with estimates, intervals, the standardised
#'   difference of the estimates and the interval overlap per coefficient.
#' @references Karr, A. F., Kohnen, C. N., Oganian, A., Reiter, J. P. and
#'   Sanil, A. P. (2006). A framework for evaluating the utility of data
#'   altered to protect confidentiality. *The American Statistician*, 60(3),
#'   224--232. \doi{10.1198/000313006X124640}
#' @examples
#' syn <- synthesize(mtcars, seed = 1)
#' ci_overlap(mtcars, syn, mpg ~ wt + hp)
#' @export
ci_overlap <- function(real, synthetic, formula, family = stats::gaussian(),
                       level = 0.95) {
  fr <- stats::glm(formula, data = real, family = family)
  fs <- stats::glm(formula, data = synthetic, family = family)
  cr <- stats::coef(summary(fr))
  cs <- stats::coef(summary(fs))
  terms <- intersect(rownames(cr), rownames(cs))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- lapply(terms, function(tm) {
    br <- cr[tm, 1L]
    sr <- cr[tm, 2L]
    bs <- cs[tm, 1L]
    ss <- cs[tm, 2L]
    lr <- br - z * sr
    ur <- br + z * sr
    ls <- bs - z * ss
    us <- bs + z * ss
    L <- max(lr, ls)
    U <- min(ur, us)
    data.frame(term = tm, estimate_real = br, estimate_synthetic = bs,
               lower_real = lr, upper_real = ur, lower_synthetic = ls,
               upper_synthetic = us, std_diff = (bs - br) / sr,
               overlap = 0.5 * ((U - L) / (ur - lr) + (U - L) / (us - ls)),
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# Disclosure risk --------------------------------------------------------------

.gower_prep <- function(df, ref) {
  lapply(names(df), function(nm) {
    v <- df[[nm]]
    r <- ref[[nm]]
    if (.is_numeric_like(r)) {
      rv <- as.numeric(r)
      rng <- if (any(!is.na(rv))) diff(range(rv, na.rm = TRUE)) else NA_real_
      if (!is.finite(rng) || rng == 0) rng <- 1
      list(num = TRUE, x = as.numeric(v) / rng)
    } else {
      list(num = FALSE, x = as.character(v))
    }
  })
}

# Minimum Gower distance from each row of A to the rows of B.
# na = "category": a missing value is treated as a category of its own
#   (missing vs missing = 0, missing vs observed = 1); every variable counts.
# na = "exclude": standard Gower (1971); a variable missing in either record
#   is left out of that comparison, and the distance is averaged over the
#   variables observed in both (NA if there are none).
.min_gower <- function(A, B, exclude_self = FALSE, block = 256L,
                       na = c("category", "exclude")) {
  na <- match.arg(na)
  nA <- length(A[[1L]]$x)
  nB <- length(B[[1L]]$x)
  p <- length(A)
  out <- numeric(nA)
  compared <- numeric(nA)
  for (s in seq(1L, nA, by = block)) {
    idx <- s:min(nA, s + block - 1L)
    D <- matrix(0, length(idx), nB)
    W <- matrix(0, length(idx), nB)
    for (j in seq_len(p)) {
      a <- A[[j]]$x[idx]
      b <- B[[j]]$x
      dj <- if (A[[j]]$num) pmin(abs(outer(a, b, "-")), 1) else (outer(a, b, "!=") + 0)
      na_a <- is.na(a)
      na_b <- is.na(b)
      if (na == "category") {
        dj[outer(na_a, na_b, "!=")] <- 1
        dj[outer(na_a, na_b, "&")] <- 0
        D <- D + dj
        W <- W + 1
      } else {
        comparable <- outer(!na_a, !na_b, "&")
        dj[!comparable] <- 0
        D <- D + dj
        W <- W + comparable
      }
    }
    D <- D / W
    D[W == 0] <- NA_real_
    if (exclude_self) D[cbind(seq_along(idx), idx)] <- Inf
    for (i in seq_along(idx)) {
      r <- D[i, ]
      if (all(is.na(r))) {
        out[idx[i]] <- NA_real_
        compared[idx[i]] <- 0
      } else {
        m <- min(r, na.rm = TRUE)
        out[idx[i]] <- m
        # among equally close records, report the one compared on most variables
        compared[idx[i]] <- max(W[i, which(r == m)]) / p
      }
    }
  }
  attr(out, "compared") <- compared
  out
}

#' Distance to closest record (disclosure risk)
#'
#' For every synthetic record (in a subsample of at most `max_rows`), the
#' Gower distance to the closest real record is computed. As a baseline the
#' same is done for real records with respect to the *other* real records.
#' Synthetic records that are much closer to real records than real records
#' are to each other, and in particular exact copies (distance 0), indicate
#' a disclosure risk.
#'
#' Numeric variables contribute their absolute difference divided by the
#' range in the real data, categorical variables 0 (equal) or 1 (different).
#' Missing values are handled according to `na`:
#'
#' * `"category"` (default): missingness is treated as a category of its own
#'   -- two missing values match, a missing and an observed value differ --
#'   and the distance is averaged over all variables. Because the
#'   synthesizer reproduces missingness patterns, and an intruder sees them,
#'   this is the appropriate view for disclosure risk.
#' * `"exclude"`: the standard Gower (1971) coefficient; a variable missing in
#'   either record is left out of the comparison and the distance is averaged
#'   over the variables observed in both. Records with no variable in common
#'   have no distance (`NA`) and are ignored.
#'
#' @inheritParams pmse
#' @param max_rows Maximum number of records used from each data set.
#' @param close_quantile Quantile of the real-to-real distances that defines
#'   a synthetic record as "too close" to a real record (see `close_share`).
#' @param na Treatment of missing values: `"category"` (default) or
#'   `"exclude"` (standard Gower); see Details.
#' @return A list with the distance vectors `synthetic` and `real`,
#'   their `quantiles`, the share of synthetic records that are exact copies
#'   of a real record on *all* variables (`exact_match_rate`), the share at
#'   distance zero on the variables compared (`zero_distance_rate`; with
#'   `na = "exclude"` this can include records that agree with a real record
#'   only on the few variables observed in both), the average share of
#'   variables on which the nearest real record was compared
#'   (`compared_share`; always 1 with `na = "category"`), the share of exact
#'   duplicates within the real data (`real_duplicate_rate`) and `ratio`,
#'   the median synthetic distance
#'   divided by the median real-to-real distance (values well below 1 are a
#'   warning sign), and `close_share`, the share of synthetic records that
#'   are closer to a real record than `close_quantile` (5\%) of real records
#'   are to their nearest real neighbour. For synthetic records that are as
#'   far from real ones as real records are from each other, `close_share`
#'   is about `close_quantile`; larger values indicate records that may be
#'   recognised as near-copies of real persons.
#' @references Gower, J. C. (1971). A general coefficient of similarity and
#'   some of its properties. *Biometrics*, 27(4), 857--871.
#'   \doi{10.2307/2528823}
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' str(dcr(iris, syn, seed = 1)[c("exact_match_rate", "ratio")])
#'
#' # data with missing values: missingness as a category, or standard Gower
#' aq <- synthesize(airquality, seed = 1)
#' dcr(airquality, aq, seed = 1)$ratio
#' dcr(airquality, aq, seed = 1, na = "exclude")$ratio
#' @export
dcr <- function(real, synthetic, vars = NULL, max_rows = 2000L, seed = NULL,
                close_quantile = 0.05, na = c("category", "exclude")) {
  close_quantile <- .check_prob(close_quantile, "close_quantile")
  na <- match.arg(na)
  p <- .check_pair(real, synthetic)
  vars <- .select_vars(p$real, vars)
  .with_seed(seed, {
    r <- .subsample(p$real[vars], max_rows)
    s <- .subsample(p$synthetic[vars], max_rows)
    R <- .gower_prep(r, p$real)
    S <- .gower_prep(s, p$real)
    d_syn <- .min_gower(S, R, na = na)
    d_real <- .min_gower(R, R, exclude_self = TRUE, na = na)
    cmp_syn <- attr(d_syn, "compared")[!is.na(d_syn)]
    cmp_real <- attr(d_real, "compared")[!is.na(d_real)]
    d_syn <- as.numeric(d_syn[!is.na(d_syn)])
    d_real <- as.numeric(d_real[!is.na(d_real)])
    if (!length(d_syn) || !length(d_real)) {
      stop("No pair of records has a variable observed in both.", call. = FALSE)
    }
    probs <- c(0.01, 0.05, 0.25, 0.5)
    q <- rbind(synthetic = stats::quantile(d_syn, probs, names = FALSE),
               real = stats::quantile(d_real, probs, names = FALSE))
    colnames(q) <- paste0("q", probs * 100)
    med_r <- stats::median(d_real)
    list(synthetic = d_syn, real = d_real, quantiles = q,
         exact_match_rate = mean(d_syn < 1e-12 & cmp_syn == 1),
         zero_distance_rate = mean(d_syn < 1e-12),
         compared_share = mean(cmp_syn),
         close_share = mean(d_syn < stats::quantile(d_real, close_quantile, names = FALSE)),
         real_duplicate_rate = mean(d_real < 1e-12 & cmp_real == 1),
         ratio = if (med_r > 0) stats::median(d_syn) / med_r else NA_real_)
  })
}

# Full comparison --------------------------------------------------------------

#' Compare synthetic data with the real data
#'
#' @description
#' `compare_synthetic()` evaluates the fidelity, utility and disclosure risk
#' of a synthetic data set in one call:
#'
#' * **marginal**: [marginal_metrics()] per variable;
#' * **association**: [association_matrix()] of both data sets and their
#'   difference;
#' * **utility**: [pmse()] and cross-validated [discriminator_auc()];
#' * **privacy**: [dcr()] distance to closest record.
#'
#' The result has `print()`, `summary()`, `plot()` and
#' `ggplot2::autoplot()` methods.
#'
#' @param real,synthetic Data frames (or objects coercible to data frames)
#'   with common columns.
#' @param vars Optional subset of variables.
#' @param metrics Which groups of metrics to compute.
#' @param max_rows Maximum rows per data set for the classifier and distance
#'   computations.
#' @param seed Optional random seed.
#' @return An object of class `sp_comparison`.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' cmp <- compare_synthetic(iris, syn, seed = 1)
#' cmp
#' plot(cmp, type = "association")
#' @export
compare_synthetic <- function(real, synthetic, vars = NULL,
                              metrics = c("marginal", "association", "utility", "privacy"),
                              max_rows = 5000L, seed = NULL) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  p <- .check_pair(real, synthetic)
  vars <- .select_vars(p$real, vars)
  .with_seed(seed, {
    out <- list(vars = vars, n_real = nrow(p$real), n_synthetic = nrow(p$synthetic),
                real = .subsample(p$real[vars], max_rows),
                synthetic = .subsample(p$synthetic[vars], max_rows))
    if ("marginal" %in% metrics) out$marginal <- marginal_metrics(p$real, p$synthetic, vars)
    if ("association" %in% metrics && length(vars) > 1L) {
      ar <- association_matrix(out$real)
      as <- association_matrix(out$synthetic)
      out$association <- list(real = ar, synthetic = as, difference = as - ar,
                              mean_abs_diff = mean(abs((as - ar)[upper.tri(ar)]), na.rm = TRUE))
    }
    if ("utility" %in% metrics) {
      out$utility <- c(pmse(p$real, p$synthetic, vars, max_rows = max_rows),
                       list(auc = discriminator_auc(p$real, p$synthetic, vars, max_rows = max_rows)))
    }
    if ("privacy" %in% metrics) {
      out$privacy <- dcr(p$real, p$synthetic, vars, max_rows = min(max_rows, 2000L))
    }
    structure(out, class = "sp_comparison")
  })
}

#' @export
print.sp_comparison <- function(x, digits = 3, ...) {
  cat("<sp_comparison> real:", x$n_real, "rows | synthetic:", x$n_synthetic,
      "rows |", length(x$vars), "variables\n")
  if (!is.null(x$marginal)) {
    cat("\nMarginal distributions (KS for numeric, TVD for categorical):\n")
    m <- x$marginal[, c("variable", "statistic", "distance", "wasserstein",
                        "missing_real", "missing_synthetic")]
    num <- vapply(m, is.numeric, logical(1))
    m[num] <- lapply(m[num], round, digits)
    print(m, row.names = FALSE)
  }
  if (!is.null(x$association)) {
    cat(sprintf("\nAssociation: mean |difference| = %.*f\n", digits,
                x$association$mean_abs_diff))
  }
  if (!is.null(x$utility)) {
    cat(sprintf("Utility: pMSE = %.*g (ratio to null = %.*f), discriminator AUC = %.*f\n",
                digits, x$utility$pmse, 2, x$utility$pmse_ratio, digits, x$utility$auc))
    if (is.finite(x$utility$auc) && x$utility$auc < 0.4) {
      cat("  Note: an AUC well below 0.5 usually means that synthetic records copy real\n",
          "  ones (see ?discriminator_auc); check the privacy metrics.\n", sep = "")
    }
  }
  if (!is.null(x$privacy)) {
    cat(sprintf("Privacy: exact copies = %s, median DCR ratio (synthetic / real) = %.*f\n",
                .fmt_pct(x$privacy$exact_match_rate, 2), 2, x$privacy$ratio))
  }
  invisible(x)
}

#' @export
summary.sp_comparison <- function(object, ...) {
  res <- data.frame(
    metric = c("mean marginal distance", "max marginal distance",
               "mean |association difference|", "pMSE ratio",
               "discriminator AUC", "exact copy rate", "DCR ratio"),
    value = c(
      if (!is.null(object$marginal)) mean(object$marginal$distance, na.rm = TRUE) else NA,
      if (!is.null(object$marginal)) max(object$marginal$distance, na.rm = TRUE) else NA,
      object$association$mean_abs_diff %||% NA,
      object$utility$pmse_ratio %||% NA,
      object$utility$auc %||% NA,
      object$privacy$exact_match_rate %||% NA,
      object$privacy$ratio %||% NA
    ),
    ideal = c("0", "0", "0", "~1", "0.5", "0", ">= 1"),
    stringsAsFactors = FALSE
  )
  res
}
