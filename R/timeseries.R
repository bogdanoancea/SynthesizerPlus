# Time-series synthesizers ---------------------------------------------------

.ts_closeness_msg <- paste(
  "'closeness', 'perturb' and 'dependence' are not supported for time-series",
  "synthesizers and were ignored. Use them with tabular data, or perturb the",
  "generated series yourself.")


#' Fit a synthesizer to a (multivariate) time series
#'
#' @description
#' Methods for `ts` objects. Four generators are available:
#'
#' * `"copula_var"` (default): each series is mapped to normal scores through
#'   its empirical marginal distribution; a vector autoregression (order
#'   chosen by AIC) is fitted to the scores, and simulated paths are mapped
#'   back through the marginal quantile functions. This reproduces the
#'   marginal distributions *and* the auto- and cross-correlation structure.
#'   Innovations are resampled from the fitted residuals, so their joint
#'   distribution is preserved.
#' * `"block"`: moving (circular) block bootstrap of the observed series.
#' * `"stationary"`: stationary bootstrap (Politis and Romano, 1994) with
#'   geometrically distributed block lengths.
#' * `"iid"`: ignores serial dependence; rows are drawn from a cross-sectional
#'   copula model (see [fit_synthesizer()]).
#'
#' Missing values are linearly interpolated before fitting (with a warning).
#'
#' *Scope.* All four methods treat the series as (approximately)
#' **stationary**: the empirical marginal distribution and the dependence
#' structure are assumed constant over time. A trend, a level shift, changing
#' variance or strong deterministic seasonality is reproduced only
#' imperfectly (e.g. a trending series is resampled around its overall
#' distribution, and a seasonal pattern only as far as the VAR order or block
#' length covers it). For such data, model the deterministic part yourself
#' (detrend, difference or deseasonalise), synthesise the stationary
#' remainder, and add the deterministic part back. Arguments that the chosen
#' method does not use trigger a warning.
#'
#' @param data A `ts` object (univariate or multivariate).
#' @param method One of `"copula_var"`, `"block"`, `"stationary"`, `"iid"`.
#' @param max_lag Maximum VAR order considered by `"copula_var"` when the
#'   order is selected by AIC.
#' @param order Optional fixed VAR order for `"copula_var"` (skips AIC
#'   selection). Useful for seasonal data on short series, where AIC tends to
#'   choose an order below the seasonal period, e.g. `order = 12` for
#'   monthly data.
#' @param block_length Block length (`"block"`) or mean block length
#'   (`"stationary"`). Default: `ceiling(n^(1/3))`.
#' @param burn_in Number of burn-in steps discarded when simulating a VAR.
#' @param seed Optional random seed making the fit reproducible (see
#'   [fit_synthesizer()]).
#' @param ... Passed to [fit_synthesizer.data.frame()] when `method = "iid"`.
#'
#' @return An object of class `sp_ts_synthesizer`. Use [generate()] to draw
#'   synthetic series of any length; the result is a `ts` with the same
#'   frequency (and, by default, the same start) as `data`.
#'
#' @references Politis, D. N. and Romano, J. P. (1994). The stationary
#'   bootstrap. *Journal of the American Statistical Association*, 89(428),
#'   1303--1313. \doi{10.1080/01621459.1994.10476870}
#'
#' @examples
#' fit <- fit_synthesizer(ldeaths)
#' syn <- generate(fit, seed = 1)
#' plot_ts(ldeaths, syn)
#'
#' fit_mv <- fit_synthesizer(EuStockMarkets[1:500, ], method = "stationary")
#' head(generate(fit_mv, n = 100, seed = 1))
#' @export
fit_synthesizer.ts <- function(data,
                               method = c("copula_var", "block", "stationary", "iid"),
                               max_lag = 10L, order = NULL, block_length = NULL,
                               burn_in = 100L, seed = NULL, ...) {
  method <- match.arg(method)
  # arguments that the chosen method does not use
  ignored <- c(
    max_lag = !missing(max_lag) && (method != "copula_var" || !is.null(order)),
    order = !missing(order) && !is.null(order) && method != "copula_var",
    burn_in = !missing(burn_in) && method != "copula_var",
    block_length = !missing(block_length) && !is.null(block_length) &&
      !method %in% c("block", "stationary")
  )
  if (any(ignored)) {
    why <- c(max_lag = if (method == "copula_var") "it is only used when 'order' is NULL"
             else "it is only used by method = \"copula_var\"",
             order = "it is only used by method = \"copula_var\"",
             burn_in = "it is only used by method = \"copula_var\"",
             block_length = "it is only used by methods \"block\" and \"stationary\"")
    for (a in names(ignored)[ignored]) {
      warning(sprintf("'%s' is ignored for method = \"%s\": %s.", a, method, why[[a]]),
              call. = FALSE)
    }
  }
  if (!is.null(seed)) {
    return(.with_seed(seed, .fit_ts(data, method, max_lag, order, block_length,
                                    burn_in, ...)))
  }
  .fit_ts(data, method, max_lag, order, block_length, burn_in, ...)
}

.fit_ts <- function(data, method, max_lag, order, block_length, burn_in, ...) {
  max_lag <- .check_count(max_lag, "max_lag")
  if (!is.null(order)) order <- .check_count(order, "order")
  burn_in <- .check_count(burn_in, "burn_in")
  if (method != "iid") {
    .warn_unused_dots(list(...), "fit_synthesizer",
                      if (any(c("closeness", "perturb", "dependence") %in% names(list(...)))) {
                        "Closeness control is not available for time series."
                      })
  }
  X <- as.matrix(data)
  if (!is.numeric(X)) stop("Time series must be numeric.", call. = FALSE)
  if (is.null(colnames(X))) {
    colnames(X) <- if (ncol(X) == 1L) "Series" else paste0("Series", seq_len(ncol(X)))
  }
  N <- nrow(X)
  if (N < 10L) stop("At least 10 observations are needed.", call. = FALSE)
  if (anyNA(X)) {
    warning("Missing values in the time series were linearly interpolated.",
            call. = FALSE)
    X <- apply(X, 2L, function(v) {
      ok <- !is.na(v)
      if (sum(ok) < 2L) stop("A series has fewer than 2 observed values.", call. = FALSE)
      stats::approx(which(ok), v[ok], xout = seq_along(v), rule = 2)$y
    })
    X <- matrix(X, nrow = N)
  }
  colnames(X) <- colnames(X) %||% colnames(as.matrix(data))
  if (is.null(colnames(X))) colnames(X) <- paste0("Series", seq_len(ncol(X)))

  obj <- list(method = method, n_train = N, start = stats::start(data),
              frequency = stats::frequency(data), names = colnames(X),
              univariate = is.null(dim(data)), data = X)

  if (method %in% c("block", "stationary")) {
    obj$block_length <- if (is.null(block_length)) {
      max(1L, as.integer(ceiling(N^(1 / 3))))
    } else {
      min(.check_count(block_length, "block_length", allow_zero = FALSE), N)
    }
  } else if (method == "iid") {
    obj$cross <- fit_synthesizer.data.frame(as.data.frame(X), ...)
    obj$data <- NULL
  } else {
    obj$marginals <- lapply(seq_len(ncol(X)), function(j) {
      .fit_marginal(X[, j], colnames(X)[j], discrete_threshold = 0L)
    })
    Z <- vapply(seq_len(ncol(X)), function(j) {
      stats::qnorm(rank(X[, j], ties.method = "average") / (N + 1))
    }, numeric(N))
    Z <- matrix(Z, nrow = N)
    obj$var <- .fit_var(Z, max_lag, order)
    obj$burn_in <- burn_in
    # stationary standard deviation of the simulated scores
    sim <- .simulate_var(obj$var, 5000L, burn_in)
    obj$z_sd <- apply(sim, 2L, stats::sd)
    obj$data <- NULL
  }
  class(obj) <- "sp_ts_synthesizer"
  obj
}

.fit_var <- function(Z, max_lag, order = NULL) {
  d <- ncol(Z)
  N <- nrow(Z)
  aic <- is.null(order)
  if (aic) {
    max_lag <- min(max_lag, floor(N / (2 * d + 2)), N - 2L)
  } else {
    if (order * d >= N - order) {
      stop(sprintf("'order' = %d is too large for %d observations of %d series.",
                   order, N, d), call. = FALSE)
    }
    max_lag <- order
  }
  fit <- if (d == 1L) {
    stats::ar(Z[, 1L], aic = aic, order.max = max(max_lag, 0L), method = "yw",
              demean = TRUE)
  } else {
    stats::ar(Z, aic = aic, order.max = max(max_lag, 0L), method = "yw",
              demean = TRUE)
  }
  p <- fit$order
  A <- array(0, dim = c(max(p, 1L), d, d))
  if (p > 0L) {
    if (d == 1L) A[seq_len(p), 1L, 1L] <- fit$ar else A[seq_len(p), , ] <- fit$ar
  }
  res <- matrix(fit$resid, ncol = d)
  res <- res[stats::complete.cases(res), , drop = FALSE]
  if (p == 0L || nrow(res) < 2L) res <- sweep(Z, 2L, colMeans(Z))
  res <- sweep(res, 2L, colMeans(res))
  list(order = p, A = A, mean = colMeans(Z), resid = res,
       init = sweep(Z, 2L, colMeans(Z)))
}

.simulate_var <- function(v, n, burn_in) {
  d <- ncol(v$resid)
  p <- v$order
  total <- n + burn_in
  idx <- sample.int(nrow(v$resid), total, replace = TRUE)
  E <- v$resid[idx, , drop = FALSE]
  Y <- matrix(0, total + max(p, 1L), d)
  if (p > 0L) {
    s <- sample.int(nrow(v$init) - p + 1L, 1L)
    Y[seq_len(p), ] <- v$init[s:(s + p - 1L), , drop = FALSE]
  }
  off <- max(p, 1L)
  for (t in seq_len(total)) {
    y <- E[t, ]
    if (p > 0L) {
      for (k in seq_len(p)) y <- y + v$A[k, , ] %*% Y[off + t - k, ]
    }
    Y[off + t, ] <- y
  }
  Y <- Y[off + burn_in + seq_len(n), , drop = FALSE]
  sweep(Y, 2L, v$mean, "+")
}

#' @rdname generate
#' @param start Start time of the generated series (default: the start of
#'   the training series).
#' @export
generate.sp_ts_synthesizer <- function(object, n = object$n_train, seed = NULL,
                                       start = object$start, ...) {
  dots <- list(...)
  if (any(c("closeness", "perturb", "dependence") %in% names(dots))) {
    warning(.ts_closeness_msg, call. = FALSE)
    dots <- dots[setdiff(names(dots), c("closeness", "perturb", "dependence"))]
  }
  .warn_unused_dots(dots, "generate")
  n <- .check_count(n, allow_zero = FALSE)
  Y <- .with_seed(seed, .generate_ts(object, n))
  colnames(Y) <- object$names
  if (object$univariate) Y <- Y[, 1L]
  stats::ts(Y, start = start, frequency = object$frequency)
}

.generate_ts <- function(object, n) {
  X <- object$data
  N <- object$n_train
  switch(object$method,
    iid = as.matrix(generate(object$cross, n = n)),
    block = {
      L <- object$block_length
      nb <- ceiling(n / L)
      starts <- sample.int(N, nb, replace = TRUE)
      idx <- as.vector(vapply(starts, function(s) ((s - 1L + seq_len(L) - 1L) %% N) + 1L,
                              integer(L)))
      X[idx[seq_len(n)], , drop = FALSE]
    },
    stationary = {
      L <- object$block_length
      idx <- integer(n)
      idx[1L] <- sample.int(N, 1L)
      if (n > 1L) {
        jump <- stats::runif(n - 1L) < 1 / L
        new_start <- sample.int(N, n - 1L, replace = TRUE)
        for (t in 2:n) {
          idx[t] <- if (jump[t - 1L]) new_start[t - 1L] else (idx[t - 1L] %% N) + 1L
        }
      }
      X[idx, , drop = FALSE]
    },
    copula_var = {
      Z <- .simulate_var(object$var, n, object$burn_in)
      U <- stats::pnorm(sweep(Z, 2L, object$z_sd, "/"))
      vapply(seq_along(object$marginals), function(j) {
        .marginal_quantile(object$marginals[[j]], U[, j])
      }, numeric(n))
    }
  ) -> Y
  matrix(Y, nrow = n)
}

#' @export
print.sp_ts_synthesizer <- function(x, ...) {
  cat("<sp_ts_synthesizer>\n")
  cat(sprintf("  trained on : %d time points, %d series (frequency %g)\n",
              x$n_train, length(x$names), x$frequency))
  cat("  method     :", x$method)
  if (x$method == "copula_var") cat(sprintf(" (VAR order %d)", x$var$order))
  if (x$method %in% c("block", "stationary")) {
    cat(sprintf(" (block length %d)", x$block_length))
  }
  cat("\n")
  invisible(x)
}
