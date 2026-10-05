# Multivariate distributions -------------------------------------------------

#' Build a correlation matrix with a given structure
#'
#' @param d Dimension.
#' @param type Structure: `"exchangeable"` (all off-diagonal entries `rho`),
#'   `"ar1"` (`rho^|i-j|`), `"toeplitz"` (first row given by `rho`, a vector
#'   of length `d - 1`), `"random"` (a random positive-definite correlation
#'   matrix) or `"identity"`.
#' @param rho Correlation parameter (see `type`).
#' @param names Optional variable names.
#' @return A `d x d` correlation matrix.
#' @examples
#' make_corr(4, "ar1", rho = 0.7)
#' make_corr(3, "toeplitz", rho = c(0.5, 0.2))
#' @export
make_corr <- function(d, type = c("exchangeable", "ar1", "toeplitz", "random", "identity"),
                      rho = 0.5, names = NULL) {
  d <- .check_count(d, "d", allow_zero = FALSE)
  type <- match.arg(type)
  R <- switch(type,
    identity = diag(d),
    exchangeable = {
      if (length(rho) != 1L || rho <= -1 / (d - 1) || rho >= 1) {
        stop("For 'exchangeable', 'rho' must lie in (-1/(d-1), 1).", call. = FALSE)
      }
      M <- matrix(rho, d, d)
      diag(M) <- 1
      M
    },
    ar1 = {
      if (length(rho) != 1L || abs(rho) >= 1) {
        stop("For 'ar1', 'rho' must lie in (-1, 1).", call. = FALSE)
      }
      rho^abs(outer(seq_len(d), seq_len(d), "-"))
    },
    toeplitz = {
      if (length(rho) != d - 1L) stop("For 'toeplitz', 'rho' must have length d - 1.", call. = FALSE)
      stats::toeplitz(c(1, rho))
    },
    random = {
      A <- matrix(stats::rnorm(d * (d + 2L)), d + 2L, d)
      stats::cov2cor(crossprod(A) + diag(stats::runif(d, 0.1, 1), d))
    }
  )
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("The resulting matrix is not positive definite.", call. = FALSE)
  }
  if (!is.null(names)) dimnames(R) <- list(names, names)
  R
}

.mat_sqrt <- function(sigma) {
  sigma <- as.matrix(sigma)
  if (!isSymmetric(unname(sigma), tol = 1e-8)) stop("'sigma' must be symmetric.", call. = FALSE)
  ch <- try(chol(sigma), silent = TRUE)
  if (!inherits(ch, "try-error")) return(ch)
  e <- eigen(sigma, symmetric = TRUE)
  if (min(e$values) < -1e-8 * max(abs(e$values))) {
    stop("'sigma' must be positive semi-definite.", call. = FALSE)
  }
  t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values)))
}

# A correlation matrix: square, finite, symmetric, unit diagonal, entries in
# [-1, 1] and positive semi-definite
.validate_corr <- function(R, arg = "corr", tol = 1e-8) {
  R <- as.matrix(R)
  if (!is.numeric(R) || nrow(R) != ncol(R) || nrow(R) < 1L) {
    stop(sprintf("'%s' must be a square numeric matrix.", arg), call. = FALSE)
  }
  if (any(!is.finite(R))) stop(sprintf("'%s' must contain finite values.", arg), call. = FALSE)
  if (!isSymmetric(unname(R), tol = tol)) stop(sprintf("'%s' must be symmetric.", arg), call. = FALSE)
  if (any(abs(diag(R) - 1) > tol)) {
    stop(sprintf("'%s' must be a correlation matrix with unit diagonal; use stats::cov2cor() to convert a covariance matrix.",
                 arg), call. = FALSE)
  }
  if (any(abs(R) > 1 + tol)) stop(sprintf("Entries of '%s' must lie in [-1, 1].", arg), call. = FALSE)
  .mat_sqrt(R)   # fails if not positive semi-definite
  R
}

# Degrees of freedom: positive number; Inf gives the Gaussian case
.check_df <- function(df) {
  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0) {
    stop("'df' must be a single positive number (Inf gives the Gaussian case).", call. = FALSE)
  }
  df
}

.mv_names <- function(M, d, prefix = "X") {
  nm <- colnames(M)
  if (is.null(nm)) nm <- paste0(prefix, seq_len(d))
  nm
}

#' Multivariate random number generators
#'
#' @description
#' Random generation for common multivariate distributions. All functions
#' return an `n x d` matrix with column names, and accept an optional `seed`.
#'
#' * `r_mvnorm()`: multivariate normal \eqn{N(\mu, \Sigma)}.
#' * `r_mvt()`: multivariate Student-*t* with `df` degrees of freedom and
#'   scale matrix `sigma`.
#' * `r_mvlnorm()`: multivariate log-normal: `exp()` of a multivariate normal
#'   with parameters `meanlog`, `sigmalog`.
#' * `r_mvskewnorm()`: multivariate skew-normal of Azzalini and Dalla Valle
#'   (1996) with location `xi`, scale matrix `omega` and shape `alpha`.
#' * `r_dirichlet()`: Dirichlet distribution on the simplex.
#' * `r_mvmixture()`: finite mixture of multivariate normals.
#'
#' @param n Number of draws.
#' @param mean,xi Mean / location vector.
#' @param sigma,omega Covariance (scale) matrix.
#' @param df Degrees of freedom (`Inf` gives the normal distribution).
#' @param meanlog,sigmalog Mean vector and covariance matrix on the log scale.
#' @param alpha Shape (skewness) vector for `r_mvskewnorm()`; concentration
#'   vector for `r_dirichlet()`.
#' @param weights Mixture weights (normalised internally).
#' @param means List of component mean vectors.
#' @param sigmas List of component covariance matrices.
#' @param return_component Should `r_mvmixture()` add the component label as
#'   attribute `"component"`?
#' @param seed Optional random seed.
#' @return A numeric matrix with `n` rows.
#' @references Azzalini, A. and Dalla Valle, A. (1996). The multivariate
#'   skew-normal distribution. *Biometrika*, 83(4), 715--726.
#'   \doi{10.1093/biomet/83.4.715}
#' @examples
#' S <- make_corr(3, "ar1", 0.6)
#' x <- r_mvnorm(1000, mean = c(0, 1, 2), sigma = S, seed = 1)
#' round(cor(x), 2)
#'
#' t3 <- r_mvt(500, df = 3, sigma = S)
#' sn <- r_mvskewnorm(500, omega = S, alpha = c(5, 0, -5))
#' dir <- r_dirichlet(5, alpha = c(1, 2, 3))
#' rowSums(dir)
#'
#' mix <- r_mvmixture(500, weights = c(0.3, 0.7),
#'                    means = list(c(-2, -2), c(2, 2)),
#'                    sigmas = list(diag(2), diag(2)))
#' @name mvdist
NULL

#' @rdname mvdist
#' @export
r_mvnorm <- function(n, mean = NULL, sigma, seed = NULL) {
  n <- .check_count(n)
  sigma <- as.matrix(sigma)
  d <- ncol(sigma)
  mean <- mean %||% rep(0, d)
  if (length(mean) != d) stop("'mean' and 'sigma' have incompatible dimensions.", call. = FALSE)
  C <- .mat_sqrt(sigma)
  out <- .with_seed(seed, matrix(stats::rnorm(n * d), n, d) %*% C)
  out <- sweep(out, 2L, mean, "+")
  colnames(out) <- .mv_names(sigma, d)
  out
}

#' @rdname mvdist
#' @export
r_mvt <- function(n, df, mean = NULL, sigma, seed = NULL) {
  .check_df(df)
  .with_seed(seed, {
    z <- r_mvnorm(n, mean = NULL, sigma = sigma)
    if (is.finite(df)) z <- z / sqrt(stats::rchisq(nrow(z), df) / df)
    d <- ncol(z)
    mean <- mean %||% rep(0, d)
    if (length(mean) != d) stop("'mean' and 'sigma' have incompatible dimensions.", call. = FALSE)
    sweep(z, 2L, mean, "+")
  })
}

#' @rdname mvdist
#' @export
r_mvlnorm <- function(n, meanlog = NULL, sigmalog, seed = NULL) {
  exp(r_mvnorm(n, mean = meanlog, sigma = sigmalog, seed = seed))
}

#' @rdname mvdist
#' @export
r_mvskewnorm <- function(n, xi = NULL, omega, alpha, seed = NULL) {
  n <- .check_count(n)
  omega <- as.matrix(omega)
  d <- ncol(omega)
  xi <- xi %||% rep(0, d)
  if (length(alpha) != d || length(xi) != d) {
    stop("'xi', 'omega' and 'alpha' have incompatible dimensions.", call. = FALSE)
  }
  w <- sqrt(diag(omega))
  Obar <- omega / outer(w, w)
  delta <- as.vector(Obar %*% alpha) / sqrt(1 + as.numeric(t(alpha) %*% Obar %*% alpha))
  big <- rbind(c(1, delta), cbind(delta, Obar))
  out <- .with_seed(seed, {
    z <- matrix(stats::rnorm(n * (d + 1L)), n, d + 1L) %*% .mat_sqrt(big)
    sgn <- ifelse(z[, 1L] > 0, 1, -1)
    z[, -1L, drop = FALSE] * sgn
  })
  out <- sweep(sweep(out, 2L, w, "*"), 2L, xi, "+")
  colnames(out) <- .mv_names(omega, d)
  out
}

#' @rdname mvdist
#' @export
r_dirichlet <- function(n, alpha, seed = NULL) {
  n <- .check_count(n)
  if (!is.numeric(alpha) || any(alpha <= 0) || length(alpha) < 2L) {
    stop("'alpha' must be a vector of at least 2 positive numbers.", call. = FALSE)
  }
  d <- length(alpha)
  g <- .with_seed(seed, matrix(stats::rgamma(n * d, shape = rep(alpha, each = n)), n, d))
  out <- g / rowSums(g)
  colnames(out) <- names(alpha) %||% paste0("X", seq_len(d))
  out
}

#' @rdname mvdist
#' @export
r_mvmixture <- function(n, weights, means, sigmas, return_component = FALSE,
                        seed = NULL) {
  n <- .check_count(n)
  k <- length(weights)
  if (!is.list(means) || !is.list(sigmas)) {
    stop("'means' and 'sigmas' must be lists (one element per component).", call. = FALSE)
  }
  if (k < 1L || length(means) != k || length(sigmas) != k) {
    stop("'weights', 'means' and 'sigmas' must have the same length.", call. = FALSE)
  }
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("'weights' must be finite, non-negative and not all zero.", call. = FALSE)
  }
  d <- length(means[[1L]])
  for (j in seq_len(k)) {
    sj <- as.matrix(sigmas[[j]])
    if (length(means[[j]]) != d || !identical(dim(sj), c(d, d))) {
      stop(sprintf("Component %d: all means must have length %d and all sigmas must be %d x %d.",
                   j, d, d, d), call. = FALSE)
    }
  }
  .with_seed(seed, {
    comp <- sample.int(k, n, replace = TRUE, prob = weights / sum(weights))
    out <- matrix(NA_real_, n, d)
    for (j in seq_len(k)) {
      idx <- which(comp == j)
      if (length(idx)) out[idx, ] <- r_mvnorm(length(idx), means[[j]], sigmas[[j]])
    }
    colnames(out) <- .mv_names(as.matrix(sigmas[[1L]]), d)
    if (return_component) attr(out, "component") <- comp
    out
  })
}

#' Sample from a copula
#'
#' @description
#' Draws uniforms on \eqn{[0,1]^d} from elliptical (Gaussian, *t*) or
#' Archimedean (Clayton, Gumbel, Frank) copulas. Archimedean copulas in
#' dimension `d > 2` are sampled with the Marshall--Olkin frailty algorithm;
#' the bivariate Frank copula also allows negative dependence (`theta < 0`),
#' sampled by conditional inversion.
#'
#' @param n Number of draws.
#' @param family Copula family.
#' @param dim Dimension (ignored when `corr` is given).
#' @param corr Correlation matrix for elliptical copulas (default:
#'   exchangeable with `rho = 0.5`). It must have a unit diagonal; a
#'   covariance matrix is rejected (convert it with [stats::cov2cor()]).
#' @param df Degrees of freedom for the *t* copula: a positive number;
#'   `Inf` gives the Gaussian copula.
#' @param theta Archimedean parameter: Clayton `theta > 0`, Gumbel
#'   `theta >= 1`, Frank `theta != 0` (`> 0` when `dim > 2`).
#' @param seed Optional random seed.
#' @return An `n x dim` matrix of uniforms.
#' @references Marshall, A. W. and Olkin, I. (1988). Families of multivariate
#'   distributions. *Journal of the American Statistical Association*,
#'   83(403), 834--841. \doi{10.1080/01621459.1988.10478671}
#'
#'   Hofert, M. (2008). Sampling Archimedean copulas. *Computational
#'   Statistics & Data Analysis*, 52(12), 5163--5174.
#'   \doi{10.1016/j.csda.2008.05.019}
#' @examples
#' u <- r_copula(1000, "clayton", dim = 3, theta = 2, seed = 1)
#' cor(u, method = "kendall")   # theoretical tau = theta / (theta + 2) = 0.5
#' @export
r_copula <- function(n, family = c("gaussian", "t", "clayton", "gumbel", "frank", "independence"),
                     dim = 2L, corr = NULL, df = 4, theta = NULL, seed = NULL) {
  n <- .check_count(n)
  family <- match.arg(family)
  if (!is.null(corr)) dim <- ncol(corr)
  d <- .check_count(dim, "dim", allow_zero = FALSE)
  .with_seed(seed, .r_copula_impl(n, family, d, corr, df, theta))
}

.r_copula_impl <- function(n, family, d, corr, df, theta) {
  if (family %in% c("clayton", "gumbel", "frank") &&
      (is.null(theta) || !is.numeric(theta) || length(theta) != 1L || !is.finite(theta))) {
    stop("'theta' must be a single finite number for Archimedean copulas.", call. = FALSE)
  }
  U <- switch(family,
    independence = matrix(stats::runif(n * d), n, d),
    gaussian = , t = {
      R <- if (is.null(corr)) make_corr(d, "exchangeable", 0.5) else .validate_corr(corr)
      if (family == "t") .check_df(df)
      cop <- list(family = family, R = R, chol = .mat_sqrt(R),
                  df = if (family == "t") df else Inf)
      .rcopula_fitted(n, cop)
    },
    clayton = {
      if (theta <= 0) stop("Clayton requires theta > 0.", call. = FALSE)
      V <- stats::rgamma(n, shape = 1 / theta)
      E <- matrix(stats::rexp(n * d), n, d)
      (1 + E / V)^(-1 / theta)
    },
    gumbel = {
      if (theta < 1) stop("Gumbel requires theta >= 1.", call. = FALSE)
      if (theta == 1) return(matrix(stats::runif(n * d), n, d))
      a <- 1 / theta
      V <- .r_pos_stable(n, a)
      E <- matrix(stats::rexp(n * d), n, d)
      exp(-(E / V)^a)
    },
    frank = {
      if (theta == 0) return(matrix(stats::runif(n * d), n, d))
      if (d == 2L) {
        u <- stats::runif(n)
        w <- stats::runif(n)
        v <- -log1p(w * expm1(-theta) / (w + (1 - w) * exp(-theta * u))) / theta
        cbind(u, v)
      } else {
        if (theta < 0) stop("Frank with dim > 2 requires theta > 0.", call. = FALSE)
        V <- .r_logseries(n, -expm1(-theta))
        E <- matrix(stats::rexp(n * d), n, d)
        -log1p(exp(-E / V) * expm1(-theta)) / theta
      }
    }
  )
  U <- matrix(U, n, d)
  colnames(U) <- colnames(corr) %||% paste0("U", seq_len(d))
  U
}

# Positive stable variate with Laplace transform exp(-t^a), 0 < a < 1
# (Kanter, 1975; Chambers, Mallows and Stuck, 1976)
.r_pos_stable <- function(n, a) {
  th <- stats::runif(n, 0, pi)
  w <- stats::rexp(n)
  (sin(a * th) / sin(th)^(1 / a)) * (sin((1 - a) * th) / w)^((1 - a) / a)
}

# Logarithmic series variate, P(V = k) = -p^k / (k log(1 - p)) (Kemp, 1981)
.r_logseries <- function(n, p) {
  h <- log1p(-p)
  u2 <- stats::runif(n)
  u1 <- stats::runif(n)
  q <- -expm1(u1 * h)
  out <- rep(1, n)
  small <- u2 < q * q
  out[small] <- floor(1 + log(u2[small]) / log(q[small]))
  mid <- !small & u2 <= q
  out[mid] <- 2
  out[u2 > p] <- 1
  out
}

#' Generate data from a copula with arbitrary marginal distributions
#'
#' @description
#' `r_mvdist()` combines a copula (see [r_copula()]) with user-specified
#' marginal distributions, returning a data frame. This is a flexible way to
#' simulate mixed-type data sets with known dependence, e.g. for simulation
#' studies or to benchmark synthesizers.
#'
#' Each element of `margins` is either
#' * a quantile function `function(p)` (e.g. from [margin_categorical()] or
#'   [margin_empirical()]),
#' * a character string naming a distribution with a `q*` function, e.g.
#'   `"norm"`, `"gamma"`, `"pois"`, or
#' * a list with element `dist` and further parameters, e.g.
#'   `list(dist = "gamma", shape = 2, rate = 1)`.
#'
#' @param n Number of rows.
#' @param margins A (named) list of marginal specifications.
#' @param copula Copula family, see [r_copula()].
#' @param ... Further arguments passed to [r_copula()] (`corr`, `df`,
#'   `theta`).
#' @param seed Optional random seed.
#' @return A data frame with one column per margin.
#' @examples
#' df <- r_mvdist(
#'   500,
#'   margins = list(
#'     income = list(dist = "lnorm", meanlog = 10, sdlog = 0.5),
#'     age = list(dist = "unif", min = 18, max = 80),
#'     sector = margin_categorical(c("agri", "industry", "services"),
#'                                 c(0.1, 0.3, 0.6))
#'   ),
#'   copula = "gaussian", corr = make_corr(3, "exchangeable", 0.4), seed = 1
#' )
#' str(df)
#' @export
r_mvdist <- function(n, margins, copula = "gaussian", ..., seed = NULL) {
  if (!is.list(margins) || !length(margins)) {
    stop("'margins' must be a non-empty list.", call. = FALSE)
  }
  d <- length(margins)
  nm <- names(margins)
  if (is.null(nm) || any(nm == "")) nm <- paste0("X", seq_len(d))
  .with_seed(seed, {
    U <- r_copula(n, family = copula, dim = d, ...)
    if (ncol(U) != d) stop("The copula dimension does not match the number of margins.", call. = FALSE)
    cols <- lapply(seq_len(d), function(j) .apply_margin(margins[[j]], U[, j]))
    names(cols) <- nm
    out <- cols
    attr(out, "row.names") <- .set_row_names(n)
    class(out) <- "data.frame"
    out
  })
}

.apply_margin <- function(spec, u) {
  if (is.function(spec)) return(spec(u))
  if (is.character(spec) && length(spec) == 1L) spec <- list(dist = spec)
  if (is.list(spec) && is.character(spec$dist)) {
    qf <- get0(paste0("q", spec$dist), mode = "function")
    if (is.null(qf)) stop("No quantile function 'q", spec$dist, "' found.", call. = FALSE)
    args <- spec[setdiff(names(spec), "dist")]
    return(do.call(qf, c(list(u), args)))
  }
  stop("Invalid margin specification.", call. = FALSE)
}

#' Marginal quantile functions for r_mvdist()
#'
#' `margin_categorical()` builds a quantile function for a categorical
#' variable; `margin_empirical()` builds the (interpolated) empirical
#' quantile function of observed data, preserving its class.
#'
#' @param levels Category labels.
#' @param probs Category probabilities (normalised internally).
#' @param ordered Return an ordered factor?
#' @param x Observed data (numeric, integer, Date, POSIXct, factor,
#'   character, logical).
#' @param interpolation `"linear"` or `"spline"`.
#' @return A function of `p` returning values.
#' @examples
#' q <- margin_categorical(c("low", "mid", "high"), c(0.2, 0.5, 0.3), ordered = TRUE)
#' table(q(runif(1000)))
#' qe <- margin_empirical(faithful$eruptions)
#' summary(qe(runif(1000)))
#' @export
margin_categorical <- function(levels, probs = NULL, ordered = FALSE) {
  levels <- as.character(levels)
  if (!length(levels) || anyNA(levels) || anyDuplicated(levels)) {
    stop("'levels' must be non-missing and unique.", call. = FALSE)
  }
  probs <- probs %||% rep(1, length(levels))
  if (!is.numeric(probs) || length(probs) != length(levels) || any(!is.finite(probs)) ||
      any(probs < 0) || sum(probs) <= 0) {
    stop("'probs' must be finite, non-negative and match 'levels' in length.", call. = FALSE)
  }
  m <- list(type = "categorical", values = levels, probs = probs / sum(probs))
  function(p) {
    .check_unit_interval(p)
    factor(.marginal_quantile(m, p), levels = levels, ordered = ordered)
  }
}

.check_unit_interval <- function(p) {
  if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1)) {
    stop("Probabilities must be numbers in [0, 1].", call. = FALSE)
  }
  invisible(p)
}

#' @rdname margin_categorical
#' @export
margin_empirical <- function(x, interpolation = c("linear", "spline")) {
  interpolation <- match.arg(interpolation)
  m <- .fit_marginal(x, "x", interpolation = interpolation, discrete_threshold = 0L)
  if (m$type == "empty") stop("'x' has no observed values.", call. = FALSE)
  proto <- x[0L]
  function(p) {
    .check_unit_interval(p)
    .restore_class(.marginal_quantile(m, p), m, proto)
  }
}
