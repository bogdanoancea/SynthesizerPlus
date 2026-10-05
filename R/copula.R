# Copula models --------------------------------------------------------------

# Fit an elliptical copula to a matrix of pseudo-observations `U` (may contain
# NA; correlations are estimated from pairwise complete observations).
# Pairs observed together in fewer than `min_pair_obs` rows (or whose
# correlation is not finite) are *unidentified*: their correlation is filled
# by the maximum-determinant completion of the identified entries (falling
# back to 0 if that fails). The counts are kept in `n_pair`; they always refer
# to the full data, never to a computational subsample.
.fit_copula <- function(U, family = c("gaussian", "t", "independence"),
                        pd_method = c("auto", "higham", "eigen", "none"),
                        df_grid = c(2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100),
                        max_rows = 5000L, min_pair_obs = 3L, tau_rows = 2000L) {
  family <- match.arg(family)
  pd_method <- match.arg(pd_method)
  d <- ncol(U)
  nm <- colnames(U)
  if (family == "independence" || d < 2L) {
    R <- diag(d)
    dimnames(R) <- list(nm, nm)
    obs <- !is.na(U)
    np <- crossprod(obs + 0)
    storage.mode(np) <- "integer"
    dimnames(np) <- list(nm, nm)
    return(list(family = "independence", R = R, df = Inf, df_estimated = NA,
                chol = R, n_pair = np, unidentified = 0L, completion = "none",
                pd_adjustment = 0))
  }

  # identification is decided on the full data: number of rows in which
  # both variables of each pair are observed
  obs <- !is.na(U)
  np <- crossprod(obs + 0)
  storage.mode(np) <- "integer"
  dimnames(np) <- list(nm, nm)

  if (family == "gaussian") {
    Z <- stats::qnorm(U)
    R <- suppressWarnings(stats::cor(Z, use = "pairwise.complete.obs"))
  } else {
    R <- .pairwise_kendall_corr(U, obs, np, tau_rows = tau_rows)
  }
  unid <- (np < min_pair_obs) | !is.finite(R)
  diag(unid) <- FALSE
  R[unid] <- NA
  diag(R) <- 1
  R <- (R + t(R)) / 2
  completion <- "none"
  if (any(unid)) {
    cmp <- .complete_corr(R, unid)
    R <- cmp$R
    completion <- cmp$method
  }
  R_before <- R
  R <- .make_pd(R, method = pd_method)
  pd_adjustment <- max(abs(R - R_before))
  dimnames(R) <- list(nm, nm)

  df <- Inf
  df_estimated <- NA
  if (family == "t") {
    est <- .estimate_t_df(U, R, df_grid = df_grid, max_rows = max_rows)
    df <- est$df
    df_estimated <- est$estimated
  }
  ch <- tryCatch(chol(R), error = function(e) {
    stop("The copula correlation matrix is not positive definite (variables are ",
         "perfectly dependent). Use pd_method = \"auto\" to repair it.", call. = FALSE)
  })
  list(family = family, R = R, df = df, df_estimated = df_estimated, chol = ch,
       n_pair = np, unidentified = as.integer(sum(unid) / 2),
       completion = completion, pd_adjustment = pd_adjustment)
}

# Latent correlations of a t copula by inverting Kendall's tau,
# rho = sin(pi/2 tau). Kendall's tau is O(n^2), so one common subsample of
# `tau_rows` rows serves all pairs. A pair that the common subsample leaves
# with fewer than min(n_pair, tau_min_rows) jointly observed rows (e.g. a
# pair with sparse overlap) is re-estimated from up to `tau_rows` of its own
# joint rows. Hence every identified pair (n_pair >= 3) is estimated, from at
# least min(n_pair, tau_min_rows) rows, whatever the random subsample.
.pairwise_kendall_corr <- function(U, obs, np, tau_rows = 2000L, tau_min_rows = 500L) {
  n <- nrow(U)
  rows <- if (n > tau_rows) sort(sample.int(n, tau_rows)) else seq_len(n)
  tau <- .kendall_matrix(U[rows, , drop = FALSE])
  if (n > tau_rows) {
    np_sub <- crossprod(obs[rows, , drop = FALSE] + 0)
    redo <- np_sub < pmin(np, tau_min_rows) & np > 0
    redo[lower.tri(redo, diag = TRUE)] <- FALSE
    ij <- which(redo, arr.ind = TRUE)
    for (k in seq_len(nrow(ij))) {
      i <- ij[k, 1L]
      j <- ij[k, 2L]
      jr <- which(obs[, i] & obs[, j])
      if (length(jr) > tau_rows) jr <- jr[sample.int(length(jr), tau_rows)]
      tau[i, j] <- tau[j, i] <- .kendall_matrix(U[jr, c(i, j), drop = FALSE])[1L, 2L]
    }
  }
  sin(pi / 2 * tau)
}

# Pairwise-complete Kendall's tau of every pair of columns of X (which may
# contain NA), for data without ties (pseudo-observations are continuous:
# ranks with random tie-breaking or the randomised distributional
# transform), where tau-a = tau-b:
#   tau_ij = sum_{k != l} sgn(x_ki - x_li) sgn(x_kj - x_lj) / (n_ij (n_ij - 1)),
# the sum running over rows where both columns are observed (a missing value
# contributes a zero sign). The sum for all pairs is the cross-product of the
# vectorised sign matrices, accumulated over blocks of rows so that memory
# stays bounded; this replaces p^2 / 2 separate O(n^2) calls of
# stats::cor(method = "kendall") with BLAS matrix products.
.kendall_matrix <- function(X, block_cells = 4e6) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  obs <- !is.na(X)
  X[!obs] <- 0
  num <- matrix(0, p, p)
  bs <- max(1L, as.integer(block_cells %/% max(1, n * p)))
  for (start in seq.int(1L, n, by = bs)) {
    B <- start:min(n, start + bs - 1L)
    M <- vapply(seq_len(p), function(j) {
      as.vector(sign(outer(X[B, j], X[, j], "-")) * outer(obs[B, j], obs[, j]))
    }, numeric(length(B) * n))
    if (!is.matrix(M)) M <- matrix(M, ncol = p)
    num <- num + crossprod(M)
  }
  nij <- crossprod(obs + 0)
  tau <- num / (nij * (nij - 1))
  tau[nij < 2] <- NA
  diag(tau) <- 1
  dimnames(tau) <- list(colnames(X), colnames(X))
  tau
}
# Maximum-determinant (maximum-entropy) positive-definite completion of a
# correlation matrix whose entries flagged in `unid` are unknown: the
# identified entries are kept and the unknown ones are chosen so that the
# inverse has zeros there, i.e. the corresponding pairs are conditionally
# independent given the remaining variables (Dempster 1972, covariance
# selection; algorithm 17.1 of Hastie, Tibshirani and Friedman 2009).
# Falls back to 0 for the unknown entries if the identified entries are not
# completable (inconsistent pairwise estimates) or the iteration fails.
# Returns list(R, method) with method "maxdet" or "zero_fallback".
.complete_corr <- function(R, unid, tol = 1e-10, max_iter = 500L) {
  zero_fill <- function() {
    R0 <- R
    R0[unid] <- 0
    list(R = R0, method = "zero_fallback")
  }
  p <- ncol(R)
  S <- R
  S[unid] <- 0
  W <- S
  ok <- tryCatch({
    for (it in seq_len(max_iter)) {
      W_old <- W
      for (j in seq_len(p)) {
        oth <- seq_len(p)[-j]
        nb <- oth[!unid[j, oth]]
        w12 <- if (length(nb)) {
          W[oth, nb, drop = FALSE] %*% solve(W[nb, nb, drop = FALSE], S[nb, j])
        } else {
          rep(0, p - 1L)
        }
        W[oth, j] <- w12
        W[j, oth] <- w12
      }
      if (max(abs(W - W_old)) < tol) break
    }
    all(is.finite(W)) &&
      min(eigen(W, symmetric = TRUE, only.values = TRUE)$values) > 0 &&
      max(abs(W[!unid] - S[!unid])) < 1e-6
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(zero_fill())
  W <- (W + t(W)) / 2
  diag(W) <- 1
  list(R = W, method = "maxdet")
}

# Profile likelihood for the degrees of freedom of a t copula on a grid (R held
# fixed at its Kendall's-tau estimate). Returns list(df, estimated): with fewer
# than `min_rows` complete rows nothing is estimated and df is fixed at
# `fallback` (estimated = FALSE).
.estimate_t_df <- function(U, R, df_grid, max_rows = 5000L, min_rows = 10L,
                           fallback = 10) {
  cc <- stats::complete.cases(U)
  Uc <- U[cc, , drop = FALSE]
  if (nrow(Uc) < min_rows) return(list(df = fallback, estimated = FALSE))
  if (nrow(Uc) > max_rows) Uc <- Uc[sample.int(nrow(Uc), max_rows), , drop = FALSE]
  Uc <- pmin(pmax(Uc, 1e-10), 1 - 1e-10)
  d <- ncol(Uc)
  Rinv <- solve(R)
  logdet <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
  ll <- vapply(df_grid, function(nu) {
    X <- stats::qt(Uc, df = nu)
    q <- rowSums((X %*% Rinv) * X)
    lmv <- lgamma((nu + d) / 2) - lgamma(nu / 2) - d / 2 * log(nu * pi) -
      0.5 * logdet - (nu + d) / 2 * log1p(q / nu)
    lmarg <- rowSums(stats::dt(X, df = nu, log = TRUE))
    sum(lmv - lmarg)
  }, numeric(1))
  list(df = df_grid[which.max(ll)], estimated = TRUE)
}

# Draw n x d uniforms from a fitted copula. `dependence` in [0, 2] scales the
# off-diagonal correlations: R_s = (1 - lambda) I + lambda R (repaired to a
# valid correlation matrix when lambda > 1).
# This keeps every margin exactly uniform (unlike rescaling the Cholesky
# factor, which changes the variance of the latent normals).
.rcopula_fitted <- function(n, cop, dependence = 1) {
  d <- ncol(cop$R)
  if (n == 0L) return(matrix(numeric(0), 0L, d))
  if (cop$family == "independence" || dependence == 0 || d < 2L) {
    return(matrix(stats::runif(n * d), n, d))
  }
  C <- cop$chol
  if (dependence != 1) {
    Rs <- (1 - dependence) * diag(d) + dependence * cop$R
    if (dependence > 1) {
      # strengthened dependence: keep a valid correlation matrix
      off <- row(Rs) != col(Rs)
      Rs[off] <- pmin(pmax(Rs[off], -0.999), 0.999)
      Rs <- .make_pd(Rs, method = "higham", eps = 1e-4)
    }
    C <- chol(Rs)
  }
  Z <- matrix(stats::rnorm(n * d), n, d) %*% C
  if (is.finite(cop$df)) {
    w <- sqrt(stats::rchisq(n, df = cop$df) / cop$df)
    stats::pt(Z / w, df = cop$df)
  } else {
    stats::pnorm(Z)
  }
}

# Ensure a correlation matrix is positive definite
.make_pd <- function(R, method = c("auto", "higham", "eigen", "none"),
                     eps = 1e-6) {
  method <- match.arg(method)
  if (method == "none") return(R)
  is_pd <- !inherits(try(chol(R), silent = TRUE), "try-error") &&
    min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > eps
  if (is_pd) return(R)
  if (method %in% c("auto", "higham")) R <- .nearest_cor(R)
  .eigen_floor(R, eps = eps)
}

# Higham (2002) alternating projections: nearest correlation matrix
.nearest_cor <- function(A, tol = 1e-9, maxit = 500L) {
  Y <- A
  dS <- matrix(0, nrow(A), ncol(A))
  for (k in seq_len(maxit)) {
    Rk <- Y - dS
    e <- eigen(Rk, symmetric = TRUE)
    X <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
    dS <- X - Rk
    Ynew <- X
    diag(Ynew) <- 1
    if (max(abs(Ynew - Y)) < tol) {
      Y <- Ynew
      break
    }
    Y <- Ynew
  }
  (Y + t(Y)) / 2
}

.eigen_floor <- function(R, eps = 1e-6) {
  e <- eigen(R, symmetric = TRUE)
  lam <- pmax(e$values, eps)
  M <- e$vectors %*% (lam * t(e$vectors))
  s <- sqrt(diag(M))
  M <- M / outer(s, s)
  diag(M) <- 1
  (M + t(M)) / 2
}
