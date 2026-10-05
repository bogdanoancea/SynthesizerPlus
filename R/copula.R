# Copula models --------------------------------------------------------------

# Fit an elliptical copula to a matrix of pseudo-observations `U` (may contain
# NA; correlations are estimated from pairwise complete observations).
.fit_copula <- function(U, family = c("gaussian", "t", "independence"),
                        pd_method = c("auto", "higham", "eigen", "none"),
                        df_grid = c(2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100),
                        max_rows = 5000L) {
  family <- match.arg(family)
  pd_method <- match.arg(pd_method)
  d <- ncol(U)
  nm <- colnames(U)
  if (family == "independence" || d < 2L) {
    R <- diag(d)
    dimnames(R) <- list(nm, nm)
    return(list(family = "independence", R = R, df = Inf, chol = R))
  }

  if (family == "gaussian") {
    Z <- stats::qnorm(U)
    R <- suppressWarnings(stats::cor(Z, use = "pairwise.complete.obs"))
  } else {
    # Kendall's tau inversion is robust for the t copula: rho = sin(pi/2 tau)
    Us <- if (nrow(U) > 2000L) U[sample.int(nrow(U), 2000L), , drop = FALSE] else U
    tau <- suppressWarnings(stats::cor(Us, method = "kendall",
                                       use = "pairwise.complete.obs"))
    R <- sin(pi / 2 * tau)
  }
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R <- (R + t(R)) / 2
  R <- .make_pd(R, method = pd_method)
  dimnames(R) <- list(nm, nm)

  df <- Inf
  if (family == "t") {
    df <- .estimate_t_df(U, R, df_grid = df_grid, max_rows = max_rows)
  }
  ch <- tryCatch(chol(R), error = function(e) {
    stop("The copula correlation matrix is not positive definite (variables are ",
         "perfectly dependent). Use pd_method = \"auto\" to repair it.", call. = FALSE)
  })
  list(family = family, R = R, df = df, chol = ch)
}

# Profile likelihood for the degrees of freedom of a t copula on a grid
.estimate_t_df <- function(U, R, df_grid, max_rows = 5000L) {
  cc <- stats::complete.cases(U)
  Uc <- U[cc, , drop = FALSE]
  if (nrow(Uc) < 10L) return(10)
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
  df_grid[which.max(ll)]
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
