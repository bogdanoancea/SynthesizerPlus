# Latent (Gaussian-copula) correlations for discrete margins ------------------
#
# For a discrete or categorical margin the copula coordinate is only known up
# to the cell (F(x-), F(x)]. The pseudo-observations fill the cell with
# independent uniform jitter (the distributional transform); correlating the
# jittered normal scores therefore *attenuates* the latent correlation, most
# strongly for binary variables such as missingness indicators. Pairs that
# involve a discrete margin are instead estimated by maximum likelihood under
# the Gaussian copula: polyserial (discrete-continuous) and polychoric
# (discrete-discrete) correlations, with thresholds fixed at the normal
# quantiles of the marginal cumulative probabilities (two-step estimator;
# Olsson 1979; Olsson, Drasgow and Dorans 1982).

# 20-point Gauss-Legendre nodes and weights on [-1, 1]
.gl20 <- local({
  x <- c(-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
         -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
         -0.5108670019508271, -0.3737060887154195, -0.2277858511416451,
         -0.0765265211334973)
  w <- c(0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
         0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
         0.1316886384491766, 0.1420961093183820, 0.1491729864726037,
         0.1527533871307258)
  list(x = c(x, -rev(x)), w = c(w, rev(w)))
})

# Bivariate standard normal CDF P(Z1 <= a, Z2 <= b) with correlation rho, for
# vectors a, b (infinite values allowed), via Plackett's identity in the
# arcsine form used by Genz (2004):
#   Phi2(a, b; rho) = Phi(a) Phi(b) +
#     1/(2 pi) int_0^asin(rho) exp(-(a^2 + b^2 - 2ab sin t) / (2 cos^2 t)) dt
.pbvnorm <- function(a, b, rho) {
  base <- stats::pnorm(a) * stats::pnorm(b)
  if (rho == 0) return(base)
  fin <- is.finite(a) & is.finite(b)
  if (!any(fin)) return(base)
  th <- asin(rho)
  t <- th / 2 * (.gl20$x + 1)                       # nodes on [0, asin(rho)]
  st <- sin(t)
  ct2 <- cos(t)^2
  af <- a[fin]
  bf <- b[fin]
  E <- exp(-(outer(af^2 + bf^2, rep(1, length(t))) - 2 * outer(af * bf, st)) /
             rep(2 * ct2, each = length(af)))
  base[fin] <- base[fin] + (th / 2) * drop(E %*% .gl20$w) / (2 * pi)
  pmin(pmax(base, 0), 1)
}

# Category index (1..K) of discrete pseudo-observations u given the cut points
# c(0, cumulative probabilities) of the margin
.cell_index <- function(u, cuts) {
  pmin(pmax(findInterval(u, cuts, rightmost.closed = TRUE), 1L), length(cuts) - 1L)
}

.rho_bounds <- c(-0.995, 0.995)

# Polyserial correlation: z continuous normal scores, k category index with
# thresholds tau (length K + 1, from -Inf to Inf). When `bin` is given (an
# integer bin index of each z, from fine quantile bins of the whole column),
# the likelihood is evaluated on the bin x category table with each bin
# represented by its mean score among these rows: each evaluation then costs
# O(bins * K) instead of O(n), and the within-bin spread of z is negligible.
.polyserial <- function(z, k, tau, bin = NULL) {
  K <- length(tau) - 1L
  if (!is.null(bin)) {
    B <- max(bin)
    cnt <- tabulate(bin, B)
    zb <- numeric(B)
    zb[cnt > 0] <- as.vector(rowsum(z, bin, reorder = TRUE))
    zb <- zb / pmax(cnt, 1L)
    N <- matrix(tabulate(bin + B * (k - 1L), B * K), B, K)
    pos <- which(N > 0)
    zz <- zb[(pos - 1L) %% B + 1L]
    kk <- (pos - 1L) %/% B + 1L
    w <- N[pos]
  } else {
    zz <- z
    kk <- k
    w <- 1
  }
  lo <- tau[kk]
  hi <- tau[kk + 1L]
  nll <- function(r) {
    s <- sqrt(1 - r^2)
    p <- stats::pnorm((hi - r * zz) / s) - stats::pnorm((lo - r * zz) / s)
    -sum(w * log(pmax(p, 1e-300)))
  }
  stats::optimize(nll, .rho_bounds, tol = 1e-5)$minimum
}

# Polychoric correlation from the K1 x K2 table of counts N, thresholds t1, t2.
# Only the corners of non-empty cells enter the likelihood, so the bivariate
# normal CDF is evaluated at those (unique) corners only; this keeps the cost
# proportional to the number of observed cells for high-cardinality margins.
.polychoric <- function(N, t1, t2) {
  K1 <- nrow(N)
  pos <- which(N > 0)
  a <- (pos - 1L) %% K1 + 1L                     # row (category of margin 1)
  b <- (pos - 1L) %/% K1 + 1L                    # column (category of margin 2)
  w <- N[pos]
  # corner ids on the (K1 + 1) x (K2 + 1) grid: (i, j) -> i + (K1 + 1) * (j - 1)
  cid <- function(i, j) i + (K1 + 1L) * (j - 1L)
  c11 <- cid(a + 1L, b + 1L)
  c01 <- cid(a, b + 1L)
  c10 <- cid(a + 1L, b)
  c00 <- cid(a, b)
  ids <- unique(c(c11, c01, c10, c00))
  gi <- (ids - 1L) %% (K1 + 1L) + 1L
  gj <- (ids - 1L) %/% (K1 + 1L) + 1L
  A <- t1[gi]
  B <- t2[gj]
  m11 <- match(c11, ids)
  m01 <- match(c01, ids)
  m10 <- match(c10, ids)
  m00 <- match(c00, ids)
  nll <- function(r) {
    F <- .pbvnorm(A, B, r)
    P <- F[m11] - F[m01] - F[m10] + F[m00]
    -sum(w * log(pmax(P, 1e-300)))
  }
  stats::optimize(nll, .rho_bounds, tol = 1e-5)$minimum
}

# Maximum-likelihood latent correlation for every pair involving at least
# one discrete column. `cuts` is a list (one element per column of U): NULL
# for continuous columns, c(0, cumsum(probs)) for discrete ones. Only pairs
# flagged in `todo` are estimated; pairs whose jointly observed rows show a
# single category of a discrete margin carry no information and are returned
# as NA (unidentified). All jointly observed rows are used; for continuous
# columns with more than `exact_max` observed values the polyserial
# likelihood is evaluated on `bins` quantile bins of the column.
.discrete_latent_corr <- function(U, cuts, todo, bins = 500L, exact_max = 1000L) {
  d <- ncol(U)
  out <- matrix(NA_real_, d, d)
  disc <- !vapply(cuts, is.null, logical(1))
  kidx <- lapply(seq_len(d), function(j) {
    if (disc[j]) .cell_index(U[, j], cuts[[j]]) else NULL
  })
  tau <- lapply(cuts, function(cc) if (is.null(cc)) NULL else stats::qnorm(cc))
  Zc <- lapply(seq_len(d), function(j) if (disc[j]) NULL else stats::qnorm(U[, j]))
  bin <- lapply(seq_len(d), function(j) {
    if (disc[j]) return(NULL)
    ok <- !is.na(U[, j])
    if (sum(ok) <= exact_max) return(NULL)
    b <- rep(NA_integer_, nrow(U))
    b[ok] <- as.integer(ceiling(rank(U[ok, j], ties.method = "first") * bins / sum(ok)))
    b
  })
  ij <- which(todo & upper.tri(todo), arr.ind = TRUE)
  for (r in seq_len(nrow(ij))) {
    i <- ij[r, 1L]
    j <- ij[r, 2L]
    rows <- which(!is.na(U[, i]) & !is.na(U[, j]))
    if (disc[i] && disc[j]) {
      K1 <- length(cuts[[i]]) - 1L
      K2 <- length(cuts[[j]]) - 1L
      N <- matrix(tabulate(kidx[[i]][rows] + K1 * (kidx[[j]][rows] - 1L), K1 * K2), K1, K2)
      if (sum(rowSums(N) > 0) < 2L || sum(colSums(N) > 0) < 2L) next
      est <- .polychoric(N, tau[[i]], tau[[j]])
    } else {
      cj <- if (disc[i]) j else i                    # continuous column
      dj <- if (disc[i]) i else j                    # discrete column
      k <- kidx[[dj]][rows]
      if (length(unique(k)) < 2L) next
      est <- .polyserial(Zc[[cj]][rows], k, tau[[dj]],
                         bin = if (is.null(bin[[cj]])) NULL else bin[[cj]][rows])
    }
    out[i, j] <- out[j, i] <- est
  }
  out
}
