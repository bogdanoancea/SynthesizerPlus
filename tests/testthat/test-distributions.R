# Statistical (not only structural) checks of the generators: empirical
# moments and Kendall's tau are compared with their theoretical values within
# Monte Carlo tolerance. Seeds are fixed, so the tests are deterministic.

ktau <- function(u, i = 1, j = 2) cor(u[, i], u[, j], method = "kendall")

test_that("r_dirichlet() reproduces the Dirichlet means and variances", {
  for (alpha in list(c(2, 3, 5), c(0.5, 1.5, 4), c(0.05, 0.1, 0.2))) {
    x <- r_dirichlet(50000, alpha, seed = 11)
    a0 <- sum(alpha)
    m <- alpha / a0
    v <- alpha * (a0 - alpha) / (a0^2 * (a0 + 1))
    expect_equal(unname(colMeans(x)), m, tolerance = 0.01)            # 1% relative
    expect_equal(unname(apply(x, 2, var)), v, tolerance = 0.05)       # 5% relative
    # covariance: -alpha_i alpha_j / (a0^2 (a0 + 1))
    expect_equal(cov(x[, 1], x[, 2]), -alpha[1] * alpha[2] / (a0^2 * (a0 + 1)),
                 tolerance = 0.08)
  }
})

test_that("r_dirichlet() does not underflow for tiny concentrations", {
  x <- r_dirichlet(200, c(1e-300, 1e-300), seed = 1)
  expect_true(all(is.finite(x)))
  expect_equal(unname(rowSums(x)), rep(1, 200))
  # the limit concentrates on the vertices of the simplex
  expect_true(all(x %in% c(0, 1)))
  expect_equal(mean(x[, 1]), 0.5, tolerance = 0.2)
  y <- r_dirichlet(500, c(1e-3, 2, 1e-200), seed = 2)
  expect_true(all(is.finite(y)))
  expect_equal(unname(rowSums(y)), rep(1, 500))
})

test_that("Gaussian and t copulas satisfy tau = 2/pi asin(rho) for every pair", {
  R <- make_corr(3, "toeplitz", rho = c(0.5, -0.2))
  for (fam in c("gaussian", "t")) {
    u <- r_copula(5000, fam, corr = R, df = 4, seed = 3)
    for (p in list(c(1, 2), c(1, 3), c(2, 3))) {
      expect_lt(abs(ktau(u, p[1], p[2]) - 2 / pi * asin(R[p[1], p[2]])), 0.03)
    }
  }
})

test_that("Archimedean copulas match tau(theta) across parameters and pairs", {
  debye1 <- function(t) integrate(function(s) s / expm1(s), 0, t)$value / t
  tau_frank <- function(th) 1 - 4 / th * (1 - debye1(th))
  cases <- list(
    list(fam = "clayton", theta = c(0.5, 2, 6),   tau = function(th) th / (th + 2)),
    list(fam = "gumbel",  theta = c(1.25, 2, 4),  tau = function(th) 1 - 1 / th),
    list(fam = "frank",   theta = c(1, 5.74, 12), tau = tau_frank)
  )
  for (cs in cases) for (th in cs$theta) {
    u <- r_copula(4000, cs$fam, dim = 3, theta = th, seed = 5)
    for (p in list(c(1, 2), c(1, 3), c(2, 3))) {
      expect_lt(abs(ktau(u, p[1], p[2]) - cs$tau(th)), 0.03)
    }
  }
  # negative dependence (bivariate Frank)
  u <- r_copula(4000, "frank", dim = 2, theta = -5.74, seed = 6)
  expect_lt(abs(ktau(u) - tau_frank(-5.74)), 0.03)
  # Gumbel theta = 1 and the independence copula: tau = 0
  expect_lt(abs(ktau(r_copula(4000, "gumbel", dim = 2, theta = 1, seed = 7))), 0.03)
  expect_lt(abs(ktau(r_copula(4000, "independence", dim = 2, seed = 7))), 0.03)
})

test_that("all copula families have uniform margins", {
  fams <- list(list("clayton", theta = 2), list("gumbel", theta = 2),
               list("frank", theta = 5), list("t", corr = make_corr(2, rho = 0.6), df = 3))
  for (f in fams) {
    u <- do.call(r_copula, c(list(3000, f[[1]], dim = 2, seed = 8), f[-1]))
    for (j in 1:2) {
      expect_gt(suppressWarnings(ks.test(u[, j], "punif")$p.value), 0.001)
    }
  }
})

test_that("r_mvt() has covariance df/(df - 2) * sigma", {
  S <- matrix(c(2, 0.6, 0.6, 1), 2)
  x <- r_mvt(60000, df = 6, mean = c(1, -1), sigma = S, seed = 9)
  expect_equal(unname(colMeans(x)), c(1, -1), tolerance = 0.03)
  expect_equal(unname(cov(x)), 6 / 4 * S, tolerance = 0.06)
})

test_that("r_mvskewnorm() matches the skew-normal mean and covariance", {
  Om <- matrix(c(1, 0.5, 0.5, 4), 2)
  alpha <- c(3, -1)
  xi <- c(1, 2)
  w <- sqrt(diag(Om))
  Ob <- Om / outer(w, w)
  delta <- drop(Ob %*% alpha) / sqrt(1 + drop(t(alpha) %*% Ob %*% alpha))
  mu <- w * delta * sqrt(2 / pi)
  x <- r_mvskewnorm(60000, xi = xi, omega = Om, alpha = alpha, seed = 10)
  expect_equal(unname(colMeans(x)), xi + mu, tolerance = 0.02)
  expect_equal(unname(cov(x)), Om - outer(mu, mu), tolerance = 0.05)
})

test_that("r_mvlnorm() matches the log-normal mean", {
  S <- matrix(c(0.25, 0.1, 0.1, 0.16), 2)
  x <- r_mvlnorm(60000, meanlog = c(0, 1), sigmalog = S, seed = 12)
  expect_equal(unname(colMeans(x)), exp(c(0, 1) + diag(S) / 2), tolerance = 0.02)
})

test_that("r_mvmixture() matches the mixture mean and covariance", {
  w <- c(0.3, 0.7)
  m <- list(c(-2, 0), c(2, 1))
  s <- list(diag(2), matrix(c(1, 0.5, 0.5, 2), 2))
  x <- r_mvmixture(60000, w, m, s, seed = 13)
  mu <- w[1] * m[[1]] + w[2] * m[[2]]
  V <- w[1] * (s[[1]] + tcrossprod(m[[1]])) + w[2] * (s[[2]] + tcrossprod(m[[2]])) -
    tcrossprod(mu)
  expect_equal(unname(colMeans(x)), mu, tolerance = 0.03)
  expect_equal(unname(cov(x)), V, tolerance = 0.04)
})
