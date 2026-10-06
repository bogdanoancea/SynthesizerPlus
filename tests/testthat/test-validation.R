# Statistical validation: data are simulated from known distributions,
# dependence structures, missingness mechanisms, strata and time-series
# processes; synthetic samples are compared with the truth (or, where the
# model can only reproduce what it sees, with the training data) within
# Monte Carlo tolerances. Seeds are fixed, so the tests are deterministic.
# The larger designs are skipped on CRAN to keep the check time short.

# Kendall's tau (tau-a with random tie-breaking, as for the training data),
# on at most 3000 rows: it is quadratic in the sample size
kendall <- function(a, b) {
  if (length(a) > 3000L) {
    i <- sample.int(length(a), 3000L)
    a <- a[i]
    b <- b[i]
  }
  SynthesizerPlus:::.kendall_matrix(cbind(rank(as.numeric(a), ties.method = "random"),
                                          rank(as.numeric(b), ties.method = "random")))[1, 2]
}

# Four variables with Gaussian-copula dependence (AR(1), rho = 0.6) and
# gamma, Poisson, 4-category and log-normal margins
mixed_design <- function(n, seed) {
  R <- make_corr(4, "ar1", 0.6)
  U <- r_copula(n, "gaussian", corr = R, seed = seed)
  list(R = R, data = data.frame(
    g = stats::qgamma(U[, 1], 2, 1),
    p = stats::qpois(U[, 2], 3),
    c = factor(c("a", "b", "c", "d")[findInterval(U[, 3], c(0, 0.1, 0.4, 0.9, 1),
                                                    rightmost.closed = TRUE)]),
    l = stats::qlnorm(U[, 4])
  ))
}

test_that("validation: marginals are recovered (continuous, count, categorical)", {
  skip_on_cran()
  d <- mixed_design(4000, 1)$data
  s <- generate(fit_synthesizer(d, seed = 1), 20000, seed = 2)
  expect_lt(suppressWarnings(ks.test(s$g, "pgamma", 2, 1)$statistic), 0.03)
  expect_lt(suppressWarnings(ks.test(s$l, "plnorm")$statistic), 0.03)
  tvd_pois <- sum(abs(tabulate(s$p + 1, 30) / nrow(s) - stats::dpois(0:29, 3))) / 2
  expect_lt(tvd_pois, 0.03)
  expect_equal(as.numeric(prop.table(table(s$c))), c(0.1, 0.3, 0.5, 0.1), tolerance = 0.03)
})

test_that("validation: latent correlations are unbiased with discrete margins", {
  skip_on_cran()
  des <- mixed_design(4000, 1)
  fit <- fit_synthesizer(des$data, seed = 1)
  # every latent correlation (continuous, count, categorical) within 0.05
  expect_lt(max(abs(copula_correlation(fit) - des$R)), 0.05)
  s <- generate(fit, 20000, seed = 2)
  d <- des$data
  # rank association of synthetic data equals that of the training data,
  # also for pairs with discrete variables (0.2.6: 0.29 instead of 0.35)
  for (p in list(c("g", "p"), c("p", "c"), c("c", "l"), c("g", "l"))) {
    expect_lt(abs(kendall(s[[p[1]]], s[[p[2]]]) - kendall(d[[p[1]]], d[[p[2]]])), 0.03)
  }
})

test_that("validation: binary-binary dependence (tetrachoric) is reproduced", {
  skip_on_cran()
  z <- r_mvnorm(5000, sigma = make_corr(2, rho = 0.7), seed = 3)
  d <- data.frame(a = z[, 1] > 0.5, b = factor(ifelse(z[, 2] > -0.3, "yes", "no")))
  s <- generate(fit_synthesizer(d, seed = 1), 40000, seed = 2)
  pt <- prop.table(table(d$a, d$b))
  ps <- prop.table(table(s$a, s$b))
  expect_lt(max(abs(ps - pt)), 0.015)
  # odds ratio within 15% on the log scale
  or <- function(t) log(t[1, 1] * t[2, 2] / (t[1, 2] * t[2, 1]))
  expect_lt(abs(or(ps) - or(pt)), 0.15)
})

test_that("validation: missing-at-random mechanisms are reproduced", {
  skip_on_cran()
  set.seed(4)
  n <- 6000
  x <- rnorm(n)
  w <- rnorm(n)
  y <- 0.6 * x + 0.3 * w + rnorm(n, sd = 0.7)
  y[runif(n) < plogis(-1 + x - w)] <- NA           # MAR on x and w
  z <- rnorm(n)
  z[runif(n) < 0.2] <- NA                          # MCAR
  d <- data.frame(x, w, y, z)
  s <- generate(fit_synthesizer(d, seed = 1), 40000, seed = 2)
  br <- c(-Inf, -1, 0, 1, Inf)
  pm <- function(D) tapply(is.na(D$y), cut(D$x, br), mean)
  expect_lt(max(abs(pm(s) - pm(d))), 0.04)           # P(missing | x)
  expect_near(cor(s$x, is.na(s$y)), cor(d$x, is.na(d$y)), 0.03)
  expect_near(cor(s$w, is.na(s$y)), cor(d$w, is.na(d$y)), 0.03)
  expect_near(mean(is.na(s$z)), 0.2, 0.015)
  expect_lt(abs(cor(s$x, is.na(s$z))), 0.02)
  # relations among the observed rows and the observed-value distribution
  expect_near(cor(s$x, s$y, use = "complete"), cor(d$x, d$y, use = "complete"), 0.03)
  expect_near(mean(s$x[!is.na(s$y)]), mean(d$x[!is.na(d$y)]), 0.03)
  expect_lt(suppressWarnings(ks.test(s$y, d$y)$statistic), 0.03)
})

test_that("validation: strata keep their own means and dependence", {
  skip_on_cran()
  set.seed(5)
  n <- 6000
  k <- sample(c("A", "B", "C"), n, TRUE, prob = c(0.6, 0.3, 0.1))
  mu <- c(A = 0, B = 3, C = -2)[k]
  rho <- c(A = 0.8, B = -0.6, C = 0)[k]
  e1 <- rnorm(n)
  d <- data.frame(k = k, u = mu + e1, v = rho * e1 + sqrt(1 - rho^2) * rnorm(n))
  s <- generate(fit_synthesizer(d, by = "k", seed = 1), 30000, seed = 2)
  expect_near(as.numeric(prop.table(table(s$k))), c(0.6, 0.3, 0.1), 0.02)
  expect_near(as.numeric(tapply(s$u, s$k, mean)), c(0, 3, -2), 0.06)
  cors <- vapply(split(s, s$k), function(q) cor(q$u, q$v), numeric(1))
  expect_near(unname(cors), c(0.8, -0.6, 0), 0.06)
})

test_that("validation: t copula recovers df and tail dependence", {
  skip_on_cran()
  R <- make_corr(3, rho = 0.5)
  U <- r_copula(5000, "t", corr = R, df = 4, seed = 6)
  fit <- fit_synthesizer(as.data.frame(stats::qnorm(U)), copula = "t", seed = 1)
  expect_true(fit$pooled$copula$df %in% c(3, 4, 5, 6))
  expect_lt(max(abs(copula_correlation(fit) - R)), 0.05)
  s <- generate(fit, 60000, seed = 2)
  us <- apply(s, 2, rank) / (nrow(s) + 1)
  ut <- r_copula(60000, "t", corr = R, df = 4, seed = 7)
  tail_dep <- function(u) mean(u[u[, 1] > 0.99, 2] > 0.99)
  expect_lt(abs(tail_dep(us) - tail_dep(ut)), 0.08)
  # a Gaussian copula has (asymptotically) no tail dependence: the
  # difference is visible
  sg <- generate(fit_synthesizer(as.data.frame(stats::qnorm(U)), seed = 1), 60000, seed = 2)
  expect_lt(tail_dep(apply(sg, 2, rank) / 60001), tail_dep(ut) - 0.05)
})

test_that("validation: VAR dynamics and margins of time series are reproduced", {
  skip_on_cran()
  set.seed(8)
  Phi <- matrix(c(0.6, 0.2, -0.1, 0.5), 2)
  TT <- 3000
  Y <- matrix(0, TT, 2)
  for (t in 2:TT) Y[t, ] <- Phi %*% Y[t - 1, ] + rnorm(2)
  # theoretical lag-0 and lag-1 correlation matrices of the VAR(1)
  G0 <- matrix(solve(diag(4) - kronecker(Phi, Phi), as.vector(diag(2))), 2)
  G1 <- Phi %*% G0
  D <- sqrt(outer(diag(G0), diag(G0)))
  lagcor <- function(X) {                   # [i, j] = Cor(Y_{t,i}, Y_{t-1,j})
    X <- scale(X)
    n <- nrow(X)
    crossprod(X[-1, ], X[-n, ]) / (n - 1)
  }
  for (m in c("copula_var", "stationary")) {
    s <- as.matrix(generate(fit_synthesizer(ts(Y), method = m, seed = 1), 20000, seed = 2))
    # reproduces the training series (sampling error of a persistent VAR of
    # length 3000 is ~0.04, so the comparison with the truth is looser)
    expect_lt(max(abs(cor(s) - cor(Y))), 0.05)
    expect_lt(max(abs(lagcor(s) - lagcor(Y))), 0.05)
    expect_lt(max(abs(cor(s) - G0 / D)), 0.1)
    expect_lt(max(abs(lagcor(s) - G1 / D)), 0.1)
    expect_lt(suppressWarnings(ks.test(s[, 1], Y[, 1])$statistic), 0.04)
  }
  # a non-Gaussian margin (exp of an AR(1)): the dependence is reproduced on
  # the latent scale
  x <- exp(stats::arima.sim(list(ar = 0.8), 3000))
  se <- generate(fit_synthesizer(ts(x), seed = 1), 20000, seed = 2)
  acf1 <- function(v) stats::acf(v, plot = FALSE, lag.max = 1)$acf[2]
  expect_near(acf1(log(se)), acf1(log(x)), 0.03)
  expect_lt(suppressWarnings(ks.test(as.numeric(se), as.numeric(x))$statistic), 0.04)
})

test_that("validation (fast, also on CRAN): MAR missingness is not attenuated", {
  set.seed(9)
  n <- 1500
  x <- rnorm(n)
  y <- 0.6 * x + rnorm(n, sd = 0.8)
  y[runif(n) < plogis(-1 + 1.5 * x)] <- NA
  d <- data.frame(x, y)
  s <- generate(fit_synthesizer(d, seed = 1), 10000, seed = 2)
  # 0.2.6 reproduced only ~60% of this association
  expect_gt(cor(s$x, is.na(s$y)), 0.85 * cor(d$x, is.na(d$y)))
})
