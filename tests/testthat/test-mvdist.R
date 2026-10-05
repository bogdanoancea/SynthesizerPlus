test_that("make_corr() builds valid correlation matrices", {
  expect_equal(make_corr(3, "identity"), diag(3))
  R <- make_corr(4, "ar1", 0.5)
  expect_equal(R[1, 3], 0.25)
  expect_equal(make_corr(3, "toeplitz", c(0.4, 0.1))[1, 3], 0.1)
  E <- make_corr(3, "exchangeable", 0.3, names = c("a", "b", "c"))
  expect_equal(colnames(E), c("a", "b", "c"))
  Rr <- make_corr(5, "random")
  expect_equal(unname(diag(Rr)), rep(1, 5))
  expect_true(min(eigen(Rr)$values) > 0)
  expect_error(make_corr(3, "exchangeable", -0.9), "exchangeable")
  expect_error(make_corr(3, "ar1", 1), "ar1")
  expect_error(make_corr(3, "toeplitz", 0.5), "length")
  expect_error(make_corr(3, "toeplitz", c(0.9, -0.9)), "positive definite")
})

test_that("r_mvnorm() has the right moments", {
  S <- matrix(c(4, 1.2, 1.2, 1), 2)
  x <- r_mvnorm(20000, mean = c(1, -1), sigma = S, seed = 1)
  expect_equal(dim(x), c(20000L, 2L))
  expect_equal(colMeans(x), c(1, -1), tolerance = 0.05, ignore_attr = TRUE)
  expect_equal(cov(x), S, tolerance = 0.05, ignore_attr = TRUE)
  expect_identical(r_mvnorm(5, sigma = S, seed = 2), r_mvnorm(5, sigma = S, seed = 2))
  # positive semi-definite covariance
  ps <- matrix(1, 2, 2)
  y <- r_mvnorm(100, sigma = ps, seed = 1)
  expect_equal(y[, 1], y[, 2])
  expect_error(r_mvnorm(5, sigma = matrix(c(1, 2, 2, 1), 2)), "semi-definite")
  expect_error(r_mvnorm(5, mean = 1:3, sigma = diag(2)), "incompatible")
  expect_error(r_mvnorm(5, sigma = matrix(c(1, 0.5, 0.2, 1), 2)), "symmetric")
})

test_that("r_mvt() has heavier tails and correct scale", {
  x <- r_mvt(50000, df = 5, sigma = diag(2), seed = 1)
  expect_equal(var(x[, 1]), 5 / 3, tolerance = 0.1)
  expect_equal(dim(r_mvt(10, df = Inf, mean = c(1, 2), sigma = diag(2))), c(10L, 2L))
  expect_error(r_mvt(10, df = -1, sigma = diag(2)), "df")
})

test_that("r_mvlnorm() is positive", {
  x <- r_mvlnorm(1000, meanlog = c(0, 1), sigmalog = diag(2) * 0.1, seed = 1)
  expect_true(all(x > 0))
})

test_that("r_mvskewnorm() is skewed in the direction of alpha", {
  skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
  x <- r_mvskewnorm(20000, omega = diag(2), alpha = c(10, 0), seed = 1)
  expect_gt(skew(x[, 1]), 0.7)
  expect_lt(abs(skew(x[, 2])), 0.1)
  y <- r_mvskewnorm(20000, omega = diag(2), alpha = c(0, -10), seed = 1)
  expect_lt(skew(y[, 2]), -0.7)
  # alpha = 0 gives the normal distribution
  z <- r_mvskewnorm(20000, xi = c(1, 1), omega = diag(2), alpha = c(0, 0), seed = 1)
  expect_equal(colMeans(z), c(1, 1), tolerance = 0.03, ignore_attr = TRUE)
  expect_error(r_mvskewnorm(5, omega = diag(2), alpha = 1), "incompatible")
})

test_that("r_dirichlet() lies on the simplex", {
  x <- r_dirichlet(5000, c(a = 1, b = 2, c = 7), seed = 1)
  expect_equal(rowSums(x), rep(1, 5000))
  expect_equal(colnames(x), c("a", "b", "c"))
  expect_equal(colMeans(x), c(0.1, 0.2, 0.7), tolerance = 0.03, ignore_attr = TRUE)
  expect_error(r_dirichlet(5, c(1, -1)), "positive")
})

test_that("r_mvmixture() samples components by weight", {
  x <- r_mvmixture(10000, weights = c(1, 3), means = list(c(-5, -5), c(5, 5)),
                   sigmas = list(diag(2), diag(2)), return_component = TRUE, seed = 1)
  comp <- attr(x, "component")
  expect_equal(mean(comp == 2), 0.75, tolerance = 0.02)
  expect_true(all(x[comp == 1, 1] < 0))
  expect_error(r_mvmixture(10, 1, list(0, 1), list(1)), "same length")
})

test_that("Gaussian and t copulas produce uniform margins", {
  for (fam in c("gaussian", "t", "independence")) {
    u <- r_copula(5000, fam, dim = 3, seed = 1)
    expect_equal(dim(u), c(5000L, 3L))
    expect_true(all(u > 0 & u < 1))
    expect_lt(ks.test(u[, 1], "punif")$statistic, 0.03)
  }
  R <- make_corr(2, "exchangeable", 0.7, names = c("x", "y"))
  u <- r_copula(5000, "gaussian", corr = R, seed = 1)
  expect_equal(colnames(u), c("x", "y"))
  expect_equal(cor(qnorm(u))[1, 2], 0.7, tolerance = 0.03)
})

test_that("Archimedean copulas reproduce Kendall's tau", {
  tau <- function(u) cor(u[, 1], u[, 2], method = "kendall")
  set.seed(1)
  # Clayton: tau = theta / (theta + 2)
  expect_equal(tau(r_copula(3000, "clayton", dim = 2, theta = 2, seed = 1)), 0.5, tolerance = 0.06)
  # Gumbel: tau = 1 - 1/theta
  expect_equal(tau(r_copula(3000, "gumbel", dim = 3, theta = 2, seed = 1)), 0.5, tolerance = 0.06)
  # Frank theta = 5.74 gives tau ~ 0.5
  expect_equal(tau(r_copula(3000, "frank", dim = 3, theta = 5.74, seed = 1)), 0.5, tolerance = 0.06)
  expect_lt(tau(r_copula(3000, "frank", dim = 2, theta = -5.74, seed = 1)), -0.4)
  for (fam in c("clayton", "gumbel", "frank")) {
    u <- r_copula(4000, fam, dim = 3, theta = 3, seed = 2)
    expect_true(all(u >= 0 & u <= 1))
    expect_lt(ks.test(u[, 2], "punif")$statistic, 0.04)
  }
  expect_equal(dim(r_copula(10, "gumbel", dim = 2, theta = 1)), c(10L, 2L))
  expect_equal(dim(r_copula(10, "frank", dim = 3, theta = 0)), c(10L, 3L))
})

test_that("r_copula() validates parameters", {
  expect_error(r_copula(10, "clayton"), "theta")
  expect_error(r_copula(10, "clayton", theta = -1), "theta > 0")
  expect_error(r_copula(10, "gumbel", theta = 0.5), "theta >= 1")
  expect_error(r_copula(10, "frank", dim = 3, theta = -1), "dim > 2")
})

test_that("r_mvdist() combines copulas and margins", {
  d <- r_mvdist(
    2000,
    margins = list(
      inc = list(dist = "lnorm", meanlog = 1),
      age = "norm",
      grp = margin_categorical(c("a", "b"), c(0.2, 0.8)),
      emp = margin_empirical(as.Date("2020-01-01") + 0:99),
      cnt = list(dist = "pois", lambda = 3)
    ),
    copula = "clayton", theta = 2, seed = 1
  )
  expect_s3_class(d, "data.frame")
  expect_equal(dim(d), c(2000L, 5L))
  expect_true(all(d$inc > 0))
  expect_s3_class(d$grp, "factor")
  expect_equal(mean(d$grp == "b"), 0.8, tolerance = 0.03)
  expect_s3_class(d$emp, "Date")
  expect_true(all(d$cnt == round(d$cnt)))
  expect_gt(cor(d$inc, d$age, method = "kendall"), 0.3)
  unnamed <- r_mvdist(5, list("norm", "unif"))
  expect_named(unnamed, c("X1", "X2"))
  expect_error(r_mvdist(5, list()), "non-empty")
  expect_error(r_mvdist(5, list("nodist")), "No quantile function")
  expect_error(r_mvdist(5, list(42)), "Invalid margin")
  expect_error(r_mvdist(5, list("norm", "norm"), corr = diag(3)), "dimension")
})

test_that("margin helpers validate input", {
  q <- margin_categorical(c("l", "m", "h"), ordered = TRUE)
  expect_true(is.ordered(q(c(0.1, 0.5, 0.9))))
  expect_equal(as.character(q(c(0.1, 0.5, 0.9))), c("l", "m", "h"))
  expect_error(margin_categorical(c("a", "b"), c(1, -1)), "non-negative")
  expect_error(margin_empirical(c(NA, NA)), "no observed")
  qe <- margin_empirical(c("x", "y", "y", "y"))
  expect_equal(mean(qe(runif(5000)) == "y"), 0.75, tolerance = 0.03)
})
