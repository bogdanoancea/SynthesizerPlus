test_that("hybrid is opt-in and the old default is unchanged", {
  set.seed(1)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  f0 <- fit_synthesizer(d, seed = 11)
  f1 <- fit_synthesizer(d, nonlinear_target = NULL, seed = 11)
  expect_identical(generate(f0, n = 50, seed = 12),
                   generate(f1, n = 50, seed = 12))
  expect_null(f0$pooled$nonlinear)
})

test_that("hybrid CART generates reproducible numeric responses", {
  skip_if_not_installed("rpart")
  set.seed(20)
  n <- 500
  d <- data.frame(x = runif(n, -2, 2), z = rnorm(n))
  d$y <- sin(3 * d$x) + (d$z > 0) * 2 + 0.2 * rnorm(n)
  f <- fit_synthesizer(d, nonlinear_target = "y",
                       nonlinear_minbucket = 12, seed = 15)
  expect_s3_class(f, "sp_synthesizer")
  expect_identical(f$pooled$nonlinear$target, "y")
  a <- generate(f, n = 800, seed = 16)
  b <- generate(f, n = 800, seed = 16)
  expect_identical(a, b)
  expect_named(a, names(d))
  expect_true(all(is.finite(a$y)))
  expect_equal(nrow(a), 800L)
  expect_lt(abs(stats::median(a$y) - stats::median(d$y)), 0.15)
  expect_error(generate(f, n = 200, dependence = 0.4), "requires dependence")
  expect_error(generate(f, n = 200, closeness = 0.7), "requires dependence")
  expect_error(generate(f, n = 200, perturb = perturbation(noise = 0.2)),
               "requires dependence")
  expect_true(is.null(fit_synthesizer(d, seed = 15)$pooled$nonlinear))
})

test_that("hybrid supports factor predictors and observation missingness", {
  skip_if_not_installed("rpart")
  set.seed(33)
  n <- 400L
  x <- rnorm(n)
  g <- factor(sample(c("alpha", "beta", "gamma"), n, TRUE))
  y <- 1.5 * (g == "gamma") + x^2 + rnorm(n, sd = 0.3)
  x[sample.int(n, 25)] <- NA_real_
  y[sample.int(n, 30)] <- NA_real_
  d <- data.frame(x = x, g = g, y = as.numeric(y))
  f <- fit_synthesizer(d, nonlinear_target = "y", seed = 4,
                       nonlinear_minbucket = 10)
  a <- generate(f, n = 1000, seed = 5)
  expect_s3_class(a$g, "factor")
  expect_true(anyNA(a$x))
  expect_true(anyNA(a$y))
  expect_true(all(is.finite(a$y[!is.na(a$y)])))
  expect_identical(levels(a$g), levels(g))
})

test_that("hybrid validates inputs and stratified models", {
  skip_if_not_installed("rpart")
  set.seed(56)
  d <- data.frame(k = rep(c("A", "B"), each = 100),
                  x = rnorm(200), y = rnorm(200))
  expect_error(fit_synthesizer(d, nonlinear_target = "k"),
               "double-precision numeric")
  expect_error(fit_synthesizer(d, nonlinear_target = "y", nonlinear_minbucket = 0),
               "nonlinear_minbucket")
  expect_error(fit_synthesizer(d, nonlinear_target = "y", nonlinear_cp = -1),
               "nonlinear_cp")
  f <- fit_synthesizer(d, by = "k", nonlinear_target = "y", seed = 12,
                       nonlinear_minbucket = 10)
  s <- generate(f, n = 300, seed = 13)
  expect_equal(nrow(s), 300L)
  expect_true(all(s$k %in% c("A", "B")))
  expect_true(all(is.finite(s$y)))
  expect_identical(generate(f, n = 300, seed = 13), s)
})
