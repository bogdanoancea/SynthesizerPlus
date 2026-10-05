test_that("copula_var reproduces marginals and autocorrelation", {
  set.seed(1)
  x <- arima.sim(list(ar = 0.8), n = 600)
  x <- ts(exp(x / 3), frequency = 4, start = c(2000, 1))
  fit <- fit_synthesizer(x)
  expect_s3_class(fit, "sp_ts_synthesizer")
  expect_output(print(fit), "VAR order")
  syn <- generate(fit, n = 5000, seed = 1)
  expect_s3_class(syn, "ts")
  expect_null(dim(syn))
  expect_equal(frequency(syn), 4)
  expect_equal(start(syn), c(2000, 1))
  a_real <- acf(x, plot = FALSE)$acf[2]
  a_syn <- acf(syn, plot = FALSE)$acf[2]
  expect_lt(abs(a_real - a_syn), 0.1)
  expect_true(all(syn >= min(x) & syn <= max(x)))
  expect_lt(abs(median(syn) - median(x)), 0.1 * sd(x))
})

test_that("multivariate copula_var keeps cross-correlation", {
  r <- diff(log(EuStockMarkets))
  fit <- fit_synthesizer(r)
  syn <- generate(fit, seed = 1)
  expect_equal(dim(syn), dim(r))
  expect_equal(colnames(syn), colnames(r))
  expect_lt(max(abs(cor(syn) - cor(r))), 0.1)
})

test_that("block and stationary bootstraps generate valid series", {
  for (m in c("block", "stationary")) {
    fit <- fit_synthesizer(EuStockMarkets, method = m, block_length = 20)
    expect_equal(fit$block_length, 20L)
    syn <- generate(fit, n = 333, seed = 1)
    expect_equal(dim(syn), c(333L, 4L))
    expect_true(all(syn[, "DAX"] %in% EuStockMarkets[, "DAX"]))
    expect_output(print(fit), "block length")
  }
  # default block length
  expect_equal(fit_synthesizer(Nile, method = "block")$block_length, 5L)
})

test_that("block bootstrap preserves serial dependence", {
  fit <- fit_synthesizer(lh, method = "stationary", block_length = 10)
  syn <- generate(fit, n = 5000, seed = 2)
  expect_lt(abs(acf(syn, plot = FALSE)$acf[2] - acf(lh, plot = FALSE)$acf[2]), 0.15)
})

test_that("iid method ignores serial dependence", {
  fit <- fit_synthesizer(ldeaths, method = "iid")
  syn <- generate(fit, n = 2000, seed = 1)
  expect_lt(abs(acf(syn, plot = FALSE)$acf[2]), 0.1)
  expect_s3_class(syn, "ts")
})

test_that("time series can have any length and a custom start", {
  fit <- fit_synthesizer(Nile)
  expect_length(generate(fit, n = 1500, seed = 1), 1500L)
  s <- generate(fit, n = 10, start = 2001)
  expect_equal(start(s)[1], 2001)
  expect_identical(generate(fit, seed = 3), generate(fit, seed = 3))
})

test_that("missing values are interpolated with a warning", {
  x <- Nile
  x[c(5, 50)] <- NA
  expect_warning(fit <- fit_synthesizer(x), "interpolated")
  expect_false(anyNA(generate(fit, seed = 1)))
})

test_that("time series input is validated", {
  expect_error(fit_synthesizer(ts(1:5)), "At least 10")
  expect_error(generate(fit_synthesizer(Nile), n = 0), "positive")
  sims <- simulate(fit_synthesizer(Nile), nsim = 2, seed = 1)
  expect_length(sims, 2L)
})
