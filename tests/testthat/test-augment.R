test_that("augment_data() appends synthetic rows to tabular data", {
  train <- mtcars[1:20, ]
  big <- augment_data(train, n = 100, seed = 1)
  expect_equal(nrow(big), 120L)
  expect_equal(big[1:20, ], train, ignore_attr = TRUE)
  expect_named(big, names(train))
  expect_equal(nrow(augment_data(train, seed = 1)), 220L)   # default n = 10 x
  lab <- augment_data(train, n = 30, label = "synthetic", seed = 1)
  expect_equal(sum(lab$synthetic), 30L)
  expect_false(any(lab$synthetic[1:20]))
  only <- augment_data(train, n = 30, keep_real = FALSE, label = "s", seed = 1)
  expect_equal(nrow(only), 30L)
  expect_true(all(only$s))
  expect_identical(augment_data(train, n = 10, seed = 2), augment_data(train, n = 10, seed = 2))
  strat <- augment_data(iris, n = 300, by = "Species", seed = 1)
  expect_equal(nrow(strat), 450L)
  m <- augment_data(as.matrix(mtcars[, 1:3]), n = 5, seed = 1)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(37L, 3L))
  expect_error(augment_data(train, label = "mpg"), "new column")
  expect_error(augment_data(train, keep_real = NA), "keep_real")
})

test_that("augment_data() returns lists of series for ts input", {
  s <- augment_data(ldeaths, nsim = 3, seed = 1)
  expect_length(s, 4L)
  expect_identical(s[[1]], ldeaths)
  expect_true(all(vapply(s, inherits, logical(1), "ts")))
  s2 <- augment_data(ldeaths, nsim = 2, n = 30, keep_real = FALSE, order = 12, seed = 1)
  expect_length(s2, 2L)
  expect_length(s2[[1]], 30L)
  expect_error(augment_data(ldeaths, nsim = 0), "positive")
})

test_that("fixed VAR order is respected and validated", {
  fit <- fit_synthesizer(ldeaths, order = 12)
  expect_equal(fit$var$order, 12L)
  syn <- generate(fit, n = 3000, seed = 1)
  expect_gt(acf(syn, lag.max = 12, plot = FALSE)$acf[13], 0.3)
  expect_error(fit_synthesizer(ldeaths, order = 40), "too large")
})

test_that("ts_windows() builds correct lag features", {
  x <- ts(1:10)
  w <- ts_windows(x, lags = 3)
  expect_named(w, c("y", "lag1", "lag2", "lag3"))
  expect_equal(nrow(w), 7L)
  expect_equal(unlist(w[1, ]), c(y = 4, lag1 = 3, lag2 = 2, lag3 = 1))
  h <- ts_windows(x, lags = 2, horizon = 3)
  expect_equal(unlist(h[1, ]), c(y = 5, lag1 = 2, lag2 = 1))
  l <- ts_windows(list(x, ts(11:15)), lags = 2)
  expect_equal(nrow(l), 8L + 3L)
  expect_equal(l$series, rep(1:2, c(8, 3)))
  expect_false(any(l$y == 11 | l$y == 12))      # no window crosses series
  mv <- ts_windows(EuStockMarkets[1:20, ], lags = 2, target = "SMI")
  expect_equal(mv$y, EuStockMarkets[3:20, "SMI"])
  expect_true(all(c("DAX_lag1", "FTSE_lag2") %in% names(mv)))
  expect_error(ts_windows(x, lags = 10), "too short")
  expect_error(ts_windows(EuStockMarkets, 2, target = "XX"), "Unknown target")
  expect_error(ts_windows(EuStockMarkets, 2, target = 9), "Invalid")
  expect_error(ts_windows(letters, 2), "numeric")
})

test_that("augment_data() passes closeness and perturb to generation", {
  train <- mtcars[1:20, ]
  near <- augment_data(train, n = 2000, keep_real = FALSE, seed = 1)
  far <- augment_data(train, n = 2000, keep_real = FALSE, closeness = 0, seed = 1)
  r <- function(d) cor(d$mpg, d$wt, method = "spearman")
  expect_gt(abs(r(near)), abs(r(far)) + 0.3)
  shifted <- augment_data(train, n = 2000, keep_real = FALSE, seed = 1,
                          perturb = perturbation(shift = c(mpg = 2)))
  expect_gt(mean(shifted$mpg), mean(train$mpg) + sd(train$mpg))
})

test_that("ts_windows() handles horizons and lists consistently", {
  s <- augment_data(ldeaths, nsim = 3, seed = 1)
  w <- ts_windows(s, lags = 6, horizon = 2)
  expect_equal(nrow(w), 4L * (72L - 6L - 2L + 1L))
  expect_equal(sort(unique(w$series)), 1:4)
  first <- w[w$series == 1, ][1, ]
  expect_equal(first$y, as.numeric(ldeaths[8]))
  expect_equal(first$lag1, as.numeric(ldeaths[6]))
  expect_equal(first$lag6, as.numeric(ldeaths[1]))
})
