# Regression tests for issues found in the external review of 0.2.0.

test_that("r_copula() rejects matrices that are not correlation matrices", {
  expect_error(r_copula(10, "gaussian", corr = diag(c(4, 1))), "unit diagonal")
  expect_error(r_copula(10, "gaussian", corr = matrix(c(1, 1.2, 1.2, 1), 2)), "\\[-1, 1\\]|semi-definite")
  expect_error(r_copula(10, "gaussian", corr = matrix(c(1, 0.5, 0.2, 1), 2)), "symmetric")
  expect_error(r_copula(10, "gaussian", corr = matrix(c(1, NA, NA, 1), 2)), "finite")
  expect_error(r_copula(10, "gaussian", corr = matrix(1:6, 2)), "square")
  bad_psd <- matrix(c(1, 0.9, -0.9, 0.9, 1, 0.9, -0.9, 0.9, 1), 3)
  expect_error(r_copula(10, "gaussian", corr = bad_psd), "semi-definite")
  # a valid matrix gives uniform margins
  u <- r_copula(20000, "gaussian", corr = make_corr(2, "exchangeable", 0.6), seed = 1)
  expect_lt(ks.test(u[, 1], "punif")$statistic, 0.015)
  # r_mvdist() goes through the same validation
  expect_error(r_mvdist(5, list("norm", "norm"), corr = diag(c(2, 1))), "unit diagonal")
})

test_that("the t copula validates its degrees of freedom", {
  for (bad in list(0, -1, NA_real_, c(3, 4), "4")) {
    expect_error(r_copula(5, "t", dim = 2, df = bad), "df", info = format(bad))
  }
  expect_error(r_mvt(5, df = NA, sigma = diag(2)), "df")
  # df = Inf is the Gaussian copula
  a <- r_copula(200, "t", dim = 2, df = Inf, seed = 3)
  b <- r_copula(200, "gaussian", dim = 2, seed = 3)
  expect_equal(unname(a), unname(b))
  expect_error(r_copula(5, "clayton", dim = 2, theta = Inf), "finite")
})

test_that("row keys distinguish NA from the string '<NA>' and never collide", {
  k <- SynthesizerPlus:::.stratum_key
  expect_length(unique(k(data.frame(x = c(NA, "<NA>")))), 2L)
  sep_clash <- data.frame(a = c("x\ry", "x"), b = c("z", "y\rz"))
  expect_length(unique(k(sep_clash)), 2L)
  dot_clash <- data.frame(a = c("1.2", "1"), b = c("3", "2.3"))
  expect_length(unique(k(dot_clash)), 2L)
  # keys are consistent across data frames, factors and characters alike
  kk <- SynthesizerPlus:::.row_keys(data.frame(f = factor(c("a", "b"))),
                                    data.frame(f = c("b", "c")))
  expect_equal(kk[[1]][2], kk[[2]][1])
  expect_false(kk[[2]][2] %in% kk[[1]])
  # numeric values are matched exactly, not via printed digits
  num <- data.frame(x = c(0.1 + 0.2, 0.3))
  expect_length(unique(k(num)), 2L)
})

test_that("stratification and disclosure risk treat NA and '<NA>' separately", {
  d <- data.frame(g = rep(c(NA, "<NA>"), each = 25), x = c(rnorm(25), rnorm(25, 10)))
  fit <- fit_synthesizer(d, by = "g")
  expect_length(fit$strata$keys, 2L)
  syn <- generate(fit, n = 2000, seed = 1)
  expect_gt(mean(syn$x[syn$g %in% "<NA>"]), 8)
  expect_lt(abs(mean(syn$x[is.na(syn$g)])), 1)
  r <- disclosure_risk(data.frame(k = c(NA, "<NA>", "a", "a"), y = 1:4),
                       data.frame(k = c(NA, "<NA>", "a", "a"), y = 1:4), keys = "k")
  expect_equal(r$n_real_uniques, 2L)
  expect_equal(r$replicated_uniques, 1)
  # stratum labels for display are built from the values
  st <- summary(fit_synthesizer(airquality, by = "Month"))$strata
  expect_equal(st$stratum, as.character(5:9))
})

test_that("Gower distance: missingness as a category or standard Gower", {
  real <- data.frame(a = c(1, NA, 10), b = c(NA, 5, 10), c = c(2, 3, 10),
                     d = c(NA, NA, 10), e = c(NA, NA, 10))
  syn <- data.frame(a = NA_real_, b = NA_real_, c = 2.9, d = NA_real_, e = NA_real_)
  g <- function(na) {
    SynthesizerPlus:::.min_gower(SynthesizerPlus:::.gower_prep(syn, real),
                                 SynthesizerPlus:::.gower_prep(real, real), na = na)
  }
  # category: vs row 2 -> a match, b mismatch, c 0.1/8, d e match -> 1.0125 / 5
  expect_equal(g("category"), (0 + 1 + 0.1 / 8 + 0 + 0) / 5)
  # standard Gower: only c is observed in both -> 0.1 / 8
  expect_equal(g("exclude"), 0.1 / 8)
  # no variable observed in both -> no distance
  none <- SynthesizerPlus:::.min_gower(
    SynthesizerPlus:::.gower_prep(data.frame(a = NA_real_, b = 1), data.frame(a = 1, b = NA_real_)),
    SynthesizerPlus:::.gower_prep(data.frame(a = 1, b = NA_real_), data.frame(a = 1, b = NA_real_)),
    na = "exclude")
  expect_true(is.na(none))
  aq <- synthesize(airquality, seed = 1)
  expect_true(is.finite(dcr(airquality, aq, seed = 1, na = "exclude")$ratio))
  expect_error(dcr(data.frame(a = c(NA, 1), b = c(1, NA)), data.frame(a = c(NA, NA), b = c(1, 1)),
                   na = "exclude"), NA)
  expect_error(dcr(airquality, aq, na = "drop"), "should be one of")
})

test_that("r_mvmixture() validates weights and component dimensions", {
  expect_error(r_mvmixture(3, c(0, 0), list(0, 1), list(matrix(1), matrix(1))), "not all zero")
  expect_error(r_mvmixture(3, c(1, NA), list(0, 1), list(matrix(1), matrix(1))), "finite")
  expect_error(r_mvmixture(3, c(1, -1), list(0, 1), list(matrix(1), matrix(1))), "non-negative")
  expect_error(r_mvmixture(10, c(1, 1), list(c(0, 0), c(0, 0, 0)), list(diag(2), diag(3))),
               "Component 2")
  expect_error(r_mvmixture(10, c(1, 1), list(c(0, 0), c(0, 0)), list(diag(2), diag(3))),
               "Component 2")
  expect_error(r_mvmixture(10, 1, c(0, 0), list(diag(2))), "lists")
  ok <- r_mvmixture(10, c(0, 1), list(c(0, 0), c(5, 5)), list(diag(2), diag(2)), seed = 1)
  expect_true(all(ok > 0))          # zero weight on the first component is allowed
})

test_that("margin helpers validate levels and probabilities", {
  expect_error(margin_categorical(c("a", "a")), "unique")
  expect_error(margin_categorical(c("a", NA)), "unique|non-missing")
  expect_error(margin_categorical(c("a", "b"), c(1, Inf)), "finite")
  q <- margin_categorical(c("a", "b"))
  expect_error(q(c(0.5, 1.2)), "\\[0, 1\\]")
  expect_error(q(NA_real_), "\\[0, 1\\]")
  expect_error(margin_empirical(1:10)(-0.1), "\\[0, 1\\]")
})

test_that("discriminator preprocessing is learned within training folds only", {
  # a variable that is constant within the training folds but not overall
  # cannot be scaled by statistics of the held-out data
  spec <- SynthesizerPlus:::.prep_fit(data.frame(x = c(1, 2, 3, NA), g = c("a", "b", "a", "a")))
  X <- SynthesizerPlus:::.prep_apply(spec, data.frame(x = c(100, NA), g = c("zzz", "b")))
  expect_equal(colnames(X), c("x..NA", "x", "g==b"))
  expect_equal(X[, "x"], c((100 - 2) / sd(c(1, 2, 3, 2)), 0))  # NA imputed by training median
  expect_equal(X[, "g==b"], c(0, 1))                            # unseen level -> reference
  # the AUC is reproducible and in a sensible range
  syn <- synthesize(iris, seed = 1)
  a1 <- discriminator_auc(iris, syn, seed = 2)
  expect_identical(a1, discriminator_auc(iris, syn, seed = 2))
  expect_true(a1 > 0.4 && a1 < 0.85)
  x <- as.data.frame(r_mvnorm(400, sigma = diag(3), seed = 1))
  y <- as.data.frame(r_mvnorm(400, sigma = diag(3), seed = 2))
  expect_lt(abs(discriminator_auc(x, y, seed = 1) - 0.5), 0.06)
  # copies of the real records push the cross-validated AUC below 0.5
  expect_lt(discriminator_auc(iris[sample(150), ], iris, seed = 2), 0.4)
})

test_that("Cramer's V cannot produce NaN from empty categories", {
  # tables are built from observed values only, so expected counts are never 0
  v <- SynthesizerPlus:::.cramers_v(factor(c("a", "b"), levels = c("a", "b", "c")),
                                    factor(c("x", "y"), levels = c("x", "y", "z")))
  expect_false(is.nan(v))
  A <- association_matrix(data.frame(f = factor(c("a", "a", "b"), levels = letters[1:5]),
                                     g = c("x", "y", "y")))
  expect_false(any(is.nan(A)))
})

test_that("date-time strings are read in the right time zone", {
  tmpl <- data.frame(t = as.POSIXct("2020-01-01", tz = "Europe/Bucharest"))
  x <- data.frame(t = c("2020-06-01 12:00:00", "2020-06-01T09:00:00Z",
                        "2020-06-01T12:00:00+03:00", "2020-06-01T05:00:00-0400",
                        "2020/06/01 12:00", NA), stringsAsFactors = FALSE)
  y <- match_types(x, tmpl)$t
  expect_equal(attr(y, "tzone"), "Europe/Bucharest")
  expect_equal(format(y[1:5], "%H:%M", tz = "Europe/Bucharest"), rep("12:00", 5))
  expect_true(is.na(y[6]))
  # the package's own CSV output (UTC with Z) round-trips in any time zone
  d <- data.frame(t = as.POSIXct(c("2021-03-28 02:30:00", "2021-10-31 03:30:00"),
                                 tz = "Europe/Bucharest"))
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  write_data(d, f)
  expect_equal(read_data(f, template = d)$t, d$t)
})
