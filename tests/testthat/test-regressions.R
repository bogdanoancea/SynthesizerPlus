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
  expect_equal(g("category"), (0 + 1 + 0.1 / 8 + 0 + 0) / 5, ignore_attr = TRUE)
  expect_equal(attr(g("category"), "compared"), 1)
  # standard Gower: only c is observed in both -> 0.1 / 8, on 1 of 5 variables
  expect_equal(g("exclude"), 0.1 / 8, ignore_attr = TRUE)
  expect_equal(attr(g("exclude"), "compared"), 1 / 5)
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

# ---- issues from the second review (0.2.1) -------------------------------------

test_that("r_dirichlet() rejects non-finite or non-positive alpha", {
  for (a in list(c(1, NA), c(1, NaN), c(1, Inf), c(1, 0), c(1, -1), 1, "a")) {
    expect_error(r_dirichlet(10, a), "finite positive", info = format(a))
  }
})

test_that("disclosure_risk() requires a finite tolerance", {
  x <- data.frame(k = 1:3, y = 1:3)
  for (tol in list(NA_real_, Inf, -1, c(0.1, 0.2), "0.1")) {
    expect_error(disclosure_risk(x, x, keys = "k", target = "y", tolerance = tol),
                 "tolerance", info = format(tol))
  }
  expect_s3_class(disclosure_risk(x, x, keys = "k", target = "y", tolerance = 0), "sp_disclosure")
})

test_that("seeds must be finite and in the integer range", {
  for (s in list(Inf, -Inf, NaN, 1e12, c(1, 2), "1")) {
    expect_error(r_copula(10, seed = s), "seed", info = format(s))
  }
  expect_error(generate(fit_synthesizer(iris), seed = Inf), "seed")
  expect_silent(r_copula(2, seed = -.Machine$integer.max))
})

test_that("make_corr() validates rho", {
  for (type in c("exchangeable", "ar1")) {
    for (r in list(NA_real_, NaN, Inf, "a", numeric(0))) {
      expect_error(make_corr(3, type, r), "rho", info = paste(type, format(r)))
    }
  }
  expect_error(make_corr(3, "toeplitz", c(NA, 0)), "rho")
  expect_equal(dim(make_corr(3, "random", rho = NA)), c(3L, 3L))   # rho unused for random
})

test_that("all mixture components are validated before sampling", {
  bad_sigma <- matrix(c(1, 2, 2, 1), 2)
  expect_error(r_mvmixture(10, c(1, 0), list(c(0, 0), c(0, 0)), list(diag(2), bad_sigma)),
               "Component 2.*semi-definite")
  expect_error(r_mvmixture(10, c(1, 0), list(c(0, 0), c("a", "b")), list(diag(2), diag(2))),
               "Component 2.*finite")
  expect_error(r_mvmixture(10, c(1, 1), list(c(0, 0), c(Inf, 0)), list(diag(2), diag(2))),
               "Component 2.*finite")
  expect_error(r_mvmixture(10, c(1, 1), list(c(0, 0), c(0, 0)),
                           list(diag(2), matrix(c(1, NA, NA, 1), 2))), "Component 2.*finite")
  expect_error(r_mvmixture(10, c(1, 1), list(c(0, 0), c(0, 0)),
                           list(diag(2), matrix(c(1, 0.5, 0.2, 1), 2))), "Component 2.*symmetric")
})

test_that("row keys handle special values consistently", {
  k <- SynthesizerPlus:::.stratum_key
  num <- k(data.frame(v = c(NA, NaN, Inf, -Inf, -0, 0, 1)))
  expect_equal(num[1], num[2])              # NA and NaN are both missing
  expect_false(num[3] == num[4])            # Inf and -Inf differ
  expect_equal(num[5], num[6])              # -0 equals 0
  expect_false(num[1] == num[6])
  chr <- k(data.frame(v = c(NA_character_, "NA", "<NA>", "")))
  expect_length(unique(chr), 4L)
  int <- k(data.frame(v = c(NA_integer_, 0L, NA_integer_)))
  expect_equal(int[1], int[3])
  lgl <- k(data.frame(v = c(NA, TRUE, FALSE, NA)))
  expect_equal(lgl[1], lgl[4])
  expect_length(unique(lgl), 3L)
  dt <- k(data.frame(v = as.Date(c(NA, "2020-01-01", "2020-01-01"))))
  expect_equal(dt[2], dt[3])
  # consistency across data frames (real vs synthetic) for the same values
  kk <- SynthesizerPlus:::.row_keys(data.frame(v = c(NaN, 1)), data.frame(v = c(NA, 1)))
  expect_equal(kk[[1]], kk[[2]])
})

test_that("exact matches under na = 'exclude' require all variables", {
  real <- data.frame(a = c(1, 2), b = c(5, 6), c = c(7, 8))
  partial <- data.frame(a = c(1, NA), b = c(NA, NA), c = c(NA, 8))   # 1 of 3 variables each
  d <- dcr(real, partial, na = "exclude")
  expect_equal(d$zero_distance_rate, 1)
  expect_equal(d$exact_match_rate, 0)
  expect_equal(d$compared_share, 1 / 3)
  full <- dcr(real, real[2:1, ], na = "exclude")
  expect_equal(full$exact_match_rate, 1)
  expect_equal(full$compared_share, 1)
  cat_d <- dcr(real, partial, na = "category")
  expect_equal(cat_d$compared_share, 1)
  expect_equal(cat_d$exact_match_rate, cat_d$zero_distance_rate)
})

test_that("comparison output flags an AUC well below 0.5", {
  cmp <- compare_synthetic(iris, iris[sample(150), ], metrics = "utility", seed = 2)
  expect_lt(cmp$utility$auc, 0.4)
  expect_output(print(cmp), "copy real")
  ok <- compare_synthetic(iris, synthesize(iris, seed = 1), metrics = "utility", seed = 2)
  expect_false(any(grepl("copy real", capture.output(print(ok)))))
})

# ---- 0.2.3: third external review -----------------------------------------

test_that("covariance / scale matrices are validated with clear messages", {
  for (f in list(function(S) r_mvnorm(10, sigma = S),
                 function(S) r_mvt(10, df = 5, sigma = S),
                 function(S) r_mvlnorm(10, sigmalog = S),
                 function(S) r_mvskewnorm(10, omega = S, alpha = c(1, 1)))) {
    expect_error(f(matrix(1:6, 2, 3)), "square")
    expect_error(f(matrix(numeric(0), 0, 0)), "square")
    expect_error(f(matrix("a", 2, 2)), "numeric")
    expect_error(f(matrix(c(1, NA, NA, 1), 2)), "finite")
    expect_error(f(matrix(c(1, Inf, Inf, 1), 2)), "finite")
    expect_error(f(matrix(c(1, 0.5, 0.2, 1), 2)), "symmetric")
    expect_error(f(matrix(c(1, 2, 2, 1), 2)), "semi-definite")
  }
  expect_error(r_mvmixture(10, 1, list(c(0, 0)), list(matrix(c(1, 0.5, 0.2, 1), 2))),
               "Component 1: 'sigma' must be symmetric")
  # singular but PSD covariances remain allowed for the normal family
  x <- r_mvnorm(50, sigma = matrix(1, 2, 2), seed = 1)
  expect_equal(x[, 1], x[, 2])
})

test_that("multivariate generators reject non-finite locations and shapes", {
  expect_error(r_mvnorm(10, mean = c(Inf, 0), sigma = diag(2)), "finite")
  expect_error(r_mvnorm(10, mean = c("a", "b"), sigma = diag(2)), "numeric")
  expect_error(r_mvnorm(10, mean = 1:3, sigma = diag(2)), "incompatible dimensions")
  expect_error(r_mvt(10, df = 5, mean = c(NA, 0), sigma = diag(2)), "finite")
  expect_error(r_mvlnorm(10, meanlog = c(NaN, 0), sigmalog = diag(2)), "finite")
  expect_error(r_mvskewnorm(10, xi = c(0, 0), omega = diag(2), alpha = c(Inf, 1)), "'alpha'.*finite")
  expect_error(r_mvskewnorm(10, xi = c(NA, 0), omega = diag(2), alpha = c(1, 1)), "'xi'.*finite")
  expect_error(r_mvskewnorm(10, omega = diag(2), alpha = 1), "'alpha'")
  expect_error(r_mvmixture(10, 1, list(c(0, Inf)), list(diag(2))), "Component 1.*finite")
  expect_error(r_mvt(-1, df = 5, sigma = diag(2)), "non-negative")
})

test_that("skew-normal requires a strictly positive scale diagonal", {
  S <- matrix(c(0, 0, 0, 1), 2, 2)          # PSD, but zero scale
  expect_error(r_mvskewnorm(10, omega = S, alpha = c(1, 1)), "positive")
})

test_that("make_corr(d = 1) is the 1 x 1 identity whatever rho", {
  expect_identical(make_corr(1, "exchangeable", rho = -5), matrix(1, 1, 1))
  expect_identical(make_corr(1, "ar1", rho = 0.5, names = "a"),
                   matrix(1, 1, 1, dimnames = list("a", "a")))
  expect_identical(make_corr(1, "random"), matrix(1, 1, 1))
})

test_that("t-copula df fallback is flagged and warned about", {
  set.seed(1)
  small <- data.frame(x = rnorm(8), y = rnorm(8))
  expect_warning(s <- fit_synthesizer(small, copula = "t"), "df was not estimated")
  expect_identical(s$pooled$copula$df, 10)
  expect_false(s$pooled$copula$df_estimated)
  expect_output(print(summary(s)), "fixed: too few complete rows")

  big <- data.frame(x = rnorm(200), y = rnorm(200))
  expect_silent(s2 <- fit_synthesizer(big, copula = "t"))
  expect_true(s2$pooled$copula$df_estimated)
  expect_true(s2$pooled$copula$df %in% c(2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100))
  expect_true(is.na(fit_synthesizer(big)$pooled$copula$df_estimated))
})

# ---- 0.2.4: fourth external review ----------------------------------------

test_that(".check_prob() rejects every non-finite value", {
  for (v in list(NaN, NA_real_, Inf, -Inf, NA)) {
    expect_error(SynthesizerPlus:::.check_prob(v, "p"), "single number")
  }
  expect_identical(SynthesizerPlus:::.check_prob(0.3, "p"), 0.3)
})

split_modules <- function(n = 300, rho = 0.8, seed = 1) {
  z <- r_mvnorm(n, sigma = make_corr(3, rho = rho), seed = seed)
  d <- data.frame(a = z[, 1], b = z[, 2], c = z[, 3])
  h <- n %/% 2
  d$a[seq_len(h)] <- NA                       # a and b never observed together
  d$b[(h + 1):n] <- NA
  d
}

test_that("unidentified pairwise correlations are completed, counted and reported", {
  d <- split_modules()
  for (fam in c("gaussian", "t")) {
    w <- testthat::capture_warnings(fit <- fit_synthesizer(d, copula = fam, missing = "drop"))
    expect_true(any(grepl("not identified", w)))
    # no complete rows at all, so the t df is the documented fallback
    if (fam == "t") expect_true(any(grepl("df was not estimated", w)))
    np <- copula_correlation(fit, what = "n_pair")
    expect_identical(np["a", "b"], 0L)
    expect_identical(np["a", "c"], 150L)
    expect_identical(unname(diag(np)), c(150L, 150L, 300L))
    expect_identical(fit$pooled$copula$unidentified, 1L)
    R <- copula_correlation(fit)
    # max-determinant completion: conditional independence of a and b given c,
    # i.e. r_ab = r_ac * r_bc, instead of the arbitrary value 0
    expect_equal(R["a", "b"], R["a", "c"] * R["b", "c"], tolerance = 1e-6)
    expect_gt(R["a", "b"], 0.4)
    expect_gt(min(eigen(R, only.values = TRUE)$values), 0)
  }
  expect_output(print(summary(suppressWarnings(fit_synthesizer(d, missing = "drop")))),
                "Unidentified pairwise correlations")
  # fully observed data: no warning, nothing unidentified
  expect_silent(f2 <- fit_synthesizer(iris))
  expect_identical(f2$pooled$copula$unidentified, 0L)
  expect_true(all(copula_correlation(f2, what = "n_pair") == 150L))
})

test_that("a pair seen together in only 2 rows does not distort identified pairs", {
  d <- split_modules()
  ref <- suppressWarnings(fit_synthesizer(d, missing = "drop"))
  z <- r_mvnorm(300, sigma = make_corr(3, rho = 0.8), seed = 1)
  d$b[151:152] <- z[151:152, 2]                # 2 joint rows: |r| = 1 trivially
  expect_warning(fit <- fit_synthesizer(d, missing = "drop"), "fewer than 3 rows")
  R0 <- copula_correlation(ref)
  R1 <- copula_correlation(fit)
  expect_equal(R1["a", "c"], R0["a", "c"], tolerance = 0.02)
  expect_equal(R1["a", "b"], R1["a", "c"] * R1["b", "c"], tolerance = 1e-6)
})

test_that("completion keeps identified entries and falls back to zero if impossible", {
  cc <- function(R, unid) SynthesizerPlus:::.complete_corr(R, unid)$R
  R <- matrix(c(1, NA, 0.5, NA, 1, 0.4, 0.5, 0.4, 1), 3)
  unid <- is.na(R)
  W <- cc(R, unid)
  expect_equal(W[1, 3], 0.5)
  expect_equal(W[2, 3], 0.4)
  expect_equal(W[1, 2], 0.2, tolerance = 1e-8)
  expect_equal(solve(W)[1, 2], 0, tolerance = 1e-8)   # zero partial correlation
  # 4 x 4 cycle a-b-c-d-a with a-c, b-d unknown
  R4 <- diag(4)
  R4[cbind(c(1, 2, 3, 4), c(2, 3, 4, 1))] <- 0.5
  R4 <- pmax(R4, t(R4))
  u4 <- matrix(FALSE, 4, 4)
  u4[cbind(c(1, 3, 2, 4), c(3, 1, 4, 2))] <- TRUE
  R4[u4] <- NA
  W4 <- cc(R4, u4)
  expect_equal(W4[!u4], R4[!u4], tolerance = 1e-8)
  expect_equal(solve(W4)[u4], rep(0, 4), tolerance = 1e-7)
  # identified entries that admit no PD completion -> unknown entries set to 0
  Rb <- matrix(c(1, 0.99, NA, 0.99, 1, -0.99, NA, -0.99, 1), 3)
  Rb[1, 3] <- Rb[3, 1] <- NA
  ub <- is.na(Rb)
  expect_equal(cc(Rb, ub)[1, 3], -0.9801, tolerance = 1e-6)   # r_13 = r_12 * r_23
  Rx <- diag(4)
  Rx[1, 2] <- Rx[2, 1] <- 0.99; Rx[2, 3] <- Rx[3, 2] <- 0.99
  Rx[1, 3] <- Rx[3, 1] <- -0.99                               # inconsistent triangle
  Rx[1, 4] <- Rx[4, 1] <- NA
  ux <- is.na(Rx)
  expect_identical(cc(Rx, ux)[1, 4], 0)
  # a variable with no identified partner is independent of all others
  Ri <- matrix(c(1, NA, NA, NA, 1, 0.6, NA, 0.6, 1), 3)
  Wi <- cc(Ri, is.na(Ri))
  expect_equal(Wi[1, ], c(1, 0, 0))
  expect_equal(Wi[2, 3], 0.6)
})

test_that("r_mvskewnorm() requires a positive-definite omega", {
  expect_error(r_mvskewnorm(5, omega = matrix(1, 2, 2), alpha = c(1, 2)), "positive definite")
  expect_error(r_mvskewnorm(5, omega = matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), 3),
                            alpha = c(1, 0, 0)), "positive definite")
  # PSD-but-singular sigma stays allowed for the normal and t families
  expect_silent(r_mvt(5, df = 4, sigma = matrix(1, 2, 2), seed = 1))
})


# ---- 0.2.5: fifth external review -----------------------------------------

test_that("t-copula n_pair uses the full data, not the Kendall subsample", {
  set.seed(123)
  n <- 5000L
  d <- data.frame(x = rnorm(n), y = rnorm(n), z = rnorm(n))
  fit <- fit_synthesizer(d, copula = "t")
  np <- copula_correlation(fit, what = "n_pair")
  expect_identical(unname(np), matrix(n, 3, 3))
  expect_identical(fit$pooled$copula$completion, "none")
})

test_that("t-copula identification and estimation do not depend on the subsample", {
  n <- 10000L
  z <- r_mvnorm(n, sigma = make_corr(3, rho = 0.7), seed = 1)
  d <- data.frame(x = z[, 1], y = z[, 2], z = z[, 3])
  d$x[-(1:3000)] <- NA
  d$y[11:7000] <- NA                        # x and y jointly observed in 10 rows only
  res <- vapply(1:6, function(s) {
    set.seed(s)
    f <- suppressWarnings(fit_synthesizer(d, copula = "t", missing = "drop"))
    c(np = copula_correlation(f, what = "n_pair")["x", "y"],
      unid = f$pooled$copula$unidentified,
      rxy = copula_correlation(f)["x", "y"])
  }, numeric(3))
  expect_true(all(res["np", ] == 10))
  expect_true(all(res["unid", ] == 0))
  # all 10 joint rows are used, so the estimate is the same for every seed
  expect_lt(diff(range(res["rxy", ])), 1e-6)
  # same pattern below the subsampling size
  small <- d[c(1:10, sample(11:n, 1500)), ]
  f <- suppressWarnings(fit_synthesizer(small, copula = "t", missing = "drop"))
  expect_identical(copula_correlation(f, what = "n_pair")["x", "y"], 10L)
})

test_that(".kendall_matrix() equals stats::cor(method = 'kendall') with missing values", {
  set.seed(4)
  U <- matrix(runif(600 * 5), 600, 5)
  U[, 2] <- U[, 1] + 0.3 * runif(600)
  U[runif(length(U)) < 0.3] <- NA
  U[1:595, 5] <- NA                                # very sparse column
  a <- SynthesizerPlus:::.kendall_matrix(U)
  b <- suppressWarnings(cor(U, method = "kendall", use = "pairwise.complete.obs"))
  expect_identical(unname(is.na(a)), unname(is.na(b)))
  expect_equal(unname(a), unname(b), tolerance = 1e-12)
  # tiny blocks give the same answer
  expect_equal(SynthesizerPlus:::.kendall_matrix(U, block_cells = 1), a, tolerance = 1e-12)
  # fewer than 2 joint rows: undefined
  V <- cbind(c(0.1, NA, NA), c(NA, 0.2, 0.3), c(0.5, 0.6, 0.7))
  expect_true(is.na(SynthesizerPlus:::.kendall_matrix(V)[1, 2]))
})

test_that("completion method and repair size are recorded and reported", {
  cc <- SynthesizerPlus:::.complete_corr
  R <- matrix(c(1, NA, 0.5, NA, 1, 0.4, 0.5, 0.4, 1), 3)
  expect_identical(cc(R, is.na(R))$method, "maxdet")
  Rx <- diag(4)
  Rx[1, 2] <- Rx[2, 1] <- 0.99; Rx[2, 3] <- Rx[3, 2] <- 0.99
  Rx[1, 3] <- Rx[3, 1] <- -0.99
  Rx[1, 4] <- Rx[4, 1] <- NA
  expect_identical(cc(Rx, is.na(Rx))$method, "zero_fallback")

  d <- split_modules()
  fit <- suppressWarnings(fit_synthesizer(d, missing = "drop"))
  expect_identical(fit$pooled$copula$completion, "maxdet")
  expect_gte(fit$pooled$copula$pd_adjustment, 0)
  expect_identical(fit_synthesizer(iris)$pooled$copula$completion, "none")
  expect_output(print(summary(fit)), "maximum-determinant completion")
  fit$pooled$copula$completion <- "zero_fallback"
  fit$pooled$copula$pd_adjustment <- 0.05
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("no consistent completion; set to 0", out)))
  expect_true(any(grepl("repair changed correlations by up to 0.05", out)))
  # the warning mentions the later repair step
  expect_warning(fit_synthesizer(d, missing = "drop"), "before any positive-definiteness repair")
})

test_that("under-represented pairs are re-estimated from at most tau_rows joint rows", {
  set.seed(5)
  n <- 3000L
  z <- r_mvnorm(n, sigma = make_corr(3, rho = 0.6), seed = 5)
  U <- apply(z, 2, function(v) rank(v) / (n + 1))
  U[sample(n, 2700), 1] <- NA                     # column 1 sparse: 300 joint rows
  obs <- !is.na(U)
  np <- crossprod(obs + 0)
  # subsample of 100 rows leaves ~10 joint rows for column 1 (< min(300, 80)),
  # so those pairs are redone on 100 of their own 300 joint rows
  R <- SynthesizerPlus:::.pairwise_kendall_corr(U, obs, np, tau_rows = 100L, tau_min_rows = 80L)
  expect_true(all(is.finite(R)))
  expect_equal(R[1, 2], 0.6, tolerance = 0.35)
  expect_equal(R, t(R))
})

test_that("the fit warning describes a zero fallback when no completion exists", {
  d <- split_modules()
  testthat::local_mocked_bindings(.complete_corr = function(R, unid) {
    R[unid] <- 0
    list(R = R, method = "zero_fallback")
  })
  expect_warning(f <- fit_synthesizer(d, missing = "drop"), "could not be completed consistently")
  expect_identical(f$pooled$copula$completion, "zero_fallback")
})

test_that("mixed completion methods across strata are reported", {
  d <- split_modules(n = 400)
  d$g <- rep(c("a", "b"), times = 200)
  calls <- 0L
  testthat::local_mocked_bindings(.complete_corr = function(R, unid) {
    calls <<- calls + 1L
    R0 <- R
    R0[unid] <- 0
    list(R = R0, method = if (calls == 2L) "zero_fallback" else "maxdet")
  })
  expect_warning(fit_synthesizer(d, by = "g", missing = "drop", min_stratum_size = 10),
                 "set to 0 in 1 model")
})
