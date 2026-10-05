test_that("census_sim is available and well-formed", {
  expect_equal(dim(census_sim), c(5000L, 8L))
  expect_true(is.ordered(census_sim$education))
  expect_false(anyNA(census_sim))
})

test_that("releasing the real data is maximal disclosure", {
  keys <- c("age", "sex", "region", "household_size")
  r <- disclosure_risk(census_sim, census_sim, keys = keys, target = c("income", "education"))
  expect_s3_class(r, "sp_disclosure")
  expect_equal(r$replicated_uniques, 1)
  expect_equal(r$exact_copies, 1)
  expect_equal(r$attribute$cap_uniques, c(1, 1))
  expect_true(all(r$attribute$match_rate == 1))
  expect_gt(r$real_unique_rate, 0.2)
  expect_output(print(r), "replicated uniques")
})

test_that("synthetic data reduce attribute disclosure to about the baseline", {
  real <- census_sim[1:2000, ]
  syn <- generate(fit_synthesizer(real, id_cols = "person_id", by = "employment"), seed = 1)
  r <- disclosure_risk(real, syn, keys = c("age", "sex", "region", "household_size"),
                       target = c("income", "employment"))
  expect_lt(r$exact_copies, 0.01)
  expect_lt(r$attribute$cap_uniques[1], 0.3)
  expect_lt(r$attribute$cap[1], r$attribute$cap_baseline[1] + 0.1)
})

test_that("disclosure_risk() handles numeric tolerance, NAs and validation", {
  real <- data.frame(k = c(1, 1, 2, 3), y = c(100, 110, 200, NA), g = c("a", "a", "b", "c"))
  syn <- data.frame(k = c(1, 2, 2, 4), y = c(104, 300, 300, 1), g = c("a", "b", "b", "c"))
  r <- disclosure_risk(real, syn, keys = "k", target = c("y", "g"), tolerance = 0.05)
  # k = 1: prediction 104, correct for 100 (4%) and 110 (5.5% -> wrong)
  expect_equal(r$attribute$cap[1], 1 / 3)
  expect_equal(r$attribute$match_rate[2], 0.75)
  expect_equal(r$attribute$cap[2], 1)
  expect_equal(r$n_real_uniques, 2L)
  expect_equal(r$replicated_uniques, 0)
  expect_null(disclosure_risk(real, syn, keys = "k")$attribute)
  expect_error(disclosure_risk(real, syn), "keys")
  expect_error(disclosure_risk(real, syn, keys = "zz"), "not found")
  expect_error(disclosure_risk(real, syn, keys = "k", target = "k"), "both")
  expect_error(disclosure_risk(real, syn, keys = "k", tolerance = -1), "tolerance")
})

test_that("noise lowers the share of near-copies", {
  real <- census_sim[1:1500, -1]
  fit <- fit_synthesizer(real)
  s1 <- generate(fit, closeness = 1, seed = 1)
  s0 <- generate(fit, closeness = 0.25, seed = 1)
  d1 <- dcr(real, s1, seed = 1)
  d0 <- dcr(real, s0, seed = 1)
  expect_gt(d1$close_share, d0$close_share)
  expect_gt(d0$ratio, d1$ratio)
  expect_equal(dcr(real, real[1:200, ], seed = 1)$close_share, 1)
  expect_error(dcr(real, s1, close_quantile = 2), "close_quantile")
})

test_that("ignore removes identifiers from the exact-copy check", {
  real <- census_sim[1:300, ]
  copy <- real
  copy$person_id <- sprintf("X%05d", seq_len(nrow(copy)))   # new ids, copied records
  with_id <- disclosure_risk(real, copy, keys = c("age", "sex"))
  without <- disclosure_risk(real, copy, keys = c("age", "sex"), ignore = "person_id")
  expect_equal(with_id$exact_copies, 0)
  expect_equal(without$exact_copies, 1)
  expect_error(disclosure_risk(real, copy, keys = "age", ignore = "age"), "cannot contain")
  expect_error(disclosure_risk(real, copy, keys = "age", ignore = "zz"), "not found")
})

test_that("CAP is computed correctly on a hand-checked example", {
  real <- data.frame(k1 = c("a", "a", "b", "b", "c"), k2 = c(1, 1, 1, 2, 1),
                     y = c("u", "v", "u", "u", "v"), stringsAsFactors = FALSE)
  syn <- data.frame(k1 = c("a", "a", "a", "b", "d"), k2 = c(1, 1, 1, 2, 1),
                    y = c("v", "v", "u", "u", "u"), stringsAsFactors = FALSE)
  r <- disclosure_risk(real, syn, keys = c("k1", "k2"), target = "y")
  # real keys: a1 (x2), b1, b2, c1; synthetic keys: a1 (x3, mode "v"), b2 ("u"), d1
  expect_equal(r$n_real_uniques, 3L)                  # b1, b2, c1
  expect_equal(r$replicated_uniques, 1 / 3)           # b2 occurs once in synthetic
  expect_equal(r$attribute$match_rate, 3 / 5)         # a1, a1, b2 matched
  expect_equal(r$attribute$cap, 2 / 3)                # a1 -> v: wrong for u, right for v; b2 -> u right
  expect_equal(r$attribute$cap_uniques, 1)            # only b2 matched among uniques
  expect_equal(r$attribute$cap_baseline, 3 / 5)       # overall mode "u"
})

test_that("disclosure risk decreases relative to releasing the real data", {
  real <- census_sim[1:1500, ]
  fit <- fit_synthesizer(real, id_cols = "person_id", by = "employment")
  syn <- generate(fit, seed = 4)
  keys <- c("age", "sex", "region", "household_size")
  r_real <- disclosure_risk(real, real, keys, c("income", "education"), ignore = "person_id")
  r_syn <- disclosure_risk(real, syn, keys, c("income", "education"), ignore = "person_id")
  expect_true(all(r_syn$attribute$cap < r_real$attribute$cap))
  expect_lt(r_syn$exact_copies, r_real$exact_copies)
  expect_lt(r_syn$replicated_uniques, r_real$replicated_uniques)
})

test_that("close_share responds to near-copies", {
  real <- census_sim[1:600, -1]
  noisy_copy <- real
  noisy_copy$income <- noisy_copy$income * 1.001
  expect_gt(dcr(real, noisy_copy, seed = 1)$close_share, 0.9)
  far <- real
  far$income <- far$income * 10
  far$age <- 90L
  expect_lt(dcr(real, far, seed = 1)$close_share, 0.05)
})

test_that("calibrate_closeness() works with a disclosure-risk metric", {
  real <- census_sim[1:600, ]
  fit <- fit_synthesizer(real, id_cols = "person_id")
  cap_income <- function(r, s) {
    disclosure_risk(r, s, keys = c("age", "sex", "region"), target = "income",
                    ignore = "person_id")$attribute$cap
  }
  cal <- suppressWarnings(calibrate_closeness(fit, real, target = 0.05, metric = cap_income,
                                              grid = c(0, 0.5, 1), refine = 1))
  expect_s3_class(cal, "sp_calibration")
  expect_true(cal$closeness >= 0 && cal$closeness <= 1)
  expect_true(all(is.finite(cal$path$value)))
})
