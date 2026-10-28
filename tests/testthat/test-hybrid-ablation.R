test_that("hybrid ablation modes retain default and are reproducible", {
  skip_if_not_installed("rpart")
  set.seed(1001)
  x <- runif(400, -2, 2)
  d <- data.frame(x = x, z = rnorm(400),
                  y = sin(3 * x) + 0.5 * x + rnorm(400, sd = 0.3))
  f0 <- fit_synthesizer(d, nonlinear_target = "y", seed = 101,
                        nonlinear_minbucket = 12)
  f1 <- fit_synthesizer(d, nonlinear_target = "y", seed = 101,
                        nonlinear_minbucket = 12,
                        nonlinear_residual = "leaf",
                        nonlinear_calibration = "quantile")
  expect_identical(generate(f0, 300, seed = 202), generate(f1, 300, seed = 202))
  expect_identical(f1$pooled$nonlinear$residual_mode, "leaf")
  for (residual in c("leaf", "global", "none")) {
    for (calibration in c("quantile", "none")) {
      f <- fit_synthesizer(d, nonlinear_target = "y", seed = 101,
                           nonlinear_minbucket = 12,
                           nonlinear_residual = residual,
                           nonlinear_calibration = calibration)
      a <- generate(f, 300, seed = 202)
      b <- generate(f, 300, seed = 202)
      expect_identical(a, b)
      expect_true(all(is.finite(a$y)))
      expect_equal(nrow(a), 300L)
      expect_identical(names(a), names(d))
    }
  }
})

test_that("conditional mean only ablation has no residual jitter", {
  skip_if_not_installed("rpart")
  set.seed(803)
  d <- data.frame(x = runif(300), y = rnorm(300))
  f <- fit_synthesizer(d, nonlinear_target = "y", seed = 14,
                       nonlinear_minbucket = 15,
                       nonlinear_residual = "none",
                       nonlinear_calibration = "none")
  syn <- generate(f, n = 300, seed = 15)
  # In this variant each response is exactly a CART prediction group mean.
  expect_true(length(unique(syn$y)) <= length(f$pooled$nonlinear$centers))
  g <- fit_synthesizer(d, nonlinear_target = "y", seed = 14,
                       nonlinear_minbucket = 15,
                       nonlinear_residual = "global",
                       nonlinear_calibration = "none")
  expect_true(length(unique(generate(g, n = 300, seed = 15)$y)) >
                length(unique(syn$y)))
})

test_that("invalid ablation settings give informative errors", {
  skip_if_not_installed("rpart")
  set.seed(54)
  d <- data.frame(x = rnorm(100), y = rnorm(100))
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_residual = "local"), "arg")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_calibration = "linear"), "arg")
})
