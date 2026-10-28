test_that("CART hybrid handles a category seen only with a missing target", {
  skip_if_not_installed("rpart")
  set.seed(7101)
  n <- 400L
  x <- rnorm(n)
  g <- sample(c("a", "b"), n, TRUE)
  y <- sin(2 * x) + (g == "a") + rnorm(n, sd = 0.3)
  g[1:15] <- "c"
  y[1:15] <- NA
  d <- data.frame(x, g, y)
  for (engine in c("cart", if (requireNamespace("mgcv", quietly = TRUE)) "gam")) {
    f <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = engine, seed = 1)
    s <- generate(f, 2000, seed = 2)
    expect_true(any(s$g == "c"))
    expect_true(all(is.finite(s$y[!is.na(s$y)])))
  }
})

test_that("hybrid calibration does not depend on the batch size", {
  skip_if_not_installed("rpart")
  set.seed(7201)
  n <- 500L
  x <- runif(n, -2, 2)
  d <- data.frame(x = x, y = x^2 + rnorm(n, sd = 0.5))
  f <- fit_synthesizer(d, nonlinear_target = "y", seed = 1)
  one <- vapply(1:60, function(s) generate(f, 1, seed = s)$y, numeric(1))
  expect_gt(length(unique(one)), 30L)            # not always the median
  small <- do.call(rbind, lapply(1:200, function(s) generate(f, 5, seed = s)))
  big <- generate(f, 1000, seed = 9)
  expect_lt(abs(sd(small$y) - sd(d$y)), 0.15 * sd(d$y))
  expect_gt(cor(small$x^2, small$y), 0.8)
  expect_true(all(small$y >= min(d$y) & small$y <= max(d$y)))
  # generation in chunks equals the same total distribution
  expect_lt(suppressWarnings(stats::ks.test(small$y, big$y)$statistic), 0.08)
})

test_that("strata without a conditional model are reported", {
  skip_if_not_installed("rpart")
  set.seed(7301)
  n <- 600L
  k <- sample(c("A", "B", "C"), n, TRUE, prob = c(.45, .45, .1))
  x <- rnorm(n)
  d <- data.frame(k, x, y = ifelse(k == "C", x^2, sin(x)) + rnorm(n, sd = .2))
  expect_warning(
    f <- fit_synthesizer(d, by = "k", nonlinear_target = "y",
                         nonlinear_minbucket = 40, seed = 1),
    "could not be fitted in 1 of 3")
  expect_identical(f$strata$conditional_fallback$k, "C")
  expect_output(print(f), "fallback")
  f2 <- suppressWarnings(fit_synthesizer(d, by = "k", nonlinear_target = "y", seed = 1))
  expect_equal(nrow(f2$strata$conditional_fallback), 0L)
})

test_that("CRPS selection rewards the conditional spread", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("mgcv")
  set.seed(7401)
  n <- 1500L
  x <- runif(n, -2, 2)
  w <- rnorm(n)
  # constant conditional mean in x; the spread jumps at x = 0
  d <- data.frame(x, w, y = 0.3 * w + rnorm(n, sd = 0.1 + 0.6 * (x > 0)))
  f <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto", seed = 1)
  sel <- f$pooled$nonlinear$selection
  expect_identical(sel$score, "crps")
  expect_identical(sel$selected, "cart")
  fm <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto",
                        nonlinear_cv_score = "mse", seed = 1)
  expect_identical(fm$pooled$nonlinear$selection$score, "mse")
  expect_equal(fm$pooled$nonlinear$selection$scores$mean_score,
               fm$pooled$nonlinear$selection$scores$mean_mse)
  expect_error(fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto",
                               nonlinear_cv_score = "rmse"), "should be one of")
})

test_that("nominal model captures non-monotone dependence on numeric predictors", {
  skip_if_not_installed("nnet")
  set.seed(7501)
  n <- 3000L
  x <- rnorm(n)
  z <- rnorm(n)
  eta <- cbind(0, 1.5 * x, -1.5 * x, 2 * abs(x) - 1.5)
  p <- exp(eta) / rowSums(exp(eta))
  y <- apply(p, 1L, function(pp) sample(c("A", "B", "C", "D"), 1L, prob = pp))
  d <- data.frame(x, z, y = factor(y))
  tvd <- function(D) {
    br <- c(-Inf, -1, 0, 1, Inf)
    r <- prop.table(table(cut(d$x, br), d$y), 1L)
    s <- prop.table(table(cut(D$x, br), D$y), 1L)
    mean(rowSums(abs(r - s)) / 2)
  }
  spl <- generate(fit_synthesizer(d, nominal_target = "y", seed = 1), 20000, seed = 2)
  lin <- generate(fit_synthesizer(d, nominal_target = "y", nominal_spline_df = 0, seed = 1),
                  20000, seed = 2)
  expect_lt(tvd(spl), 0.06)
  expect_lt(tvd(spl), tvd(lin))
  f <- fit_synthesizer(d, nominal_target = "y", seed = 1)
  expect_true(length(f$pooled$nominal$splines) >= 2L)   # x and z
  expect_identical(f$pooled$nominal$spline_df, 3L)       # AIC keeps splines
  expect_named(f$pooled$nominal$selection, c("spline_df", "converged", "aic"))
  expect_error(fit_synthesizer(d, nominal_target = "y", nominal_spline_df = 1),
               "nominal_spline_df")
})

test_that("nominal spline model converges on well-separated classes", {
  skip_if_not_installed("nnet")
  f <- fit_synthesizer(iris, nominal_target = "Species", seed = 42)
  expect_identical(f$pooled$nominal$model$convergence, 0L)
  expect_identical(f$pooled$nominal$spline_df, 0L)       # AIC prefers linear
  s <- generate(f, 300, seed = 43)
  expect_true(all(levels(iris$Species) %in% s$Species))
  # the classes stay separated by petal length
  expect_gt(mean(s$Petal.Length[s$Species == "setosa"] < 2.5), 0.9)
})
