test_that("GAM hybrid is opt-in, reproducible, and preserves CART defaults", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("rpart")
  set.seed(4101)
  x <- runif(450, -2, 2)
  d <- data.frame(x = x, z = rnorm(450),
                  group = factor(sample(c("A", "B"), 450, TRUE)))
  d$y <- sin(2.2 * d$x) + 0.8 * d$z + 0.4 * (d$group == "B") +
    rnorm(450, sd = 0.2)
  a <- fit_synthesizer(d, nonlinear_target = "y", seed = 4102)
  b <- fit_synthesizer(d, nonlinear_target = "y",
                       nonlinear_engine = "cart", seed = 4102)
  expect_identical(generate(a, 200, seed = 4103),
                   generate(b, 200, seed = 4103))
  f <- fit_synthesizer(d, nonlinear_target = "y",
                       nonlinear_engine = "gam", seed = 4102)
  expect_identical(f$pooled$nonlinear$engine, "gam")
  one <- generate(f, 500, seed = 4104)
  expect_identical(one, generate(f, 500, seed = 4104))
  expect_named(one, names(d))
  expect_equal(nrow(one), 500L)
  expect_true(all(is.finite(one$y)))
  expect_s3_class(one$group, "factor")
  expect_true(abs(stats::median(one$y) - stats::median(d$y)) < 0.2)
})

test_that("GAM hybrid handles missingness, factors, and optional interactions", {
  skip_if_not_installed("mgcv")
  set.seed(4201)
  n <- 500L
  d <- data.frame(x = runif(n, -2, 2), z = rnorm(n),
                  group = factor(sample(c("A", "B", "C"), n, TRUE)))
  d$y <- d$x * d$z + sin(2 * d$x) + 0.5 * (d$group == "C") +
    rnorm(n, sd = 0.3)
  d$x[sample.int(n, 24)] <- NA_real_
  d$y[sample.int(n, 35)] <- NA_real_
  for (residual in c("leaf", "global", "none")) {
    f <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "gam",
                         nonlinear_residual = residual,
                         nonlinear_gam_interactions = TRUE, seed = 4202)
    s <- generate(f, 500, seed = 4203)
    expect_true(all(is.finite(s$y[!is.na(s$y)])))
    expect_true(anyNA(s$y))
    expect_true(anyNA(s$x))
    expect_identical(levels(s$group), levels(d$group))
    expect_identical(s, generate(f, 500, seed = 4203))
  }
})

test_that("GAM option validation is specific and baseline remains unchanged", {
  skip_if_not_installed("mgcv")
  set.seed(4301)
  d <- data.frame(x = runif(120), y = rnorm(120))
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_engine = "unknown"), "arg")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_engine = "gam", nonlinear_gam_k = 2), "nonlinear_gam_k")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_engine = "gam", nonlinear_local_k = 0), "nonlinear_local_k")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
                nonlinear_engine = "gam", nonlinear_gam_interactions = NA),
               "nonlinear_gam_interactions")
  fa <- fit_synthesizer(d, seed = 4302)
  fb <- fit_synthesizer(d, nonlinear_target = NULL,
                        nonlinear_engine = "gam", seed = 4302)
  expect_identical(generate(fa, 100, seed = 4303),
                   generate(fb, 100, seed = 4303))
})
