test_that("automatic engine selection uses reproducible training-only CV", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("mgcv")
  set.seed(5101)
  n <- 260L
  x <- runif(n, -2, 2)
  z <- rnorm(n)
  d <- data.frame(x = x, z = z,
                  group = factor(sample(c("A", "B", "C"), n, TRUE)),
                  y = sin(2 * x) + 0.7 * x * z + 0.3 * rnorm(n))
  f1 <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto",
                        nonlinear_cv_folds = 3L, nonlinear_cv_rule = "one_se",
                        nonlinear_minbucket = 12L, seed = 5102)
  f2 <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto",
                        nonlinear_cv_folds = 3L, nonlinear_cv_rule = "one_se",
                        nonlinear_minbucket = 12L, seed = 5102)
  sel <- f1$pooled$nonlinear$selection
  expect_true(sel$selected %in% c("cart", "gam", "gam_ti"))
  expect_identical(sel$rule, "one_se")
  expect_equal(sel$folds, 3L)
  expect_named(sel$scores, c("candidate", "engine", "interactions",
                             "complexity", "n_valid", "mean_score", "se_score",
                             "mean_mse", "se_mse", "issue"))
  expect_identical(sel$score, "crps")
  expect_equal(nrow(sel$fold_mse), 3L)
  expect_equal(nrow(sel$fold_scores), 3L)
  expect_identical(sel$selected, f2$pooled$nonlinear$selection$selected)
  expect_equal(sel$fold_mse, f2$pooled$nonlinear$selection$fold_mse)
  if (identical(sel$reason, "cv")) {
    valid <- sel$scores[is.finite(sel$scores$mean_score), , drop = FALSE]
    best <- valid[which.min(valid$mean_score), , drop = FALSE]
    eligible <- valid[valid$mean_score <= best$mean_score + best$se_score + 1e-12,
                      , drop = FALSE]
    expect_equal(sel$selected,
                 eligible$candidate[which.min(eligible$complexity)])
  }
  # Auto must refit the selected engine on all training rows, with no
  # hidden transformation or change to the copula-fit pathway.
  manual <- fit_synthesizer(d, nonlinear_target = "y",
    nonlinear_engine = sel$engine,
    nonlinear_gam_interactions = sel$interactions,
    nonlinear_minbucket = 12L, seed = 5102)
  s1 <- generate(f1, 250L, seed = 5103)
  expect_identical(s1, generate(manual, 250L, seed = 5103))
  expect_identical(s1, generate(f1, 250L, seed = 5103))
  expect_identical(s1, generate(f2, 250L, seed = 5103))
  expect_equal(nrow(s1), 250L)
  expect_true(all(is.finite(s1$y)))
})

test_that("selector supports the minimum-score rule and missing targets", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("mgcv")
  set.seed(5201)
  n <- 210L
  d <- data.frame(x = runif(n, -2, 2), z = rnorm(n),
                  category = factor(sample(c("A", "B"), n, TRUE)))
  d$y <- 0.6 * d$x + 0.5 * d$z + rnorm(n, sd = .3)
  d$y[sample.int(n, 10L)] <- NA_real_
  f <- fit_synthesizer(d, nonlinear_target = "y", nonlinear_engine = "auto",
                       nonlinear_cv_rule = "min", nonlinear_cv_folds = 3L,
                       nonlinear_minbucket = 10L, seed = 5202)
  selection <- f$pooled$nonlinear$selection
  expect_identical(selection$rule, "min")
  valid <- selection$scores[is.finite(selection$scores$mean_score), , drop = FALSE]
  expect_true(nrow(valid) >= 1L)
  expect_equal(selection$selected, valid$candidate[which.min(valid$mean_score)])
  generated <- generate(f, 200, seed = 5203)
  expect_true(anyNA(generated$y))
  expect_true(all(is.finite(generated$y[!is.na(generated$y)])))
})

test_that("automatic selector falls back to CART in small strata", {
  skip_if_not_installed("rpart")
  set.seed(5301)
  n <- 165L
  d <- data.frame(grp = factor(c(rep("small", 24L), rep("big", 141L))),
                  x = runif(n), y = rnorm(n))
  f <- fit_synthesizer(d, by = "grp", nonlinear_target = "y",
                       nonlinear_engine = "auto", nonlinear_cv_folds = 5L,
                       nonlinear_minbucket = 10L, seed = 5302)
  choices <- lapply(f$strata$models, function(z) z$nonlinear$selection)
  expect_true(any(vapply(choices, function(z)
    identical(z$reason, "insufficient_observed_responses_for_cv"), logical(1))))
  s <- generate(f, 200L, seed = 5303)
  expect_equal(nrow(s), 200L)
  expect_true(all(is.finite(s$y)))
})

test_that("CV argument validation and unchanged defaults", {
  set.seed(5401)
  d <- data.frame(x = runif(140), y = rnorm(140))
  a <- fit_synthesizer(d, seed = 5402)
  b <- fit_synthesizer(d, nonlinear_target = NULL,
                       nonlinear_engine = "auto", seed = 5402)
  expect_identical(generate(a, 50, seed = 5403),
                   generate(b, 50, seed = 5403))
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
               nonlinear_engine = "auto", nonlinear_cv_folds = 1L),
               "nonlinear_cv_folds")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
               nonlinear_engine = "auto", nonlinear_cv_folds = 11L),
               "nonlinear_cv_folds")
  expect_error(fit_synthesizer(d, nonlinear_target = "y",
               nonlinear_engine = "auto", nonlinear_cv_rule = "other"),
               "arg")
})
