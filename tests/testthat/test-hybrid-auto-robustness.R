test_that("auto selection remains auditable for weak and high-dimensional cases", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("rpart")
  set.seed(628)
  n <- 160L
  df <- as.data.frame(matrix(rnorm(n * 10L), nrow = n, ncol = 10L))
  names(df) <- paste0("x", seq_len(10L))
  df$y <- .2 * df$x1 + rnorm(n)
  fit <- fit_synthesizer(df, nonlinear_target = "y", nonlinear_engine = "auto",
                         nonlinear_minbucket = 10L, nonlinear_cv_folds = 3L,
                         seed = 12)
  choice <- fit$pooled$nonlinear$selection
  expect_true(choice$selected %in% c("cart", "gam", "gam_ti"))
  expect_s3_class(choice$scores, "data.frame")
  expect_equal(nrow(choice$scores), 3L)
  # Interactions on >3 smooth predictors must fail explicitly, never silently.
  expect_true(!is.finite(choice$scores$mean_mse[choice$scores$candidate == "gam_ti"]))
  expect_equal(nrow(generate(fit, n = 50L, seed = 23)), 50L)
})
