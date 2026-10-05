test_that("marginal_metrics() is zero for identical data and positive otherwise", {
  m0 <- marginal_metrics(iris, iris)
  expect_equal(nrow(m0), 5L)
  expect_equal(m0$distance, rep(0, 5))
  expect_equal(m0$wasserstein[1:4], rep(0, 4))
  expect_equal(m0$statistic, c(rep("KS", 4), "TVD"))
  shifted <- iris
  shifted$Sepal.Length <- shifted$Sepal.Length + 1
  shifted$Species <- factor("setosa", levels = levels(iris$Species))
  m1 <- marginal_metrics(iris, shifted)
  expect_gt(m1$distance[1], 0.4)
  expect_equal(m1$distance[5], 2 / 3, tolerance = 1e-8)
  expect_equal(nrow(marginal_metrics(iris, iris, vars = "Species")), 1L)
  aq <- marginal_metrics(airquality, airquality)
  expect_equal(aq$missing_real[1], mean(is.na(airquality$Ozone)))
})

test_that("association_matrix() handles mixed types", {
  A <- association_matrix(iris)
  expect_equal(dim(A), c(5L, 5L))
  expect_equal(unname(diag(A)), rep(1, 5))
  expect_true(isSymmetric(A))
  expect_gt(A["Petal.Length", "Species"], 0.9)
  d <- data.frame(a = factor(rep(c("x", "y"), 50)), b = factor(rep(c("x", "y"), 50)),
                  c = factor(sample(c("u", "v"), 100, TRUE)))
  B <- association_matrix(d)
  expect_gt(B["a", "b"], 0.95)
  expect_lt(B["a", "c"], 0.3)
  expect_equal(association_matrix(data.frame(z = 1:3)), matrix(1, dimnames = list("z", "z")))
  # constant and mostly missing columns
  e <- data.frame(k = rep(1, 10), g = factor(rep("a", 10)), n = c(1, 2, rep(NA, 8)),
                  x = 1:10)
  E <- association_matrix(e)
  expect_true(is.na(E["n", "x"]))
  expect_equal(E["g", "x"], 0)
})

test_that("pmse() detects differences", {
  set.seed(1)
  real <- data.frame(x = rnorm(1000), g = sample(c("a", "b"), 1000, TRUE))
  same <- data.frame(x = rnorm(1000), g = sample(c("a", "b"), 1000, TRUE))
  diff <- data.frame(x = rnorm(1000, 1), g = sample(c("a", "b"), 1000, TRUE))
  p0 <- pmse(real, same, seed = 1)
  p1 <- pmse(real, diff, seed = 1)
  expect_named(p0, c("pmse", "pmse_null", "pmse_ratio", "n_params"))
  q1 <- pmse(real, diff, seed = 1, model = "quadratic")
  expect_gt(q1$n_params, p1$n_params)
  expect_lt(p0$pmse_ratio, 5)
  expect_gt(p1$pmse_ratio, 50)
  expect_gt(p1$pmse, p0$pmse)
})

test_that("discriminator_auc() is near 0.5 for exchangeable data", {
  set.seed(2)
  real <- data.frame(x = rnorm(800), y = runif(800))
  same <- data.frame(x = rnorm(800), y = runif(800))
  diff <- data.frame(x = rnorm(800, 2), y = runif(800))
  expect_lt(abs(discriminator_auc(real, same, seed = 1) - 0.5), 0.06)
  expect_gt(discriminator_auc(real, diff, seed = 1), 0.85)
  expect_error(discriminator_auc(real, same, folds = 1), "at least 2")
})

test_that("classifier features handle missing values and many levels", {
  d <- data.frame(x = c(NA, rnorm(99)), g = sample(letters, 100, TRUE),
                  k = 1, stringsAsFactors = FALSE)
  X <- SynthesizerPlus:::.prep_features(d, max_levels = 5)
  expect_true("x..NA" %in% colnames(X))
  expect_equal(sum(grepl("^g==", colnames(X))), 4L)
  expect_false("k" %in% colnames(X))
  expect_false(anyNA(X))
  expect_true(is.finite(pmse(d, d[sample(100), ])$pmse))
})

test_that("ci_overlap() is 1 for identical data", {
  res <- ci_overlap(mtcars, mtcars, mpg ~ wt + hp)
  expect_equal(res$overlap, rep(1, 3))
  expect_equal(res$std_diff, rep(0, 3))
  syn <- synthesize(mtcars, seed = 1)
  res2 <- ci_overlap(mtcars, syn, am ~ wt, family = binomial())
  expect_equal(nrow(res2), 2L)
  expect_true(all(res2$overlap <= 1))
})

test_that("dcr() flags copies of real records", {
  copy <- iris[sample(150, 100), ]
  d <- dcr(iris, copy, seed = 1)
  expect_equal(d$exact_match_rate, 1)
  expect_length(d$synthetic, 100L)
  expect_length(d$real, 150L)
  expect_equal(dim(d$quantiles), c(2L, 4L))
  far <- data.frame(Sepal.Length = 100, Sepal.Width = 100, Petal.Length = 100,
                    Petal.Width = 100, Species = factor("x"))
  far <- far[rep(1, 10), ]
  expect_gt(dcr(iris, far)$ratio, 10)
  # missing values in Gower distance
  aq <- dcr(airquality, airquality, seed = 1)
  expect_equal(aq$exact_match_rate, 1)
})

test_that("compare_synthetic() bundles all metrics", {
  syn <- synthesize(iris, seed = 1)
  cmp <- compare_synthetic(iris, syn, seed = 1)
  expect_s3_class(cmp, "sp_comparison")
  expect_named(cmp, c("vars", "n_real", "n_synthetic", "real", "synthetic",
                      "marginal", "association", "utility", "privacy"))
  expect_output(print(cmp), "Utility")
  s <- summary(cmp)
  expect_equal(nrow(s), 7L)
  expect_true(all(is.finite(s$value)))
  expect_gt(cmp$utility$auc, 0.3)
  expect_lt(cmp$utility$auc, 0.7)
  only <- compare_synthetic(iris, syn, metrics = "marginal")
  expect_null(only$utility)
  expect_true(is.na(summary(only)$value[4]))
  expect_identical(compare_synthetic(iris, syn, seed = 3)$utility,
                   compare_synthetic(iris, syn, seed = 3)$utility)
})

test_that("compare_synthetic() accepts vectors, matrices and subsets", {
  cmp <- compare_synthetic(rnorm(100), rnorm(100), metrics = c("marginal", "utility"))
  expect_equal(cmp$marginal$variable, "x")
  m <- as.matrix(mtcars)
  expect_s3_class(compare_synthetic(m, m, metrics = "marginal"), "sp_comparison")
  sub <- compare_synthetic(iris, iris, vars = c("Sepal.Length", "Species"),
                           metrics = "association")
  expect_equal(dim(sub$association$real), c(2L, 2L))
})

test_that("metric inputs are validated", {
  expect_error(compare_synthetic(iris, mtcars), "no columns in common")
  expect_error(marginal_metrics(iris[1, ], iris), "at least 2 rows")
  expect_error(marginal_metrics(iris, iris, vars = "zz"), "not found")
  expect_error(pmse(mean, iris), "data frames")
})
