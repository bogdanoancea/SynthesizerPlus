test_that("fit_synthesizer() returns a well-formed object", {
  fit <- fit_synthesizer(iris)
  expect_s3_class(fit, "sp_synthesizer")
  expect_equal(fit$n_train, 150L)
  expect_equal(fit$model_cols, names(iris))
  expect_output(print(fit), "sp_synthesizer")
  s <- summary(fit)
  expect_s3_class(s, "summary.sp_synthesizer")
  expect_equal(nrow(s$columns), 5L)
  expect_output(print(s), "Copula")
})

test_that("generate() preserves names, classes, levels and size", {
  d <- mixed_data()
  fit <- fit_synthesizer(d)
  syn <- generate(fit, n = 250, seed = 1)
  expect_s3_class(syn, "data.frame")
  expect_equal(nrow(syn), 250L)
  expect_named(syn, names(d))
  expect_equal(lapply(syn, class), lapply(d, class))
  expect_equal(levels(syn$fac), levels(d$fac))
  expect_equal(levels(syn$ord), levels(d$ord))
  expect_true(is.ordered(syn$ord))
  expect_equal(attr(syn$time, "tzone"), "UTC")
  expect_true(all(syn$chr %in% d$chr))
  expect_false("unused" %in% syn$fac)
})

test_that("default n equals the training size and n = 0 works", {
  fit <- fit_synthesizer(mtcars)
  expect_equal(nrow(generate(fit)), nrow(mtcars))
  z <- generate(fit, n = 0)
  expect_equal(nrow(z), 0L)
  expect_named(z, names(mtcars))
})

test_that("continuous and discrete values stay within the observed support", {
  d <- mixed_data()
  syn <- generate(fit_synthesizer(d), n = 2000, seed = 2)
  expect_true(all(syn$num >= min(d$num) & syn$num <= max(d$num)))
  expect_true(all(syn$int >= min(d$int) & syn$int <= max(d$int)))
  expect_true(all(syn$cnt %in% d$cnt))
  expect_true(all(syn$date >= min(d$date) & syn$date <= max(d$date)))
  expect_true(all(syn$time >= min(d$time) & syn$time <= max(d$time)))
})

test_that("marginal distributions are reproduced", {
  set.seed(3)
  x <- data.frame(a = rgamma(2000, 2), b = sample(letters[1:4], 2000, TRUE,
                                                  prob = c(0.5, 0.3, 0.15, 0.05)))
  syn <- generate(fit_synthesizer(x), n = 20000, seed = 4)
  expect_lt(abs(mean(syn$a) - mean(x$a)), 0.05)
  expect_lt(abs(sd(syn$a) - sd(x$a)), 0.05)
  p_real <- prop.table(table(x$b))
  p_syn <- prop.table(table(syn$b))
  expect_lt(max(abs(p_real - p_syn)), 0.02)
})

test_that("rare categories are not inflated", {
  x <- data.frame(f = factor(c(rep("a", 999), "b")), y = rnorm(1000))
  syn <- generate(fit_synthesizer(x), n = 10000, seed = 1)
  expect_lt(mean(syn$f == "b"), 0.005)
})

test_that("rank correlation between numeric variables is preserved", {
  S <- make_corr(3, "exchangeable", 0.7)
  x <- as.data.frame(r_mvnorm(3000, sigma = S, seed = 1))
  syn <- generate(fit_synthesizer(x), n = 10000, seed = 2)
  diff <- cor(syn, method = "spearman") - cor(x, method = "spearman")
  expect_lt(max(abs(diff)), 0.05)
})

test_that("categorical variables take part in the dependence structure", {
  set.seed(1)
  x <- rnorm(3000)
  d <- data.frame(x = x, g = factor(ifelse(x > 0, "hi", "lo")))
  syn <- generate(fit_synthesizer(d), seed = 1)
  concordance <- mean((syn$x > 0) == (syn$g == "hi"))
  expect_gt(concordance, 0.7)
})

test_that("dependence = 0 removes the dependence but keeps the marginals", {
  S <- make_corr(2, "exchangeable", 0.9)
  x <- as.data.frame(r_mvnorm(3000, sigma = S, seed = 1))
  fit <- fit_synthesizer(x)
  ind <- generate(fit, n = 10000, dependence = 0, seed = 1)
  half <- generate(fit, n = 10000, dependence = 0.5, seed = 1)
  full <- generate(fit, n = 10000, dependence = 1, seed = 1)
  r <- function(d) cor(d[[1]], d[[2]], method = "spearman")
  expect_lt(abs(r(ind)), 0.05)
  expect_true(r(ind) < r(half) && r(half) < r(full))
  expect_lt(abs(sd(half[[1]]) - sd(x[[1]])), 0.05)
  expect_error(generate(fit, dependence = 2), "dependence")
})

test_that("missing values are modelled", {
  aq <- airquality
  fit <- fit_synthesizer(aq)
  syn <- generate(fit, n = 20000, seed = 1)
  expect_lt(abs(mean(is.na(syn$Ozone)) - mean(is.na(aq$Ozone))), 0.02)
  expect_lt(abs(mean(is.na(syn$Solar.R)) - mean(is.na(aq$Solar.R))), 0.01)
  expect_false(anyNA(syn$Wind))

  syn_ind <- generate(fit_synthesizer(aq, missing = "independent"), n = 20000, seed = 1)
  expect_lt(abs(mean(is.na(syn_ind$Ozone)) - mean(is.na(aq$Ozone))), 0.02)

  syn_drop <- generate(fit_synthesizer(aq, missing = "drop"), seed = 1)
  expect_false(anyNA(syn_drop))
})

test_that("joint missing-data model reproduces association of missingness", {
  set.seed(2)
  n <- 4000
  z <- rnorm(n)
  y <- rnorm(n)
  y[z > 0.5] <- NA # y is missing for large z
  d <- data.frame(z = z, y = y)
  syn <- generate(fit_synthesizer(d), n = n, seed = 3)
  rate_hi <- mean(is.na(syn$y[syn$z > 1]))
  rate_lo <- mean(is.na(syn$y[syn$z < -1]))
  expect_gt(rate_hi, rate_lo + 0.4)
})

test_that("all-missing columns are reproduced as missing", {
  d <- data.frame(a = rnorm(20), b = NA_real_, c = factor(rep(NA, 20), levels = "x"))
  syn <- generate(fit_synthesizer(d), seed = 1)
  expect_true(all(is.na(syn$b)))
  expect_true(all(is.na(syn$c)))
  expect_s3_class(syn$c, "factor")
})

test_that("stratified models keep strata consistent", {
  fit <- fit_synthesizer(iris, by = "Species")
  syn <- generate(fit, n = 3000, seed = 1)
  m_syn <- tapply(syn$Petal.Length, syn$Species, mean)
  m_real <- tapply(iris$Petal.Length, iris$Species, mean)
  expect_lt(max(abs(m_syn - m_real)), 0.1)
  expect_true(all(syn$Petal.Length[syn$Species == "setosa"] <=
                    max(iris$Petal.Length[iris$Species == "setosa"])))
  prop <- generate(fit, n = 300, stratum_sizes = "proportional", seed = 1)
  expect_equal(as.vector(table(prop$Species)), c(100L, 100L, 100L))
})

test_that("small strata fall back to the pooled model", {
  d <- data.frame(g = c(rep("big", 100), rep("small", 5)), x = rnorm(105))
  fit <- fit_synthesizer(d, by = "g", min_stratum_size = 10)
  expect_null(fit$strata$models[["small"]])
  expect_output(print(summary(fit)), "pooled")
  syn <- generate(fit, n = 1000, seed = 1)
  expect_true(all(c("big", "small") %in% syn$g))
})

test_that("multi-column strata work", {
  fit <- fit_synthesizer(mtcars, by = c("am", "vs"), min_stratum_size = 5)
  syn <- generate(fit, seed = 1)
  combos <- unique(paste(syn$am, syn$vs))
  expect_true(all(combos %in% unique(paste(mtcars$am, mtcars$vs))))
})

test_that("id columns are regenerated", {
  d <- data.frame(id = sprintf("P%03d", 1:50), x = rnorm(50), stringsAsFactors = FALSE)
  syn <- generate(fit_synthesizer(d, id_cols = "id"), n = 20)
  expect_false(any(syn$id %in% d$id))
  expect_equal(anyDuplicated(syn$id), 0L)
  d2 <- data.frame(id = 1:50, x = rnorm(50))
  expect_equal(generate(fit_synthesizer(d2, id_cols = "id"), n = 5)$id, 1:5)
})

test_that("seeds give reproducible output and do not touch the global RNG", {
  fit <- fit_synthesizer(iris)
  expect_identical(generate(fit, seed = 9), generate(fit, seed = 9))
  set.seed(42)
  a <- runif(1)
  set.seed(42)
  generate(fit, seed = 1)
  b <- runif(1)
  expect_identical(a, b)
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  generate(fit, seed = 1)
  expect_false(exists(".Random.seed", envir = globalenv()))
  set.seed(1)
})

test_that("vectors, matrices and single columns are supported", {
  v <- synthesize(rnorm(50), n = 10, seed = 1)
  expect_type(v, "double")
  expect_length(v, 10L)
  f <- synthesize(factor(c("a", "b", "b")), n = 5, seed = 1)
  expect_s3_class(f, "factor")
  m <- synthesize(as.matrix(mtcars[, 1:3]), n = 4, seed = 1)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(4L, 3L))
  expect_equal(colnames(m), names(mtcars)[1:3])
  m2 <- synthesize(matrix(rnorm(40), 20), seed = 1)
  expect_equal(colnames(m2), c("V1", "V2"))
  expect_s3_class(synthesize(Sys.Date() + 0:9, n = 3), "Date")
})

test_that("t copula and independence copula are available", {
  ft <- fit_synthesizer(mtcars, copula = "t")
  expect_true(is.finite(ft$pooled$copula$df))
  expect_equal(nrow(generate(ft, n = 10)), 10L)
  fi <- fit_synthesizer(mtcars, copula = "independence")
  expect_equal(copula_correlation(fi), diag(11), ignore_attr = TRUE)
})

test_that("t copula estimates heavy joint tails", {
  set.seed(1)
  u <- r_copula(3000, "t", corr = make_corr(2, "exchangeable", 0.5), df = 3, seed = 1)
  fit <- fit_synthesizer(as.data.frame(qnorm(u)), copula = "t")
  expect_lte(fit$pooled$copula$df, 6)
})

test_that("copula_correlation() returns a valid correlation matrix", {
  fit <- fit_synthesizer(airquality, by = "Month")
  R <- copula_correlation(fit)
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(unname(diag(R)), rep(1, 5))
  expect_gt(nrow(copula_correlation(fit, include_missing = TRUE)), 5L)
  expect_true(is.matrix(copula_correlation(fit, stratum = "5")))
  expect_error(copula_correlation(fit, stratum = "13"), "Unknown stratum")
  expect_error(copula_correlation(1), "sp_synthesizer")
})

test_that("singular correlation matrices are repaired", {
  x <- data.frame(a = 1:50 + 0.1 * rnorm(50))
  x$b <- x$a
  x$c <- -x$a
  for (m in c("auto", "higham", "eigen")) {
    fit <- fit_synthesizer(x, pd_method = m)
    expect_true(min(eigen(fit$pooled$copula$R)$values) > 0)
    expect_equal(nrow(generate(fit, n = 10)), 10L)
  }
})

test_that("spline interpolation and knots work", {
  set.seed(1)
  x <- data.frame(v = rexp(5000))
  fs <- fit_synthesizer(x, interpolation = "spline")
  expect_true(all(generate(fs, n = 1000, seed = 1)$v >= 0))
  fk <- fit_synthesizer(x, knots = 100)
  expect_length(fk$pooled$marginals$v$quantiles, 100L)
  expect_equal(mean(generate(fk, n = 20000, seed = 1)$v), mean(x$v), tolerance = 0.05)
  expect_lt(object.size(fk), object.size(fit_synthesizer(x)))
})

test_that("discrete_threshold switches between discrete and continuous models", {
  d <- data.frame(k = rep(1:5, 20))
  expect_equal(fit_synthesizer(d)$pooled$marginals$k$type, "discrete")
  expect_equal(fit_synthesizer(d, discrete_threshold = 2)$pooled$marginals$k$type,
               "continuous")
})

test_that("input validation gives informative errors", {
  expect_error(fit_synthesizer(iris[1, ]), "at least 2 rows")
  expect_error(fit_synthesizer(iris, by = "nope"), "not found")
  expect_error(fit_synthesizer(iris, id_cols = "Species", by = "Species"), "both")
  expect_error(fit_synthesizer(iris, id_cols = names(iris)), "No columns")
  expect_error(fit_synthesizer(list(a = 1)), "Cannot fit")
  d <- data.frame(x = 1:3)
  d$l <- list(1, 2, 3)
  expect_error(fit_synthesizer(d), "Unsupported")
  fit <- fit_synthesizer(iris)
  expect_error(generate(fit, n = -1), "non-negative")
  expect_error(generate(fit, n = 1.5), "integer")
  expect_error(generate(fit, seed = "a"), "seed")
  expect_error(fit_synthesizer(iris, knots = 1), "knots")
})

test_that("synthesize() fits, generates and can write to a file", {
  s <- synthesize(iris, n = 30, seed = 1)
  expect_equal(nrow(s), 30L)
  expect_identical(s, synthesize(iris, n = 30, seed = 1))
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  res <- synthesize(iris, n = 20, seed = 1, file = f)
  expect_true(file.exists(f))
  expect_identical(readRDS(f), res)
  expect_equal(nrow(synthesize(airquality, by = "Month", seed = 1)), 153L)
})

test_that("simulate() returns multiple synthetic data sets", {
  fit <- fit_synthesizer(mtcars)
  sims <- simulate(fit, nsim = 3, seed = 1, n = 10)
  expect_length(sims, 3L)
  expect_true(all(vapply(sims, nrow, integer(1)) == 10L))
  expect_false(identical(sims[[1]], sims[[2]]))
  expect_identical(sims, simulate(fit, nsim = 3, seed = 1, n = 10))
})

test_that("print methods for marginals work", {
  m <- fit_synthesizer(airquality)$pooled$marginals$Ozone
  expect_output(print(m), "missing")
})

test_that("classed numeric columns such as difftime keep class and units", {
  d <- data.frame(dur = as.difftime(runif(100, 1, 60), units = "mins"), x = rnorm(100))
  fit <- fit_synthesizer(d)
  s <- generate(fit, n = 500, seed = 1)
  expect_s3_class(s$dur, "difftime")
  expect_equal(units(s$dur), "mins")
  expect_true(all(as.numeric(s$dur) >= 1 & as.numeric(s$dur) <= 60))
  sh <- generate(fit, n = 5000, seed = 1, perturb = perturbation(shift = c(dur = 1)))
  expect_s3_class(sh$dur, "difftime")
  expect_gt(mean(as.numeric(sh$dur)), mean(as.numeric(d$dur)) + 10)
  m <- marginal_metrics(d, s)
  expect_equal(m$type[1], "numeric")
})
