test_that("perturbation() validates and prints", {
  p <- perturbation(noise = 0.5, swap = c(Species = 0.2), shift = c(a = 1, b = -1))
  expect_s3_class(p, "sp_perturbation")
  expect_output(print(p), "Species = 0.2")
  expect_error(perturbation(noise = -1), "noise")
  expect_error(perturbation(swap = 2), "swap")
  expect_error(perturbation(scale = 0), "scale")
  expect_error(perturbation(temperature = -1), "temperature")
  expect_error(perturbation(dependence = 3), "dependence")
  expect_error(perturbation(shift = c(1, 2)), "named")
  expect_error(perturbation(shift = NA_real_), "shift")
})

test_that("closeness = 1 and an empty perturbation reproduce the default output", {
  fit <- fit_synthesizer(iris)
  base <- generate(fit, seed = 1)
  expect_identical(generate(fit, closeness = 1, seed = 1), base)
  expect_identical(generate(fit, perturb = perturbation(), seed = 1), base)
  expect_error(generate(fit, closeness = 1.5), "closeness")
  expect_error(generate(fit, perturb = list(noise = 1)), "perturbation")
  expect_error(generate(fit, perturb = perturbation(shift = c(zz = 1))), "not in the data")
})

test_that("lower closeness moves synthetic data further from the real data", {
  fit <- fit_synthesizer(iris)
  d <- sapply(c(1, 0.5, 0), function(cl) {
    dcr(iris, generate(fit, closeness = cl, seed = 1), seed = 1)$ratio
  })
  expect_true(all(diff(d) > 0))
  a <- sapply(c(1, 0), function(cl) {
    discriminator_auc(iris, generate(fit, closeness = cl, seed = 1), seed = 1)
  })
  expect_gt(a[2], a[1] + 0.1)
  # closeness 0: independent, but classes and levels preserved
  s0 <- generate(fit, closeness = 0, seed = 1, n = 3000)
  expect_equal(lapply(s0, class), lapply(iris, class))
  expect_lt(abs(cor(s0$Petal.Length, s0$Petal.Width)), 0.1)
})

test_that("noise preserves mean and variance but creates new values", {
  set.seed(1)
  d <- data.frame(x = rgamma(500, 2), k = sample(1:3, 500, TRUE))
  fit <- fit_synthesizer(d)
  s <- generate(fit, n = 30000, perturb = perturbation(noise = 1), seed = 1)
  expect_equal(mean(s$x), mean(d$x), tolerance = 0.03)
  expect_equal(sd(s$x), sd(d$x), tolerance = 0.05)
  expect_true(all(s$x >= 0))            # sign bound
  near <- generate(fit, n = 500, seed = 2)
  far <- generate(fit, n = 500, perturb = perturbation(noise = 1), seed = 2)
  expect_gt(dcr(d, far, seed = 1)$ratio, dcr(d, near, seed = 1)$ratio)
  o <- generate(fit, n = 5000, perturb = perturbation(noise = 1, bounds = "observed"), seed = 1)
  expect_true(all(o$x >= min(d$x) & o$x <= max(d$x)))
  expect_true(all(s$k == round(s$k)))    # discrete variable stays whole-numbered
  n <- generate(fit, n = 5000, perturb = perturbation(noise = 3, bounds = "none"), seed = 1)
  expect_lt(min(n$x), 0)
})

test_that("integer and date variables stay integer-valued under noise", {
  d <- data.frame(i = 1:200, dt = as.Date("2020-01-01") + 0:199)
  s <- generate(fit_synthesizer(d), perturb = perturbation(noise = 0.7), seed = 1)
  expect_type(s$i, "integer")
  expect_s3_class(s$dt, "Date")
  expect_true(all(as.numeric(s$dt) == round(as.numeric(s$dt))))
})

test_that("swap keeps marginals and weakens association", {
  set.seed(2)
  x <- rnorm(3000)
  d <- data.frame(x = x, g = factor(ifelse(x > 0, "hi", "lo")))
  fit <- fit_synthesizer(d)
  s0 <- generate(fit, n = 20000, seed = 1)
  s1 <- generate(fit, n = 20000, perturb = perturbation(swap = 0.8), seed = 1)
  conc <- function(s) mean((s$x > 0) == (s$g == "hi"))
  expect_lt(conc(s1), conc(s0) - 0.15)
  expect_equal(mean(s1$g == "hi"), mean(d$g == "hi"), tolerance = 0.03)
})

test_that("distribution-level perturbations act as specified", {
  fit <- fit_synthesizer(iris)
  s <- generate(fit, n = 20000, seed = 1,
                perturb = perturbation(shift = c(Sepal.Length = 2),
                                       scale = c(Sepal.Width = 2)))
  expect_equal(mean(s$Sepal.Length) - mean(iris$Sepal.Length), 2 * sd(iris$Sepal.Length),
               tolerance = 0.05)
  expect_equal(sd(s$Sepal.Width) / sd(iris$Sepal.Width), 2, tolerance = 0.05)
  expect_equal(mean(s$Petal.Length), mean(iris$Petal.Length), tolerance = 0.02)

  g <- data.frame(g = factor(rep(c("a", "b", "c"), c(70, 20, 10))))
  fg <- fit_synthesizer(g)
  hot <- prop.table(table(generate(fg, n = 20000, perturb = perturbation(temperature = 3), seed = 1)$g))
  cold <- prop.table(table(generate(fg, n = 20000, perturb = perturbation(temperature = 0.5), seed = 1)$g))
  expect_lt(hot[["a"]], 0.55)
  expect_gt(cold[["a"]], 0.85)

  xy <- as.data.frame(r_mvnorm(3000, sigma = make_corr(2, "exchangeable", 0.5), seed = 1))
  fxy <- fit_synthesizer(xy)
  r <- function(dep) cor(generate(fxy, n = 20000, perturb = perturbation(dependence = dep), seed = 1))[1, 2]
  expect_gt(r(1.6), 0.75)
  expect_lt(abs(r(0)), 0.03)
})

test_that("closeness and perturb combine, also in strata and helpers", {
  fit <- fit_synthesizer(iris, by = "Species")
  s <- generate(fit, closeness = 0.5, perturb = perturbation(shift = 1), seed = 1)
  expect_equal(nrow(s), 150L)
  expect_gt(mean(s$Sepal.Length), mean(iris$Sepal.Length))
  p <- SynthesizerPlus:::.combine_perturbation(SynthesizerPlus:::.closeness_perturbation(0.5),
                                               perturbation(noise = 0.75, swap = 0.5))
  expect_equal(p$noise, sqrt(0.75^2 + 0.75^2))
  expect_equal(p$swap, 1 - 0.75 * 0.5)
  expect_error(SynthesizerPlus:::.combine_perturbation(perturbation(noise = c(a = 1)),
                                                       perturbation(noise = c(b = 1))),
               "per-variable")
  a <- augment_data(mtcars, n = 50, closeness = 0.3, seed = 1)
  expect_equal(nrow(a), 82L)
  syn <- synthesize(mtcars, closeness = 0.5, seed = 1)
  expect_equal(dim(syn), dim(mtcars))
  expect_s3_class(synthesize(ldeaths, seed = 1), "ts")
})

test_that("calibrate_closeness() reaches a target", {
  fit <- fit_synthesizer(iris)
  cal <- calibrate_closeness(fit, iris, target = 4, metric = "dcr",
                             grid = seq(0, 1, 0.25), refine = 3)
  expect_s3_class(cal, "sp_calibration")
  expect_lt(abs(cal$achieved - 4), 0.6)
  expect_true(cal$closeness > 0 && cal$closeness < 1)
  expect_true(nrow(cal$path) %in% c(5L + 3L, 5L + 4L))
  expect_output(print(cal), "closeness")
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_s3_class(plot(cal), "ggplot")
  grDevices::dev.off()
  unlink(f)
  targets <- c(auc = 0.8, pmse = 5, marginal = 0.08)
  for (m in names(targets)) {
    cm <- suppressWarnings(calibrate_closeness(fit, iris, target = targets[[m]], metric = m,
                                               grid = c(0, 0.5, 1), refine = 1))
    expect_true(is.finite(cm$achieved))
  }
  expect_warning(calibrate_closeness(fit, iris, target = 100, metric = "dcr",
                                     grid = c(0, 1), refine = 0), "outside")
  expect_error(calibrate_closeness(1, iris, 1), "sp_synthesizer")
  expect_error(calibrate_closeness(fit, iris, "a"), "target")
  expect_error(calibrate_closeness(fit, iris, 1, grid = 0.5), "grid")
  fv <- fit_synthesizer(rnorm(100))
  expect_s3_class(suppressWarnings(
    calibrate_closeness(fv, data.frame(x = rnorm(100)), 0.2, metric = "marginal",
                        grid = c(0, 1), refine = 1)), "sp_calibration")
})

test_that("quadratic discriminator detects changes in correlation", {
  xy <- as.data.frame(r_mvnorm(1000, sigma = make_corr(2, "exchangeable", 0.8), seed = 1))
  ind <- as.data.frame(r_mvnorm(1000, sigma = diag(2), seed = 2))
  names(ind) <- names(xy)
  expect_lt(discriminator_auc(xy, ind, model = "linear", seed = 1), 0.6)
  expect_gt(discriminator_auc(xy, ind, model = "quadratic", seed = 1), 0.7)
})

test_that("calibrate_closeness() accepts a custom metric", {
  fit <- fit_synthesizer(iris)
  share <- function(real, syn) dcr(real, syn, seed = 1)$close_share
  cal <- suppressWarnings(calibrate_closeness(fit, iris, target = 0.05, metric = share,
                                              grid = c(0, 0.5, 1), refine = 2))
  expect_equal(cal$metric, "custom")
  expect_true(is.finite(cal$achieved))
  bad <- function(real, syn) c(1, 2)
  expect_error(calibrate_closeness(fit, iris, 1, metric = bad, grid = c(0, 1)), "single number")
})

test_that("closeness knob maps to the documented perturbation", {
  p <- SynthesizerPlus:::.closeness_perturbation(0.4)
  expect_equal(p$noise, 1.5 * 0.6)
  expect_equal(p$swap, 0.5 * 0.6)
  expect_equal(p$temperature, 1.6)
  expect_equal(p$dependence, 0.4)
  expect_true(SynthesizerPlus:::.is_identity_perturbation(SynthesizerPlus:::.closeness_perturbation(1)))
})

test_that("dependence argument and perturbation dependence multiply", {
  xy <- as.data.frame(r_mvnorm(3000, sigma = make_corr(2, "exchangeable", 0.6), seed = 1))
  fit <- fit_synthesizer(xy)
  r <- function(...) cor(generate(fit, n = 20000, seed = 1, ...))[1, 2]
  expect_equal(r(dependence = 0.5, perturb = perturbation(dependence = 2)), r(), tolerance = 0.03)
  expect_lt(r(dependence = 0.5), r() - 0.2)
})

test_that("calibration path and plot are consistent", {
  fit <- fit_synthesizer(mtcars)
  cal <- calibrate_closeness(fit, mtcars, target = 3, metric = "dcr", grid = c(0, 0.5, 1),
                             refine = 2)
  expect_true(all(diff(cal$path$closeness) > 0))
  expect_equal(cal$achieved, cal$path$value[cal$path$closeness == cal$closeness])
  expect_error(calibrate_closeness(fit, mtcars, 3, refine = -1), "refine")
})

test_that("numeric perturbations also act on discrete numeric variables", {
  d <- data.frame(k = rep(1:5, 40), x = rnorm(200))
  fit <- fit_synthesizer(d)
  expect_equal(fit$pooled$marginals$k$type, "discrete")
  s <- generate(fit, n = 20000, seed = 1, perturb = perturbation(shift = c(k = 1)))
  expect_equal(mean(s$k) - mean(d$k), sd(d$k), tolerance = 0.1)
  expect_true(all(s$k == round(s$k)))
  n <- generate(fit, n = 20000, seed = 1, perturb = perturbation(noise = c(k = 1)))
  expect_equal(mean(n$k), mean(d$k), tolerance = 0.02)
  expect_equal(sd(n$k), sd(d$k), tolerance = 0.1)
  # within strata, where variables often become discrete
  st <- fit_synthesizer(airquality, by = "Month", min_stratum_size = 10)
  sh <- generate(st, seed = 1, perturb = perturbation(shift = c(Temp = 1)))
  expect_gt(mean(sh$Temp), mean(airquality$Temp) + 3)
  # non-integer discrete values are left exactly unchanged without numeric settings
  r <- data.frame(v = rep(c(0.5, 1.25, 2.75), 10))
  expect_true(all(generate(fit_synthesizer(r), perturb = perturbation(swap = 0.5), seed = 1)$v %in%
                    c(0.5, 1.25, 2.75)))
})
