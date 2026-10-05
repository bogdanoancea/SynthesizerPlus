# Cross-cutting behaviour: argument checking, silent-failure guards and
# consistency of the API across data types.

test_that("misspelled or misplaced arguments are reported", {
  expect_warning(fit_synthesizer(iris, copla = "t"), "Unused argument.*copla")
  expect_warning(fit_synthesizer(iris, closeness = 0.5), "arguments of generate")
  fit <- fit_synthesizer(iris)
  expect_warning(generate(fit, n = 5, closness = 0.5), "Unused argument.*closness")
  expect_warning(fit_synthesizer(ldeaths, maxlag = 3), "Unused argument.*maxlag")
  expect_warning(fit_synthesizer(ldeaths, method = "iid", foo = 1), "Unused argument.*foo")
  expect_warning(augment_data(mtcars, n = 5, bye = "cyl"), "Unused argument.*bye")
  expect_warning(synthesize(iris, n = 5, dependnce = 0.5), "Unused argument.*dependnce")
})

test_that("closeness controls are rejected loudly for time series", {
  fit <- fit_synthesizer(ldeaths)
  expect_warning(s <- generate(fit, closeness = 0.5, seed = 1), "not supported for time-series")
  expect_identical(s, generate(fit, seed = 1))
  expect_warning(generate(fit, perturb = perturbation(noise = 1)), "not supported")
  expect_warning(synthesize(ldeaths, closeness = 0.5, seed = 1), "not supported")
  expect_warning(augment_data(ldeaths, nsim = 1, closeness = 0.5), "tabular data only")
  expect_warning(simulate(fit, nsim = 1, closeness = 0.5), "not supported")
  expect_silent(synthesize(ldeaths, seed = 1))
})

test_that("small continuous samples are modelled as continuous", {
  set.seed(1)
  x <- rnorm(15)
  fit <- fit_synthesizer(x)
  expect_equal(fit$pooled$marginals$x$type, "continuous")
  syn <- generate(fit, n = 200, seed = 1)
  expect_gt(mean(!syn %in% x), 0.9)                # new values, not resampled ones
  # integer-valued and repeating variables stay discrete
  expect_equal(fit_synthesizer(data.frame(k = c(1:10, 1:5)))$pooled$marginals$k$type, "discrete")
  expect_equal(fit_synthesizer(data.frame(r = rep(c(0.5, 1.25, 2.75), 4)))$pooled$marginals$r$type,
               "discrete")
  expect_equal(fit_synthesizer(data.frame(k = 1:15))$pooled$marginals$k$type, "discrete")
  # a single distinct value is always reproduced
  expect_true(all(generate(fit_synthesizer(data.frame(c = rep(2.5, 5))), n = 10)$c == 2.5))
})

test_that("perturbing strata or identifiers gives an explanatory error", {
  expect_error(generate(fit_synthesizer(iris, by = "Species"),
                        perturb = perturbation(swap = c(Species = 0.5))),
               "a 'by' variable")
  expect_error(generate(fit_synthesizer(data.frame(id = 1:5, x = rnorm(5)), id_cols = "id"),
                        perturb = perturbation(noise = c(id = 1))),
               "identifier column")
})

test_that("all input types round-trip through fit/generate/closeness", {
  inputs <- list(
    vector = rnorm(40),
    factor = factor(sample(c("a", "b", "c"), 40, TRUE)),
    character = sample(c("x", "y"), 40, TRUE),
    logical = sample(c(TRUE, FALSE, NA), 40, TRUE),
    date = as.Date("2024-01-01") + sample(0:300, 40),
    matrix = as.matrix(mtcars[1:4]),
    data.frame = mtcars
  )
  for (nm in names(inputs)) {
    x <- inputs[[nm]]
    fit <- fit_synthesizer(x)
    for (cl in c(1, 0.5, 0)) {
      s <- generate(fit, n = 7, closeness = cl, seed = 1)
      expect_equal(class(s), class(x), info = paste(nm, cl))
      expect_equal(NROW(s), 7L, info = nm)
    }
  }
})

test_that("augment_data() extends vectors", {
  v <- augment_data(rnorm(20), n = 30, seed = 1)
  expect_type(v, "double")
  expect_length(v, 50L)
  f <- augment_data(factor(c("a", "b", "b")), n = 5, seed = 1)
  expect_s3_class(f, "factor")
  expect_length(f, 8L)
  expect_error(augment_data(list(1)), "Cannot augment")
})

test_that("closeness and perturbation work with every model option", {
  opts <- list(list(copula = "t"), list(copula = "independence"),
               list(missing = "independent"), list(missing = "drop"),
               list(interpolation = "spline"), list(knots = 20),
               list(by = "Month", min_stratum_size = 10))
  for (o in opts) {
    fit <- do.call(fit_synthesizer, c(list(airquality), o))
    s <- generate(fit, closeness = 0.4, seed = 1,
                  perturb = perturbation(shift = c(Temp = 1), temperature = 2))
    expect_equal(dim(s), dim(airquality), info = paste(names(o), collapse = ","))
    expect_gt(mean(s$Temp), mean(airquality$Temp))
  }
})

test_that("per-variable perturbations leave other variables untouched", {
  fit <- fit_synthesizer(iris)
  a <- generate(fit, seed = 3)
  b <- generate(fit, seed = 3, perturb = perturbation(scale = c(Sepal.Width = 3)))
  expect_identical(a$Petal.Length, b$Petal.Length)
  expect_identical(a$Species, b$Species)
  expect_false(identical(a$Sepal.Width, b$Sepal.Width))
})

test_that("bounds options are respected for every perturbation", {
  d <- data.frame(x = rexp(500) + 1)
  fit <- fit_synthesizer(d)
  p <- function(b) generate(fit, n = 5000, seed = 1,
                            perturb = perturbation(noise = 1, shift = -2, scale = 3, bounds = b))$x
  expect_true(all(p("observed") >= min(d$x) & p("observed") <= max(d$x)))
  expect_true(all(p("sign") >= 0))
  expect_lt(min(p("none")), 0)
  neg <- fit_synthesizer(data.frame(x = -rexp(200)))
  expect_true(all(generate(neg, n = 2000, perturb = perturbation(shift = 3), seed = 1)$x <= 0))
})

test_that("generation is reproducible with every control", {
  fit <- fit_synthesizer(census_sim[1:500, ], id_cols = "person_id", by = "sex")
  args <- list(closeness = 0.7, perturb = perturbation(noise = c(income = 0.4), swap = 0.1),
               dependence = 0.8, seed = 11)
  expect_identical(do.call(generate, c(list(fit), args)), do.call(generate, c(list(fit), args)))
})
