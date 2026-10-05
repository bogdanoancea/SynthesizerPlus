render <- function(p) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  print(p)
  invisible(TRUE)
}

syn_iris <- synthesize(iris, seed = 1)

test_that("plot_marginals() returns ggplots for each variable kind", {
  both <- plot_marginals(iris, syn_iris)
  expect_s3_class(both, "sp_plot_list")
  expect_named(both, c("numeric", "categorical"))
  expect_true(render(both))
  for (k in c("density", "histogram", "ecdf")) {
    p <- plot_marginals(iris, syn_iris, vars = c("Sepal.Length", "Petal.Width"), kind = k)
    expect_s3_class(p, "ggplot")
    expect_true(render(p))
  }
  only_real <- plot_marginals(iris, vars = "Sepal.Length")
  expect_s3_class(only_real, "ggplot")
  expect_true(render(plot_marginals(rnorm(50), rnorm(50))))
  expect_error(plot_marginals(iris, vars = "zz"), "not found")
})

test_that("plot_categorical() lumps rare levels and shows missing values", {
  d <- data.frame(g = c(letters, NA)[sample(27, 300, TRUE)])
  p <- plot_categorical(d, d, max_levels = 5)
  expect_s3_class(p, "ggplot")
  vals <- unique(p$data$value)
  expect_true("(other)" %in% vals)
  expect_lte(length(vals), 5L)
  expect_true(render(p))
  expect_error(plot_categorical(mtcars), "No categorical")
})

test_that("plot_association() draws real, synthetic and difference panels", {
  p <- plot_association(iris, syn_iris)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$panel),
               c("Real", "Synthetic", "Difference (synthetic - real)"))
  expect_true(render(p))
  p1 <- plot_association(mtcars, show_values = FALSE)
  expect_equal(levels(p1$data$panel), "Real")
  expect_true(render(p1))
  expect_error(plot_association(iris, vars = "Species"), "two variables")
})

test_that("plot_pairs(), plot_qq() and plot_dcr() render", {
  p <- plot_pairs(iris, syn_iris, vars = c("Sepal.Length", "Petal.Length"))
  expect_s3_class(p, "ggplot")
  expect_true(render(p))
  expect_true(render(plot_pairs(as.data.frame(r_copula(200, "clayton", theta = 2)))))
  expect_error(plot_pairs(iris, vars = "Sepal.Length"), "two numeric")
  q <- plot_qq(iris, syn_iris)
  expect_s3_class(q, "ggplot")
  expect_true(render(q))
  expect_error(plot_qq(data.frame(g = letters), data.frame(g = letters)), "No numeric")
  d <- plot_dcr(dcr(iris, syn_iris, seed = 1))
  expect_s3_class(d, "ggplot")
  expect_true(render(d))
  expect_error(plot_dcr(list()), "dcr")
})

test_that("time-series plots render", {
  fit <- fit_synthesizer(mdeaths)
  syn <- generate(fit, seed = 1)
  expect_true(render(plot_ts(mdeaths, syn)))
  expect_true(render(plot_ts(EuStockMarkets)))
  expect_true(render(plot_ts(rnorm(20), rnorm(20))))
  a <- plot_acf(mdeaths, syn, lag_max = 12)
  expect_s3_class(a, "ggplot")
  expect_equal(max(a$data$lag), 12)
  expect_true(render(a))
})

test_that("autoplot() and plot() work for comparisons", {
  cmp <- compare_synthetic(iris, syn_iris, seed = 1)
  for (tp in c("marginals", "association", "qq", "pairs", "dcr")) {
    p <- ggplot2::autoplot(cmp, type = tp)
    expect_true(inherits(p, "ggplot") || inherits(p, "sp_plot_list"))
  }
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_invisible(plot(cmp, type = "association"))
  grDevices::dev.off()
  unlink(f)
  no_priv <- compare_synthetic(iris, syn_iris, metrics = "marginal")
  expect_error(ggplot2::autoplot(no_priv, type = "dcr"), "privacy")
})
