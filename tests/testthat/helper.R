mixed_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  data.frame(
    num = x,
    int = as.integer(round(50 + 10 * x + rnorm(n))),
    cnt = rpois(n, 2),
    fac = factor(ifelse(x + rnorm(n, sd = 0.5) > 0, "high", "low"),
                 levels = c("low", "high", "unused")),
    ord = factor(sample(c("S", "M", "L"), n, TRUE), levels = c("S", "M", "L"),
                 ordered = TRUE),
    chr = sample(c("a", "b", "c"), n, TRUE, prob = c(0.6, 0.3, 0.1)),
    lgl = runif(n) < 0.3,
    date = as.Date("2020-01-01") + sample(0:700, n, TRUE),
    time = as.POSIXct("2021-03-01 12:00:00", tz = "UTC") + round(runif(n) * 1e6),
    stringsAsFactors = FALSE
  )
}

# absolute-difference expectation (testthat 3e's expect_equal() tolerance is
# relative); works for vectors and matrices
expect_near <- function(object, expected, tol) {
  diff <- max(abs(as.numeric(object) - as.numeric(expected)))
  expect(is.finite(diff) && diff <= tol,
         sprintf("max |difference| = %.4g exceeds %.4g", diff, tol))
  invisible(object)
}
