# Edge cases and rarely used branches.

test_that("stratum sizes: proportional rounding and empty strata", {
  fit <- fit_synthesizer(iris, by = "Species")
  p7 <- generate(fit, n = 7, stratum_sizes = "proportional", seed = 1)
  expect_equal(nrow(p7), 7L)
  expect_true(all(table(p7$Species) %in% c(2L, 3L)))
  expect_equal(SynthesizerPlus:::.largest_remainder(10, c(0.55, 0.3, 0.15)), c(6L, 3L, 1L))
  expect_equal(sum(SynthesizerPlus:::.largest_remainder(7, rep(1 / 3, 3))), 7L)
  # tiny n with many strata: some strata get zero records
  s2 <- generate(fit, n = 2, seed = 5)
  expect_equal(nrow(s2), 2L)
  expect_lte(length(unique(s2$Species)), 2L)
})

test_that("identifier columns keep their type", {
  d <- data.frame(fid = factor(sprintf("F%02d", 1:20)), nid = as.numeric(1:20),
                  lid = 1:20 > 10, x = rnorm(20))
  syn <- generate(fit_synthesizer(d, id_cols = c("fid", "nid", "lid")), n = 12)
  expect_s3_class(syn$fid, "factor")
  expect_equal(syn$nid, as.numeric(1:12))
  expect_equal(syn$lid, 1:12)
  expect_equal(anyDuplicated(as.character(syn$fid)), 0L)
})

test_that("printing and summaries cover all model variants", {
  ft <- fit_synthesizer(mtcars, copula = "t", by = "am", id_cols = "carb")
  expect_output(print(ft), "strata")
  expect_output(print(ft), "id columns")
  expect_output(print(summary(ft)), "df =")
  expect_output(print(summary(fit_synthesizer(data.frame(a = NA_real_, b = 1:4)))), "-")
  expect_error(copula_correlation(fit_synthesizer(iris), stratum = "x"), "not stratified")
})

test_that("input validation of fit_synthesizer()", {
  expect_error(fit_synthesizer(data.frame(a = 1:3)[, 0, drop = FALSE]), "at least 2 rows|1 column")
  bad <- data.frame(a = 1:3, b = 1:3)
  names(bad) <- c("a", "a")
  expect_error(fit_synthesizer(bad), "unique")
  expect_error(fit_synthesizer(iris, by = 1), "character vector")
  expect_error(SynthesizerPlus:::.as_plain_df(1:3), "data frame")
  expect_error(SynthesizerPlus:::.fit_marginal(complex(real = 1:3), "z"), "unsupported class")
  i64 <- structure(c(1, 2, 3), class = "integer64")
  expect_error(fit_synthesizer(data.frame(x = 1:3, i = i64)), "Unsupported column class in: i")
})

test_that("copula estimation corner cases", {
  # t copula with very few complete rows uses a default df
  d <- data.frame(a = c(rnorm(8), rep(NA, 20)), b = rnorm(28))
  expect_equal(fit_synthesizer(d, copula = "t", missing = "drop")$pooled$copula$df, 10)
  # pd_method = "none" on a valid matrix
  expect_s3_class(fit_synthesizer(iris, pd_method = "none"), "sp_synthesizer")
  expect_error(fit_synthesizer(data.frame(a = 1:30, b = 1:30, c = -(1:30)), pd_method = "none"),
               "not positive definite")
  # t df estimation subsamples large inputs
  big <- as.data.frame(r_mvt(6000, df = 4, sigma = make_corr(2, "exchangeable", 0.5), seed = 1))
  expect_lte(fit_synthesizer(big, copula = "t")$pooled$copula$df, 8)
  # a single-observation continuous marginal
  m <- SynthesizerPlus:::.fit_marginal(c(1.5, NA), "x", discrete_threshold = 0L)
  expect_equal(SynthesizerPlus:::.marginal_quantile(m, c(0.1, 0.9)), c(1.5, 1.5))
})

test_that("time-series corner cases", {
  expect_error(fit_synthesizer(ts(matrix(letters[1:20], 10))), "numeric")
  x <- ts(c(NA, 1, rep(NA, 10)))
  expect_error(suppressWarnings(fit_synthesizer(x)), "fewer than 2")
  # white noise: VAR order 0, innovations are the centred data
  wn <- ts(rnorm(200))
  f <- fit_synthesizer(wn)
  expect_equal(f$var$order, 0L)
  expect_length(generate(f, n = 50, seed = 1), 50L)
})

test_that("metric helpers handle degenerate inputs", {
  expect_true(is.na(SynthesizerPlus:::.ks_stat(c(NA, NA), 1:3)))
  expect_true(is.na(SynthesizerPlus:::.wasserstein(numeric(0), 1:3)))
  expect_equal(SynthesizerPlus:::.wasserstein(rep(1, 5), rep(2, 5)), 1)   # zero sd -> unscaled
  expect_true(is.na(SynthesizerPlus:::.tvd(c(NA, NA), c(NA, NA))))
  expect_equal(SynthesizerPlus:::.cramers_v(rep("a", 5), letters[1:5]), 0)
  expect_equal(SynthesizerPlus:::.corr_ratio(rep(1, 6), rep(c("a", "b"), 3)), 0)
  expect_equal(SynthesizerPlus:::.cramers_v(c("a", "b", "a", "b"), c("x", "x", "y", "y")), 0)
  expect_true(is.na(SynthesizerPlus:::.auc(c(1, 1), c(0.2, 0.4))))
  expect_equal(dim(SynthesizerPlus:::.prep_features(data.frame(a = rep(NA_real_, 3)))), c(3L, 0L))
  X <- SynthesizerPlus:::.prep_features(as.data.frame(matrix(rnorm(250), 10)))
  expect_equal(ncol(SynthesizerPlus:::.add_quadratic(X)) - ncol(X), 20L + 190L)
  expect_identical(SynthesizerPlus:::.add_quadratic(
    SynthesizerPlus:::.prep_features(data.frame(g = c("a", "b", "a")))), 
    SynthesizerPlus:::.prep_features(data.frame(g = c("a", "b", "a"))))
  expect_true(is.matrix(association_matrix(as.matrix(mtcars[1:3]))))
  far <- data.frame(x = c(1, 1), k = c(5, 5))
  expect_equal(dcr(far, far)$exact_match_rate, 1)
})

test_that("disclosure_risk() without real uniques and with all-missing targets", {
  real <- data.frame(k = rep(1:3, each = 2), y = c(1, 2, 3, 4, 5, 6))
  r <- disclosure_risk(real, real, keys = "k", target = "y")
  expect_equal(r$n_real_uniques, 0L)
  expect_equal(r$replicated_uniques, 0)
  expect_true(is.na(r$attribute$cap_uniques))
  r2 <- disclosure_risk(data.frame(k = 1:3, y = NA_real_), data.frame(k = 1:3, y = NA_real_),
                        keys = "k", target = "y")
  expect_true(is.na(r2$attribute$cap))
})

test_that("plot helpers accept matrices, ts and comparison objects", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plot_marginals(as.matrix(mtcars[1:3]), max_vars = 2), "ggplot")
  expect_equal(length(unique(plot_marginals(mtcars, max_vars = 3)$data$variable)), 3L)
  expect_s3_class(plot_qq(EuStockMarkets[1:100, ], EuStockMarkets[101:200, ]), "ggplot")
  expect_error(plot_marginals(list(1)), "data frames")
  expect_error(plot_marginals(data.frame(a = 1)[0], NULL), "No variables")
  cmp <- compare_synthetic(iris, synthesize(iris, seed = 1), seed = 1)
  expect_s3_class(plot_dcr(cmp), "ggplot")
})

test_that("rarely used I/O branches", {
  skip_if_not_installed("hdf5r")
  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f))
  write_data(mtcars, f)
  write_data(iris, f)                         # overwrite an existing HDF5 file
  expect_equal(read_data(f), iris)
  d <- data.frame(t = as.difftime(c(1, 2), units = "mins"), z = complex(real = 1:2))
  write_data(d, f)            # classed numeric stored as number, unsupported class as text
  back <- read_data(f)
  expect_type(back$t, "double")
  expect_type(back$z, "character")
  expect_equal(match_types(back, d)$t, d$t)
  h5 <- hdf5r::H5File$new(f, mode = "w")
  g <- h5$create_group("data")
  g$create_dataset("a", robj = 1:3)
  g$create_dataset("b", robj = 1:2)
  h5$close_all()
  expect_error(read_data(f), "different lengths")
})

test_that("match_types() handles classed and unsupported templates; generate_to_file() checks paths", {
  x <- data.frame(a = c("1", "2"), z = c("1+0i", "2+0i"), stringsAsFactors = FALSE)
  tmpl <- data.frame(a = as.difftime(c(1, 2), units = "secs"), z = complex(real = 1:2))
  y <- match_types(x, tmpl)
  expect_equal(y$a, as.difftime(c(1, 2), units = "secs"))
  expect_identical(y$z, x$z)
  f <- tempfile()
  writeLines("x", f)
  on.exit(unlink(f))
  expect_error(generate_to_file(fit_synthesizer(iris), 10, f, format = "rds"), "not a directory")
  expect_equal(SynthesizerPlus:::.canonical_ext("rdata"), "rda")
  expect_equal(SynthesizerPlus:::.canonical_ext("excel"), "xlsx")
  expect_equal(SynthesizerPlus:::.canonical_ext("parquet"), "parquet")
})

test_that("unused-dots helper handles unnamed arguments", {
  expect_warning(SynthesizerPlus:::.warn_unused_dots(list(1), "f"), "<unnamed>")
  expect_null(SynthesizerPlus:::.warn_unused_dots(list(), "f"))
})

test_that("multivariate generators validate dimensions", {
  expect_error(r_mvt(5, df = 3, mean = 1:3, sigma = diag(2)), "incompatible")
})
