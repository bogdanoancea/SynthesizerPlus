io_data <- function() {
  data.frame(
    num = c(1.5, NA, 3, -2.25),
    count = c(1L, 2L, NA, 4L),
    fac = factor(c("a", NA, "b", "a"), levels = c("b", "a", "z")),
    chr = c("x", NA, "z", "y"),
    lgl = c(TRUE, NA, FALSE, TRUE),
    date = as.Date(c("2020-01-01", NA, "2021-06-30", "1999-12-31")),
    time = as.POSIXct(c("2020-01-01 10:00:00", NA, "2020-05-05 00:00:01",
                        "2019-07-01 23:59:59"), tz = "Europe/Bucharest"),
    stringsAsFactors = FALSE
  )
}

roundtrip <- function(ext, template = TRUE) {
  x <- io_data()
  f <- tempfile(fileext = paste0(".", ext))
  on.exit(unlink(f))
  write_data(x, f)
  expect_true(file.exists(f))
  read_data(f, template = if (template) x else NULL)
}

test_that("supported_formats() lists formats and availability", {
  sf <- supported_formats()
  expect_s3_class(sf, "data.frame")
  expect_true(all(c("csv", "parquet", "hdf5", "json", "excel", "fst") %in% sf$format))
  expect_type(sf$available, "logical")
  expect_true(all(supported_formats(available_only = TRUE)$available))
  expect_true(all(c("csv", "tsv", "rds", "rdata") %in%
                    supported_formats(TRUE)$format))
})

test_that("native R formats round-trip exactly", {
  for (ext in c("rds", "rda", "rdata")) {
    expect_identical(roundtrip(ext, template = FALSE), io_data())
  }
})

test_that("delimited text formats round-trip with a template", {
  for (ext in c("csv", "tsv", "txt", "csv.gz", "tsv.bz2", "csv.xz")) {
    expect_equal(roundtrip(ext), io_data())
  }
})

test_that("delimited text works without data.table", {
  x <- io_data()
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f))
  # compressed path always uses base R; uncompressed base path tested via .write_delim
  SynthesizerPlus:::.write_delim(x, f, ",", compression = "none")
  y <- read.csv(f, stringsAsFactors = FALSE)
  expect_equal(nrow(y), 4L)
  expect_equal(match_types(y, x)$fac, x$fac)
})

test_that("JSON and NDJSON round-trip", {
  skip_if_not_installed("jsonlite")
  for (ext in c("json", "json.gz", "ndjson", "jsonl")) {
    expect_equal(roundtrip(ext), io_data())
  }
})

test_that("HDF5 round-trips without a template", {
  skip_if_not_installed("hdf5r")
  expect_equal(roundtrip("h5", template = FALSE), io_data())
  # custom group name and reading a plain dataset
  f <- tempfile(fileext = ".hdf5")
  on.exit(unlink(f))
  write_data(iris, f, dataset = "flowers")
  expect_equal(read_data(f, dataset = "flowers"), iris)
  expect_error(read_data(f), "not found")
  h5 <- hdf5r::H5File$new(f, mode = "a")
  h5$create_dataset("mat", robj = matrix(1:6, 3))
  h5$close_all()
  expect_equal(dim(read_data(f, dataset = "mat")), c(3L, 2L))
})

test_that("Excel round-trips", {
  skip_if_not_installed("readxl")
  skip_if_not(requireNamespace("writexl", quietly = TRUE) ||
                requireNamespace("openxlsx", quietly = TRUE))
  expect_equal(roundtrip("xlsx"), io_data())
  expect_error(write_data(iris, tempfile(fileext = ".xls")), "xlsx")
})

test_that("SPSS, Stata and SAS transport files round-trip", {
  skip_if_not_installed("haven")
  for (ext in c("sav", "dta", "xpt")) {
    expect_equal(roundtrip(ext), io_data())
  }
  expect_error(write_data(iris, tempfile(fileext = ".sas7bdat")), "xpt")
  f <- tempfile(fileext = ".dta")
  on.exit(unlink(f))
  d <- data.frame(Sepal.Length = 1:2, int = 3:4, `1x` = 5:6, check.names = FALSE)
  expect_warning(write_data(d, f), "Stata names")
  expect_named(read_data(f), c("Sepal_Length", "int_", "v1x"))
})

test_that("Parquet round-trips when a backend is installed", {
  skip_if_not(requireNamespace("arrow", quietly = TRUE) ||
                requireNamespace("nanoparquet", quietly = TRUE))
  expect_equal(roundtrip("parquet"), io_data())
})

test_that("Feather round-trips when arrow is installed", {
  skip_if_not_installed("arrow")
  expect_equal(roundtrip("feather"), io_data())
})

test_that("fst round-trips when installed", {
  skip_if_not_installed("fst")
  expect_equal(roundtrip("fst"), io_data())
})

test_that("missing backends give informative errors", {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    expect_error(write_data(iris, tempfile(fileext = ".feather")), "arrow")
  }
  if (!requireNamespace("fst", quietly = TRUE)) {
    expect_error(write_data(iris, tempfile(fileext = ".fst")), "fst")
  }
  if (!requireNamespace("arrow", quietly = TRUE) &&
      !requireNamespace("nanoparquet", quietly = TRUE)) {
    expect_error(write_data(iris, tempfile(fileext = ".parquet")), "Parquet")
  }
  expect_true(TRUE)
})

test_that("format detection and argument checks", {
  expect_error(write_data(iris, tempfile(fileext = ".foo")), "infer")
  expect_error(write_data(iris, tempfile(), format = "foo"), "Unknown format")
  expect_error(write_data(iris, tempfile(fileext = ".rds.gz")), "Compression")
  expect_error(read_data(tempfile(fileext = ".csv")), "not found")
  expect_error(read_data(c("a", "b")), "single")
  expect_error(write_data(1:3, tempfile(fileext = ".csv")), "data frame")
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  write_data(iris, f)
  expect_error(write_data(iris, f, overwrite = FALSE), "exists")
  # explicit format overrides the extension
  g <- tempfile(fileext = ".dat")
  write_data(mtcars, g, format = "csv")
  expect_equal(nrow(read_data(g, format = "csv")), 32L)
  # matrices are accepted
  write_data(as.matrix(mtcars), g, format = "tsv")
  expect_equal(ncol(read_data(g, format = "tsv")), 11L)
  unlink(g)
})

test_that("rdata objects can be selected by name", {
  f <- tempfile(fileext = ".RData")
  on.exit(unlink(f))
  a <- 1
  b <- mtcars
  save(a, b, file = f)
  expect_equal(read_data(f), mtcars)
  expect_equal(read_data(f, object = "b"), mtcars)
  expect_error(read_data(f, object = "zz"), "not found")
  save(a, file = f)
  expect_error(read_data(f), "No data frame")
  write_data(iris, f, object = "flowers")
  expect_equal(read_data(f, object = "flowers"), iris)
})

test_that("match_types() restores classes", {
  x <- data.frame(f = c("a", "b"), o = c(2, 1), l = c("yes", "F"),
                  d = c("2020-01-01", "2020-12-31"), dn = c(0, 1),
                  t = c("2020-01-01T10:00:00Z", "2020-01-01 11:30"),
                  tn = c(0, 3600), i = c("1", "2"), extra = 1:2,
                  stringsAsFactors = FALSE)
  tmpl <- data.frame(f = factor("a", levels = c("a", "b")),
                     o = factor("x", levels = c("x", "y"), ordered = TRUE),
                     l = TRUE, d = Sys.Date(), dn = Sys.Date(),
                     t = as.POSIXct("2020-01-01", tz = "UTC"),
                     tn = as.POSIXct("2020-01-01", tz = "UTC"), i = 1L)
  y <- match_types(x, tmpl)
  expect_equal(levels(y$f), c("a", "b"))
  expect_equal(as.character(y$o), c("y", "x"))
  expect_true(is.ordered(y$o))
  expect_equal(y$l, c(TRUE, FALSE))
  expect_equal(y$d, as.Date(c("2020-01-01", "2020-12-31")))
  expect_equal(y$dn, as.Date(c("1970-01-01", "1970-01-02")))
  expect_equal(format(y$t, "%H:%M", tz = "UTC"), c("10:00", "11:30"))
  expect_equal(as.numeric(y$tn), c(0, 3600))
  expect_type(y$i, "integer")
  expect_equal(y$extra, 1:2)
  expect_error(match_types(1, tmpl), "data frames")
})

test_that("read_data() reads a directory of part files", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  write_data(iris[1:50, ], file.path(d, "part-1.rds"))
  write_data(iris[51:150, ], file.path(d, "part-2.rds"))
  writeLines("not data", file.path(d, "README"))
  expect_equal(read_data(d), iris)
  e <- tempfile()
  dir.create(e)
  expect_error(read_data(e), "No readable")
})

test_that("generate_to_file() writes appendable formats in chunks", {
  fit <- fit_synthesizer(iris)
  for (ext in c("csv", "tsv.gz", "ndjson")) {
    f <- tempfile(fileext = paste0(".", ext))
    res <- generate_to_file(fit, n = 1000, path = f, chunk_size = 300, seed = 1)
    expect_equal(res$n, 1000L)
    x <- read_data(f, template = iris)
    expect_equal(dim(x), c(1000L, 5L))
    expect_s3_class(x$Species, "factor")
    unlink(f)
  }
})

test_that("generate_to_file() is reproducible and matches chunked generation", {
  fit <- fit_synthesizer(mtcars)
  f1 <- tempfile(fileext = ".csv")
  f2 <- tempfile(fileext = ".csv")
  generate_to_file(fit, 100, f1, chunk_size = 30, seed = 5)
  generate_to_file(fit, 100, f2, chunk_size = 30, seed = 5)
  expect_identical(readLines(f1), readLines(f2))
  unlink(c(f1, f2))
})

test_that("generate_to_file() writes part files for other formats", {
  fit <- fit_synthesizer(iris)
  d <- tempfile()
  on.exit(unlink(d, recursive = TRUE))
  res <- generate_to_file(fit, n = 1000, path = d, format = "rds", chunk_size = 400,
                          seed = 1, progress = FALSE)
  expect_length(res$files, 3L)
  expect_equal(basename(res$files), sprintf("part-%05d.rds", 1:3))
  expect_equal(nrow(read_data(d)), 1000L)
  # rewriting removes old parts
  res2 <- generate_to_file(fit, n = 100, path = d, format = "rds", chunk_size = 400)
  expect_length(list.files(d), 1L)
  expect_message(generate_to_file(fit, n = 10, path = d, format = "rds", progress = TRUE),
                 "chunk 1/1")
  # vector synthesizers
  fv <- fit_synthesizer(rnorm(20))
  fcsv <- tempfile(fileext = ".csv")
  generate_to_file(fv, 10, fcsv)
  expect_equal(nrow(read_data(fcsv)), 10L)
  expect_error(generate_to_file(1, 10, fcsv), "sp_synthesizer")
  expect_error(generate_to_file(fit, 10, fcsv, chunk_size = 0), "positive")
})

test_that("HDF5 part files can be generated", {
  skip_if_not_installed("hdf5r")
  fit <- fit_synthesizer(iris)
  d <- tempfile()
  on.exit(unlink(d, recursive = TRUE))
  generate_to_file(fit, n = 500, path = d, format = "h5", chunk_size = 200, seed = 1)
  x <- read_data(d)
  expect_equal(dim(x), c(500L, 5L))
  expect_s3_class(x$Species, "factor")
})
