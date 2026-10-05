# Guards against S3 methods that work inside the package but are not
# registered, and therefore silently fail for users calling from their session.

test_that("every S3 method defined in the package is registered", {
  ns <- asNamespace("SynthesizerPlus")
  reg <- getNamespaceInfo("SynthesizerPlus", "S3methods")
  registered <- paste(reg[, 1], reg[, 2], sep = ".")
  generics <- c("fit_synthesizer", "generate", "augment_data", "print", "summary",
                "plot", "simulate", "autoplot")
  fns <- ls(ns, all.names = TRUE)
  methods <- fns[vapply(fns, function(f) any(startsWith(f, paste0(generics, "."))),
                        logical(1))]
  expect_setequal(intersect(methods, registered), methods)
  expect_false(any(startsWith(getNamespaceExports("SynthesizerPlus"), ".")))
})

test_that("dispatch works from the user's global environment", {
  call_global <- function(expr) eval(expr, envir = globalenv())
  expect_s3_class(call_global(quote(SynthesizerPlus::fit_synthesizer(datasets::ldeaths))),
                  "sp_ts_synthesizer")
  expect_s3_class(call_global(quote(SynthesizerPlus::fit_synthesizer(datasets::iris))),
                  "sp_synthesizer")
  expect_s3_class(call_global(quote(SynthesizerPlus::fit_synthesizer(as.matrix(datasets::mtcars)))),
                  "sp_synthesizer")
  s <- call_global(quote(SynthesizerPlus::augment_data(datasets::ldeaths, nsim = 1, seed = 1)))
  expect_s3_class(s[[2]], "ts")
  fit <- fit_synthesizer(datasets::ldeaths)
  expect_s3_class(call_global(bquote(SynthesizerPlus::generate(.(fit), seed = 1))), "ts")
})
