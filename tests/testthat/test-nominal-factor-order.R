testthat::test_that("unordered factor level permutations leave fitted distribution invariant", {
  set.seed(3101)
  n <- 420L
  z <- stats::rnorm(n)
  group <- ifelse(z < -0.3, "B", ifelse(z < 0.8, "A", "C"))
  risk <- ifelse(z + stats::rnorm(n, sd = 0.9) > 0, "Yes", "No")
  d <- data.frame(
    x = z,
    y = 2 * z + 1.2 * (group == "C") + stats::rnorm(n),
    group = factor(group, levels = c("A", "B", "C", "UNUSED")),
    risk = factor(risk, levels = c("No", "Yes"))
  )
  d$y[seq(1L, n, 13L)] <- NA_real_
  d2 <- d
  d2$group <- factor(as.character(d$group), levels = c("C", "A", "UNUSED", "B"))
  d2$risk <- factor(as.character(d$risk), levels = c("Yes", "No"))

  f1 <- SynthesizerPlus::fit_synthesizer(d, seed = 1234, missing = "joint")
  f2 <- SynthesizerPlus::fit_synthesizer(d2, seed = 1234, missing = "joint")
  # Verify fitted marginals/correlation matrices are the same, not just a metric.
  testthat::expect_identical(f1$pooled$marginals$group$values,
                             f2$pooled$marginals$group$values)
  testthat::expect_identical(f1$pooled$marginals$risk$values,
                             f2$pooled$marginals$risk$values)
  testthat::expect_equal(f1$pooled$copula$R, f2$pooled$copula$R, tolerance = 1e-10)

  s1 <- SynthesizerPlus::generate(f1, n = 700, seed = 4321)
  s2 <- SynthesizerPlus::generate(f2, n = 700, seed = 4321)
  for (nm in names(d)) {
    testthat::expect_identical(as.character(s1[[nm]]), as.character(s2[[nm]]))
  }
  testthat::expect_identical(levels(s1$group), levels(d$group))
  testthat::expect_identical(levels(s2$group), levels(d2$group))
  testthat::expect_false("UNUSED" %in% as.character(s1$group))
})

testthat::test_that("ordered categories keep their statistically meaningful order", {
  d <- data.frame(
    y = seq_len(50L) + stats::rnorm(50L),
    severity = ordered(rep(c("low", "medium", "high", "medium", "low"), 10L),
                       levels = c("low", "medium", "high"))
  )
  fit <- SynthesizerPlus::fit_synthesizer(d, seed = 110)
  testthat::expect_identical(fit$pooled$marginals$severity$values,
                             c("low", "medium", "high"))
  syn <- SynthesizerPlus::generate(fit, n = 100L, seed = 111)
  testthat::expect_true(is.ordered(syn$severity))
  testthat::expect_identical(levels(syn$severity), levels(d$severity))
})
