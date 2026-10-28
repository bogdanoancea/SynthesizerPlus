test_that("nominal model is conditional, unordered and reproducible", {
    skip_if_not_installed("nnet")
    set.seed(615)
    n <- 900L
    x <- rnorm(n)
    z <- rnorm(n)
    score <- cbind(0, 1.1 * x - .4 * z, -0.8 * x + 0.7 * z)
    pr <- exp(score - apply(score, 1L, max))
    pr <- pr / rowSums(pr)
    lab <- apply(pr, 1L, function(p) sample.int(3L, 1L, prob = p))
    d <- data.frame(x = x, z = z, category = factor(c("red", "blue", "green")[lab]))
    f <- fit_synthesizer(d, nominal_target = "category", seed = 12)
    expect_false("category" %in% colnames(f$pooled$copula$R))
    a <- generate(f, n = 3000L, seed = 29)
    b <- generate(f, n = 3000L, seed = 29)
    expect_identical(a, b)
    expect_s3_class(a$category, "factor")
    expect_identical(levels(a$category), levels(d$category))
    expect_equal(sort(unique(as.character(a$category))), sort(unique(as.character(d$category))))
    expect_true(all(is.finite(a$x)))
    # Distinct conditional class probabilities were specified in the DGP.
    expect_gt(abs(mean(a$category[a$x > 1] == "blue") -
                      mean(a$category[a$x < -1] == "blue")), 0.1)
})

test_that("nominal synthesis is equivariant under category relabelling", {
    skip_if_not_installed("nnet")
    set.seed(616)
    d <- data.frame(x = rnorm(900), y = factor(sample(c("a", "b", "c"), 900, TRUE)))
    d2 <- d
    d2$y <- factor(c(a = "zeta", b = "mu", c = "alpha")[as.character(d$y)],
                   levels = c("alpha", "mu", "zeta"))
    f1 <- fit_synthesizer(d, nominal_target = "y", seed = 21)
    f2 <- fit_synthesizer(d2, nominal_target = "y", seed = 21)
    expect_equal(unname(copula_correlation(f1)), unname(copula_correlation(f2)), tolerance = 1e-10)
    s1 <- generate(f1, 2000L, seed = 31)
    s2 <- generate(f2, 2000L, seed = 31)
    expect_equal(s1$x, s2$x)
    # Probabilities must not depend on lexical ordering. Tiny numerical
    # differences from multinomial parameterisation are allowed.
    original <- as.character(s1$y)
    inverse <- c(zeta = "a", mu = "b", alpha = "c")
    mapped <- unname(inverse[as.character(s2$y)])
    expect_lt(mean(original != mapped), 0.04)
})

test_that("nominal model preserves target missingness and validates arguments", {
    skip_if_not_installed("nnet")
    set.seed(617)
    x <- rnorm(1200)
    y <- factor(ifelse(x + rnorm(1200) > 0, "A", "B"))
    y[seq.int(1, 1200, by = 8)] <- NA
    d <- data.frame(x, y)
    f <- fit_synthesizer(d, nominal_target = "y", missing = "joint", seed = 2)
    s <- generate(f, 4000L, seed = 3)
    # Missingness is generated stochastically; test its probability, not
    # exact equality of empirical proportions.  A four-SE bound is conservative
    # under the fitted Bernoulli missingness model.
    p_expected <- mean(is.na(d$y))
    p_observed <- mean(is.na(s$y))
    se <- sqrt(p_expected * (1 - p_expected) / nrow(s))
    expect_lt(abs(p_observed - p_expected), 4 * se)
    expect_error(fit_synthesizer(d, nominal_target = "x"), "unordered factor")
    expect_error(fit_synthesizer(d, nominal_target = "y", nonlinear_target = "x"),
                 "mutually exclusive")
    expect_error(generate(f, 100L, dependence = .5), "Conditional synthesis")
})

test_that("Student-t degrees-of-freedom grid matches reference likelihood", {
    set.seed(618)
    R <- matrix(c(1, .4, .1, .4, 1, .2, .1, .2, 1), 3, 3)
    U <- r_copula(250, "t", corr = R, df = 4, seed = 619)
    candidates <- c(2, 3, 4, 5, 8, 12)
    X <- pmin(pmax(U, 1e-10), 1 - 1e-10)
    inv <- solve(R)
    det <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
    ref <- vapply(candidates, function(nu) {
        z <- stats::qt(X, df = nu)
        q <- rowSums((z %*% inv) * z)
        lmv <- lgamma((nu + 3) / 2) - lgamma(nu / 2) - 3 / 2 * log(nu * pi) -
            .5 * det - (nu + 3) / 2 * log1p(q / nu)
        sum(lmv - rowSums(stats::dt(z, df = nu, log = TRUE)))
    }, numeric(1))
    faster <- SynthesizerPlus:::.estimate_t_df(U, R, candidates, max_rows = 500L)
    expect_identical(as.numeric(faster$df), as.numeric(candidates[which.max(ref)]))
    fit <- fit_synthesizer(as.data.frame(U), copula = "t", t_profile = TRUE, seed = 1)
    expect_named(fit$pooled$copula$profile,
                 c("kendall_seconds", "discrete_latent_seconds",
                   "correlation_repair_seconds", "df_grid_seconds",
                   "other_seconds", "total_seconds"))
    expect_true(all(is.finite(fit$pooled$copula$profile)))
    expect_gte(fit$pooled$copula$profile[["total_seconds"]], 0)
    control <- fit_synthesizer(as.data.frame(U), copula = "t", seed = 1)
    expect_equal(control$pooled$copula$R, fit$pooled$copula$R)
    expect_equal(control$pooled$copula$df, fit$pooled$copula$df)
})

test_that("sparse nominal stratum retains baseline copula", {
    skip_if_not_installed("nnet")
    set.seed(622)
    n <- 300L
    d <- data.frame(k = factor(c(rep("A", 240), rep("B", 60))),
                    x = rnorm(n),
                    group = factor(sample(c("red", "blue", "green"), n, TRUE)))
    # Force stratum B to a single observed outcome: no multinomial model.
    d$group[d$k == "B"] <- "red"
    expect_warning(
        fit <- fit_synthesizer(d, by = "k", nominal_target = "group", seed = 23),
        "nominal conditional model could not be fitted in 1 of 2")
    expect_identical(as.character(fit$strata$conditional_fallback$k), "B")
    mods <- fit$strata$models
    sparse <- mods[[which(as.character(fit$strata$values$k) == "B")]]
    expect_null(sparse$nominal)
    expect_true("group" %in% colnames(sparse$copula$R))
    syn <- generate(fit, 500, seed = 24)
    expect_equal(nrow(syn), 500L)
})
