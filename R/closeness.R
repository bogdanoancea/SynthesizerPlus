# Controlling how close synthetic data are to the real data -------------------

#' Perturbations of the synthetic data distribution
#'
#' @description
#' `perturbation()` describes how synthetic data should *deviate* from the
#' fitted model. It is passed to [generate()] (argument `perturb`) and has
#' two groups of settings.
#'
#' **Record-level** settings move synthetic records away from the real
#' records while (approximately) keeping the distributions:
#'
#' * `noise` -- smoothing of numeric variables. Each value is replaced by
#'   \eqn{\mu + (x - \mu + h s \varepsilon) / \sqrt{1 + h^2}} with
#'   \eqn{\varepsilon \sim N(0, 1)}, where \eqn{h} is `noise` and \eqn{\mu},
#'   \eqn{s} are the mean and standard deviation of the variable. Mean and
#'   variance are preserved; the distribution becomes smoother and values no
#'   longer coincide with observed ones.
#' * `swap` -- probability that a value of a categorical or discrete variable
#'   is replaced by an independent draw from its marginal distribution. The
#'   marginal distribution is preserved; the association with other variables
#'   is weakened for the affected records.
#'
#' **Distribution-level** settings deliberately change the distributions,
#' e.g. for stress tests or to simulate a domain shift:
#'
#' * `shift` -- shift of numeric variables, in units of their standard
#'   deviation;
#' * `scale` -- multiplier of the spread of numeric variables around their
#'   median;
#' * `temperature` -- reshapes the probabilities of categorical and discrete
#'   variables, \eqn{p_k \propto p_k^{1/T}}: `T > 1` makes them more uniform,
#'   `T < 1` concentrates them on the most frequent values;
#' * `dependence` -- strength of the dependence between variables: 0 gives
#'   independent variables, 1 the fitted dependence, values above 1 (up to 2)
#'   strengthen it.
#'
#' Numeric variables that only take whole numbers (counts, ages, dates)
#' stay whole numbers; they are rounded stochastically after perturbation, so
#' that shifts and noise are not biased by rounding.
#'
#' `noise`, `swap`, `shift`, `scale` and `temperature` can be single numbers
#' (applied to every eligible variable) or named vectors that apply to the
#' named variables only.
#'
#' @param noise Non-negative smoothing bandwidth for numeric variables.
#' @param swap Probability in \[0, 1\] of independent re-drawing for
#'   categorical and discrete variables.
#' @param shift Shift of numeric variables in standard deviations.
#' @param scale Positive spread multiplier for numeric variables.
#' @param temperature Positive temperature for categorical and discrete
#'   probabilities.
#' @param dependence Dependence strength in \[0, 2\].
#' @param bounds How to keep perturbed numeric values in range:
#'   `"sign"` (default; variables that are never negative, or never positive,
#'   keep that sign), `"observed"` (clip to the observed range) or `"none"`.
#' @return An object of class `sp_perturbation`.
#' @seealso [generate()], [calibrate_closeness()]
#' @examples
#' fit <- fit_synthesizer(iris)
#'
#' # records further from real ones, same distributions
#' far <- generate(fit, perturb = perturbation(noise = 0.5, swap = 0.2), seed = 1)
#'
#' # a shifted scenario: longer petals, more uniform species, weaker dependence
#' p <- perturbation(shift = c(Petal.Length = 1), temperature = 3, dependence = 0.5)
#' p
#' scen <- generate(fit, perturb = p, seed = 1)
#' tapply(scen$Petal.Length, scen$Species, mean)
#' @export
perturbation <- function(noise = 0, swap = 0, shift = 0, scale = 1, temperature = 1,
                         dependence = 1, bounds = c("sign", "observed", "none")) {
  bounds <- match.arg(bounds)
  chk <- function(x, name, ok, msg) {
    if (!is.numeric(x) || !length(x) || anyNA(x) || !all(ok(x))) {
      stop(sprintf("'%s' must be %s.", name, msg), call. = FALSE)
    }
    if (length(x) > 1L && (is.null(names(x)) || any(names(x) == ""))) {
      stop(sprintf("'%s' must be a single number or a named vector.", name), call. = FALSE)
    }
    x
  }
  structure(list(
    noise = chk(noise, "noise", function(v) v >= 0 & is.finite(v), "non-negative"),
    swap = chk(swap, "swap", function(v) v >= 0 & v <= 1, "in [0, 1]"),
    shift = chk(shift, "shift", is.finite, "finite"),
    scale = chk(scale, "scale", function(v) v > 0 & is.finite(v), "positive"),
    temperature = chk(temperature, "temperature", function(v) v > 0 & is.finite(v), "positive"),
    dependence = chk(dependence, "dependence", function(v) length(v) == 1L & v >= 0 & v <= 2,
                     "a single number in [0, 2]"),
    bounds = bounds
  ), class = "sp_perturbation")
}

#' @export
print.sp_perturbation <- function(x, ...) {
  fmt <- function(lab, v) {
    if (length(v) == 1L && is.null(names(v))) return(paste(lab, format(signif(v, 3))))
    sprintf("%s (%s)", lab, paste(sprintf("%s = %s", names(v), format(signif(v, 3))),
                                  collapse = ", "))
  }
  cat("<sp_perturbation>\n")
  cat("  record level      :", fmt("noise", x$noise), "|", fmt("swap", x$swap), "\n")
  cat("  distribution level:", fmt("shift", x$shift), "|", fmt("scale", x$scale), "|",
      fmt("temperature", x$temperature), "|", fmt("dependence", x$dependence), "\n")
  cat("  bounds            :", x$bounds, "\n")
  invisible(x)
}

# Map the closeness knob (1 = as close as possible) to a perturbation
.closeness_perturbation <- function(closeness) {
  d <- 1 - closeness
  perturbation(noise = 1.5 * d, swap = 0.5 * d, temperature = 1 + d,
               dependence = closeness)
}

# Combine two perturbations: noise adds in quadrature, swaps compound,
# temperatures and dependence multiply; shift/scale come from `b`
.combine_perturbation <- function(a, b) {
  if (is.null(b)) return(a)
  if (!inherits(b, "sp_perturbation")) {
    stop("'perturb' must be created with perturbation().", call. = FALSE)
  }
  comb <- function(x, y, f) {
    if (length(x) == 1L && is.null(names(x))) return(f(x, y))
    if (length(y) == 1L && is.null(names(y))) return(f(y, x))
    stop("Cannot combine two per-variable settings.", call. = FALSE)
  }
  b$noise <- comb(a$noise, b$noise, function(s, v) sqrt(s^2 + v^2))
  b$swap <- comb(a$swap, b$swap, function(s, v) 1 - (1 - s) * (1 - v))
  b$temperature <- comb(a$temperature, b$temperature, function(s, v) s * v)
  b$dependence <- min(2, a$dependence * b$dependence)
  b
}

# value of a (possibly named) setting for variable `nm`
.pv <- function(x, nm, default) {
  if (length(x) == 1L && is.null(names(x))) return(x)
  if (nm %in% names(x)) return(x[[nm]])
  default
}

.is_identity_perturbation <- function(p) {
  all(p$noise == 0) && all(p$swap == 0) && all(p$shift == 0) && all(p$scale == 1) &&
    all(p$temperature == 1) && p$dependence == 1
}

.check_pert_names <- function(p, cols, by = NULL, id_cols = NULL) {
  for (f in c("noise", "swap", "shift", "scale", "temperature")) {
    nm <- names(p[[f]])
    bad <- setdiff(nm, cols)
    if (length(bad)) {
      why <- ifelse(bad %in% by, " (a 'by' variable: strata cannot be perturbed)",
                    ifelse(bad %in% id_cols, " (an identifier column)", " (not in the data)"))
      stop(sprintf("Cannot perturb variable(s) in '%s': %s", f,
                   paste0(bad, why, collapse = ", ")), call. = FALSE)
    }
  }
  invisible(NULL)
}

# Draw values for one marginal under a perturbation
.perturbed_values <- function(m, u, p) {
  nm <- m$name
  if (m$type == "empty") return(.marginal_quantile(m, u))
  int_valued <- isTRUE(m$integer_valued)
  if (m$type %in% c("categorical", "discrete")) {
    probs0 <- m$probs
    temp <- .pv(p$temperature, nm, 1)
    if (temp != 1) {
      w <- m$probs^(1 / temp)
      m$probs <- w / sum(w)
    }
    v <- .marginal_quantile(m, u)
    sw <- .pv(p$swap, nm, 0)
    if (sw > 0) {
      idx <- which(stats::runif(length(u)) < sw)
      if (length(idx)) v[idx] <- .marginal_quantile(m, stats::runif(length(idx)))
    }
    if (m$type == "categorical") return(v)
    # summary statistics of the (unperturbed) discrete distribution
    vals <- m$values
    pr <- probs0
    mu <- sum(vals * pr)
    s <- sqrt(sum(pr * (vals - mu)^2))
    med <- vals[which(cumsum(pr) >= 0.5)[1L]]
    lo <- min(vals)
    hi <- max(vals)
  } else {
    m$integer_valued <- FALSE
    v <- .marginal_quantile(m, u)
    q <- m$quantiles
    mu <- mean(q)
    s <- if (length(q) > 1L) stats::sd(q) else 0
    med <- stats::median(q)
    lo <- min(q)
    hi <- max(q)
  }
  # numeric perturbations (continuous and discrete numeric variables)
  h <- .pv(p$noise, nm, 0)
  sc <- .pv(p$scale, nm, 1)
  sh <- .pv(p$shift, nm, 0)
  if (h == 0 && sc == 1 && sh == 0) {
    return(if (int_valued) round(v) else v)
  }
  if (!is.finite(s)) s <- 0
  if (h > 0 && s > 0) v <- mu + (v - mu + h * s * stats::rnorm(length(v))) / sqrt(1 + h^2)
  if (sc != 1) v <- med + sc * (v - med)
  if (sh != 0) v <- v + sh * s
  if (p$bounds == "observed") {
    v <- pmin(pmax(v, lo), hi)
  } else if (p$bounds == "sign") {
    if (lo >= 0) v <- pmax(v, 0)
    if (hi <= 0) v <- pmin(v, 0)
  }
  # stochastic rounding keeps perturbed whole-number variables unbiased
  if (int_valued) v <- floor(v + stats::runif(length(v)))
  v
}

#' Calibrate the closeness of synthetic to real data
#'
#' @description
#' Finds the value of `closeness` (see [generate()]) for which synthetic data
#' reach a target value of a distance metric. The metric is evaluated on a
#' grid of closeness values with common random numbers, the grid interval
#' containing the target is located, and the solution is refined by
#' bisection, followed by a final linear interpolation step. Metric values
#' are themselves estimates, so the achieved value varies somewhat with the
#' random seed.
#'
#' Available metrics (all computed against `data`, which should be the data
#' the synthesizer was fitted on):
#'
#' * `"auc"` -- cross-validated discriminator AUC ([discriminator_auc()]);
#'   0.5 means indistinguishable, larger values mean more distant;
#' * `"dcr"` -- median distance to the closest real record relative to the
#'   real-to-real nearest-neighbour distance ([dcr()]); larger means further
#'   from real records;
#' * `"pmse"` -- pMSE ratio ([pmse()]) with the quadratic propensity model;
#' * `"marginal"` -- mean marginal distance (KS / TVD,
#'   [marginal_metrics()]);
#' * a **function** `function(real, synthetic)` returning a single number,
#'   for any other criterion, e.g. a disclosure-risk measure from
#'   [disclosure_risk()].
#'
#' @param object A fitted `sp_synthesizer`.
#' @param data The real data the synthesizer was fitted on.
#' @param target Target value of the metric.
#' @param metric One of `"auc"`, `"dcr"`, `"pmse"`, `"marginal"`, or a function
#'   `function(real, synthetic)` returning a number.
#' @param n Number of synthetic records generated per evaluation (default:
#'   size of `data`).
#' @param grid Closeness values evaluated first.
#' @param refine Number of bisection steps after the grid search.
#' @param max_rows Maximum rows used by the metric computations.
#' @param seed Random seed (the same seed is used at every closeness value).
#' @param ... Passed to [generate()] (e.g. `perturb` for directional
#'   perturbations kept fixed during calibration).
#' @return An object of class `sp_calibration` with elements `closeness`
#'   (the calibrated value), `achieved` (metric value at that closeness),
#'   `target`, `metric` and `path` (all evaluations). It has `print()` and
#'   `plot()` methods.
#' @examples
#' fit <- fit_synthesizer(iris)
#' cal <- calibrate_closeness(fit, iris, target = 0.75, metric = "auc",
#'                            grid = seq(0, 1, 0.25), refine = 2)
#' cal
#' syn <- generate(fit, closeness = cal$closeness, seed = 1)
#' @export
calibrate_closeness <- function(object, data, target,
                                metric = c("auc", "dcr", "pmse", "marginal"),
                                n = NULL, grid = seq(0, 1, by = 0.1), refine = 4L,
                                max_rows = 2000L, seed = 1, ...) {
  if (!inherits(object, "sp_synthesizer")) {
    stop("'object' must be an sp_synthesizer.", call. = FALSE)
  }
  metric_fun <- NULL
  if (is.function(metric)) {
    metric_fun <- metric
    metric <- "custom"
  } else {
    metric <- match.arg(metric)
  }
  if (!is.numeric(target) || length(target) != 1L || is.na(target)) {
    stop("'target' must be a single number.", call. = FALSE)
  }
  grid <- sort(unique(grid))
  if (length(grid) < 2L || any(grid < 0 | grid > 1)) {
    stop("'grid' must contain at least two values in [0, 1].", call. = FALSE)
  }
  refine <- .check_count(refine, "refine")
  n <- n %||% NROW(data)
  data <- if (is.data.frame(data)) data else as.data.frame(data)

  evaluate <- function(cl) {
    syn <- generate(object, n = n, closeness = cl, seed = seed, ...)
    if (!is.data.frame(syn)) syn <- as.data.frame(syn)
    if (object$input == "vector") names(syn) <- names(data)[1L]
    switch(metric,
      custom = {
        v <- metric_fun(data, syn)
        if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
          stop("A custom 'metric' must return a single number.", call. = FALSE)
        }
        v
      },
      auc = discriminator_auc(data, syn, max_rows = max_rows, seed = seed),
      dcr = dcr(data, syn, max_rows = max_rows, seed = seed)$ratio,
      pmse = pmse(data, syn, max_rows = max_rows, seed = seed, model = "quadratic")$pmse_ratio,
      marginal = mean(marginal_metrics(data, syn)$distance, na.rm = TRUE)
    )
  }
  vals <- vapply(grid, evaluate, numeric(1))
  path <- data.frame(closeness = grid, value = vals)

  dev <- vals - target
  cross <- which(dev[-1L] * dev[-length(dev)] <= 0)
  if (!length(cross)) {
    best <- which.min(abs(dev))
    warning(sprintf("Target %s = %g is outside the attainable range [%.3g, %.3g]; returning the closest grid value.",
                    metric, target, min(vals), max(vals)), call. = FALSE)
    cl <- grid[best]
    achieved <- vals[best]
  } else {
    i <- cross[1L]
    lo <- grid[i]
    hi <- grid[i + 1L]
    flo <- dev[i]
    fhi <- dev[i + 1L]
    for (k in seq_len(refine)) {
      mid <- (lo + hi) / 2
      fm <- evaluate(mid) - target
      path <- rbind(path, data.frame(closeness = mid, value = fm + target))
      if (fm * flo <= 0) {
        hi <- mid
        fhi <- fm
      } else {
        lo <- mid
        flo <- fm
      }
    }
    # final step: linear interpolation within the bracket
    if (refine > 0L && fhi != flo) {
      xi <- lo + (hi - lo) * (0 - flo) / (fhi - flo)
      if (xi > lo && xi < hi) {
        path <- rbind(path, data.frame(closeness = xi, value = evaluate(xi)))
      }
    }
    path <- path[order(path$closeness), , drop = FALSE]
    best <- which.min(abs(path$value - target))
    cl <- path$closeness[best]
    achieved <- path$value[best]
  }
  rownames(path) <- NULL
  structure(list(closeness = cl, achieved = achieved, target = target,
                 metric = metric, path = path),
            class = "sp_calibration")
}

#' @export
print.sp_calibration <- function(x, ...) {
  cat("<sp_calibration>\n")
  cat(sprintf("  metric    : %s (target %g)\n", x$metric, x$target))
  cat(sprintf("  closeness : %.3f  (achieved %s = %.3f)\n", x$closeness, x$metric, x$achieved))
  cat("  use: generate(fit, closeness = <value>)\n")
  invisible(x)
}

#' @export
plot.sp_calibration <- function(x, ...) {
  lab <- c(auc = "discriminator AUC", dcr = "DCR ratio", pmse = "pMSE ratio",
           marginal = "mean marginal distance", custom = "metric")[[x$metric]]
  p <- ggplot2::ggplot(x$path, ggplot2::aes(x = .data$closeness, y = .data$value)) +
    ggplot2::geom_hline(yintercept = x$target, colour = "#898781", linetype = "dashed") +
    ggplot2::geom_line(colour = .sp_cols[["Synthetic"]], linewidth = 0.7) +
    ggplot2::geom_point(colour = .sp_cols[["Synthetic"]], size = 2) +
    ggplot2::annotate("point", x = x$closeness, y = x$achieved, size = 4, shape = 21,
                      fill = "white", colour = "#0b0b0b") +
    ggplot2::labs(x = "closeness", y = lab,
                  title = sprintf("Calibration: %s = %g at closeness %.2f", lab,
                                  x$target, x$closeness)) +
    .theme_sp()
  print(p)
  invisible(p)
}
