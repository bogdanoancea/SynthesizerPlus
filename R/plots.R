# Visualisation ----------------------------------------------------------------

# Colour-vision-deficiency-checked palette: real = blue, synthetic = orange
.sp_cols <- c(Real = "#2a78d6", Synthetic = "#eb6834")
.sp_div <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")

.theme_sp <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "#52514e"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#0b0b0b"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

.long_source <- function(real, synthetic, vars, as_numeric = TRUE) {
  parts <- list()
  srcs <- list(Real = real, Synthetic = synthetic)
  for (s in names(srcs)) {
    df <- srcs[[s]]
    if (is.null(df)) next
    for (v in vars) {
      val <- df[[v]]
      val <- if (as_numeric) as.numeric(val) else as.character(val)
      parts[[length(parts) + 1L]] <- data.frame(variable = v, value = val, source = s,
                                                stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, parts)
  out$variable <- factor(out$variable, levels = vars)
  out$source <- factor(out$source, levels = c("Real", "Synthetic"))
  out
}

.as_df_pair <- function(real, synthetic) {
  to_df <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "ts")) x <- as.matrix(x)
    if (is.matrix(x)) x <- as.data.frame(x)
    if (is.atomic(x) && is.null(dim(x))) x <- data.frame(x = x)
    if (!is.data.frame(x)) stop("Data must be data frames.", call. = FALSE)
    x
  }
  list(real = to_df(real), synthetic = to_df(synthetic))
}

.pick_vars <- function(df, vars, numeric, max_vars) {
  if (is.null(vars)) {
    keep <- vapply(df, function(v) if (numeric) .is_numeric_like(v) else !.is_numeric_like(v),
                   logical(1))
    vars <- names(df)[keep]
    if (length(vars) > max_vars) vars <- vars[seq_len(max_vars)]
  } else {
    .check_cols(vars, names(df), "vars")
  }
  vars
}

#' Plot marginal distributions of real and synthetic data
#'
#' `plot_marginals()` overlays the distributions of numeric variables
#' (density, histogram or empirical CDF) and shows side-by-side category
#' proportions for categorical variables. When the selected variables
#' contain both kinds, a list of two plots (class `sp_plot_list`) is
#' returned; printing it draws both.
#'
#' @param real A data frame (or vector/matrix) with the real data.
#' @param synthetic Optional synthetic data with the same columns.
#' @param vars Variables to plot (default: all, up to `max_vars` per kind).
#' @param kind For numeric variables: `"density"`, `"histogram"` or `"ecdf"`.
#' @param max_vars Maximum number of variables of each kind shown by default.
#' @param max_levels Maximum categories per categorical variable (the rest
#'   are pooled into "(other)").
#' @param bins Number of histogram bins.
#' @return A `ggplot` object or an `sp_plot_list`.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' plot_marginals(iris, syn)
#' plot_marginals(iris, syn, vars = "Sepal.Length", kind = "ecdf")
#' @export
plot_marginals <- function(real, synthetic = NULL, vars = NULL,
                           kind = c("density", "histogram", "ecdf"),
                           max_vars = 12L, max_levels = 15L, bins = 30L) {
  kind <- match.arg(kind)
  d <- .as_df_pair(real, synthetic)
  if (is.null(vars)) {
    nv <- .pick_vars(d$real, NULL, TRUE, max_vars)
    cv <- .pick_vars(d$real, NULL, FALSE, max_vars)
  } else {
    .check_cols(vars, names(d$real), "vars")
    isnum <- vapply(d$real[vars], .is_numeric_like, logical(1))
    nv <- vars[isnum]
    cv <- vars[!isnum]
  }
  plots <- list()
  if (length(nv)) plots$numeric <- .plot_numeric(d, nv, kind, bins)
  if (length(cv)) plots$categorical <- plot_categorical(d$real, d$synthetic, cv,
                                                         max_levels = max_levels)
  if (!length(plots)) stop("No variables to plot.", call. = FALSE)
  if (length(plots) == 1L) return(plots[[1L]])
  structure(plots, class = "sp_plot_list")
}

.plot_numeric <- function(d, vars, kind, bins) {
  long <- .long_source(d$real, d$synthetic, vars)
  long <- long[!is.na(long$value), , drop = FALSE]
  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, colour = .data$source,
                                          fill = .data$source))
  g <- switch(kind,
    density = g + ggplot2::geom_density(alpha = 0.15, linewidth = 0.7),
    histogram = g + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(.data$density)), bins = bins,
      position = "identity", alpha = 0.35, linewidth = 0.2),
    ecdf = g + ggplot2::stat_ecdf(linewidth = 0.7, geom = "step")
  )
  g + ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::scale_colour_manual(values = .sp_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = .sp_cols, drop = FALSE) +
    ggplot2::labs(x = NULL, y = if (kind == "ecdf") "cumulative probability" else "density",
                  title = "Marginal distributions") +
    .theme_sp()
}

#' @rdname plot_marginals
#' @export
plot_categorical <- function(real, synthetic = NULL, vars = NULL, max_levels = 15L) {
  d <- .as_df_pair(real, synthetic)
  vars <- .pick_vars(d$real, vars, FALSE, 12L)
  if (!length(vars)) stop("No categorical variables to plot.", call. = FALSE)
  long <- .long_source(d$real, d$synthetic, vars, as_numeric = FALSE)
  long$value[is.na(long$value)] <- "(missing)"
  # lump rare levels per variable based on the real data
  for (v in vars) {
    rows <- long$variable == v
    tab <- sort(table(long$value[rows & long$source == "Real"]), decreasing = TRUE)
    if (length(tab) > max_levels) {
      keep <- names(tab)[seq_len(max_levels - 1L)]
      long$value[rows & !long$value %in% keep] <- "(other)"
    }
  }
  agg <- stats::aggregate(list(count = rep(1, nrow(long))),
                          by = list(variable = long$variable, value = long$value,
                                    source = long$source), FUN = sum)
  tot <- stats::ave(agg$count, agg$variable, agg$source, FUN = sum)
  agg$proportion <- agg$count / tot
  ggplot2::ggplot(agg, ggplot2::aes(x = .data$value, y = .data$proportion,
                                    fill = .data$source)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(padding = 0.15, preserve = "single"),
                      width = 0.8) +
    ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::scale_fill_manual(values = .sp_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(z) paste0(round(100 * z), "%")) +
    ggplot2::labs(x = NULL, y = "share of records", title = "Categorical distributions") +
    .theme_sp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' @export
print.sp_plot_list <- function(x, ...) {
  for (p in x) print(p)
  invisible(x)
}

#' Plot association matrices of real and synthetic data
#'
#' Heatmaps of the mixed-type [association_matrix()] of the real data, the
#' synthetic data and their difference (synthetic minus real), on a common
#' diverging scale.
#'
#' @inheritParams plot_marginals
#' @param show_values Print the values in the cells? (Default: when there
#'   are at most 8 variables.)
#' @return A `ggplot` object.
#' @examples
#' syn <- synthesize(mtcars, seed = 1)
#' plot_association(mtcars, syn, vars = c("mpg", "cyl", "disp", "hp", "wt"))
#' @export
plot_association <- function(real, synthetic = NULL, vars = NULL, show_values = NULL) {
  d <- .as_df_pair(real, synthetic)
  vars <- vars %||% names(d$real)
  .check_cols(vars, names(d$real), "vars")
  if (length(vars) < 2L) stop("At least two variables are needed.", call. = FALSE)
  show_values <- show_values %||% (length(vars) <= 8L)
  mats <- list(Real = association_matrix(d$real, vars))
  if (!is.null(d$synthetic)) {
    mats$Synthetic <- association_matrix(d$synthetic, vars)
    mats$`Difference (synthetic - real)` <- mats$Synthetic - mats$Real
  }
  long <- do.call(rbind, lapply(names(mats), function(nm) {
    M <- mats[[nm]]
    data.frame(row = factor(rep(vars, times = length(vars)), levels = rev(vars)),
               col = factor(rep(vars, each = length(vars)), levels = vars),
               value = as.vector(M), panel = nm, stringsAsFactors = FALSE)
  }))
  long$panel <- factor(long$panel, levels = names(mats))
  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::facet_wrap(~panel, nrow = 1L) +
    ggplot2::scale_fill_gradient2(low = .sp_div[["low"]], mid = .sp_div[["mid"]],
                                  high = .sp_div[["high"]], midpoint = 0,
                                  limits = c(-1, 1), na.value = "grey90",
                                  name = "association") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL, title = "Pairwise association") +
    .theme_sp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.title = ggplot2::element_text(),
                   panel.grid = ggplot2::element_blank())
  if (show_values) {
    g <- g + ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(.data$value), "",
                                                            sprintf("%.2f", .data$value))),
                                size = 2.6, colour = "#0b0b0b")
  }
  g
}

#' Scatter-plot matrix of real and/or synthetic data
#'
#' Pairwise scatter plots of numeric variables, with real and synthetic
#' records side by side (columns of the panel grid are variables, the colour
#' identifies the source). Useful to inspect multivariate structure, also for
#' data generated with [r_mvdist()] or [mvdist] functions.
#'
#' @inheritParams plot_marginals
#' @param max_points Maximum number of points drawn per data set.
#' @param alpha Point transparency.
#' @return A `ggplot` object.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' plot_pairs(iris, syn, vars = c("Sepal.Length", "Petal.Length", "Petal.Width"))
#'
#' x <- as.data.frame(r_copula(1000, "clayton", dim = 3, theta = 3, seed = 1))
#' plot_pairs(x)
#' @export
plot_pairs <- function(real, synthetic = NULL, vars = NULL, max_points = 1500L,
                       alpha = 0.35) {
  d <- .as_df_pair(real, synthetic)
  vars <- .pick_vars(d$real, vars, TRUE, 5L)
  if (length(vars) < 2L) stop("At least two numeric variables are needed.", call. = FALSE)
  srcs <- list(Real = d$real, Synthetic = d$synthetic)
  parts <- list()
  for (s in names(srcs)) {
    df <- srcs[[s]]
    if (is.null(df)) next
    df <- .subsample(df[vars], max_points)
    for (i in vars) for (j in vars) {
      if (i == j) next
      parts[[length(parts) + 1L]] <- data.frame(
        x = as.numeric(df[[i]]), y = as.numeric(df[[j]]),
        xvar = i, yvar = j, source = s, stringsAsFactors = FALSE)
    }
  }
  long <- do.call(rbind, parts)
  long$xvar <- factor(long$xvar, levels = vars)
  long$yvar <- factor(long$yvar, levels = vars)
  long$source <- factor(long$source, levels = c("Real", "Synthetic"))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$source)) +
    ggplot2::geom_point(alpha = alpha, size = 0.8, na.rm = TRUE) +
    ggplot2::facet_grid(yvar ~ xvar, scales = "free") +
    ggplot2::scale_colour_manual(values = .sp_cols, drop = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggplot2::labs(x = NULL, y = NULL, title = "Pairwise scatter plots") +
    .theme_sp()
}

#' Quantile-quantile plots of synthetic against real data
#'
#' For each numeric variable, quantiles of the synthetic data are plotted
#' against those of the real data; points on the diagonal indicate matching
#' distributions.
#'
#' @inheritParams plot_marginals
#' @param probs Probabilities at which quantiles are compared.
#' @return A `ggplot` object.
#' @examples
#' syn <- synthesize(mtcars, seed = 1)
#' plot_qq(mtcars, syn, vars = c("mpg", "hp", "wt"))
#' @export
plot_qq <- function(real, synthetic, vars = NULL, max_vars = 12L,
                    probs = seq(0.01, 0.99, by = 0.01)) {
  d <- .as_df_pair(real, synthetic)
  vars <- .pick_vars(d$real, vars, TRUE, max_vars)
  if (!length(vars)) stop("No numeric variables to plot.", call. = FALSE)
  long <- do.call(rbind, lapply(vars, function(v) {
    data.frame(variable = v,
               real = stats::quantile(as.numeric(d$real[[v]]), probs, na.rm = TRUE, names = FALSE),
               synthetic = stats::quantile(as.numeric(d$synthetic[[v]]), probs, na.rm = TRUE,
                                           names = FALSE),
               stringsAsFactors = FALSE)
  }))
  long$variable <- factor(long$variable, levels = vars)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$real, y = .data$synthetic)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#898781", linetype = "dashed") +
    ggplot2::geom_point(colour = .sp_cols[["Synthetic"]], size = 1.2) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::labs(x = "real quantiles", y = "synthetic quantiles", title = "Q-Q plots") +
    .theme_sp()
}

#' Plot real and synthetic time series
#'
#' `plot_ts()` draws the real and synthetic series one above the other for
#' each component; `plot_acf()` compares their autocorrelation functions.
#'
#' @param real A `ts` object (or numeric vector/matrix).
#' @param synthetic A synthetic `ts` with the same number of series.
#' @param lag_max Maximum lag for the autocorrelation function.
#' @return A `ggplot` object.
#' @examples
#' fit <- fit_synthesizer(ldeaths)
#' syn <- generate(fit, seed = 1)
#' plot_ts(ldeaths, syn)
#' plot_acf(ldeaths, syn)
#' @export
plot_ts <- function(real, synthetic = NULL) {
  to_long <- function(x, src) {
    if (is.null(x)) return(NULL)
    tt <- if (inherits(x, "ts")) as.numeric(stats::time(x)) else seq_len(NROW(x))
    M <- as.matrix(x)
    nm <- colnames(M) %||% if (ncol(M) == 1L) "Series" else paste0("Series", seq_len(ncol(M)))
    data.frame(time = rep(tt, ncol(M)), value = as.vector(M),
               series = rep(nm, each = nrow(M)), source = src,
               stringsAsFactors = FALSE)
  }
  long <- rbind(to_long(real, "Real"), to_long(synthetic, "Synthetic"))
  long$source <- factor(long$source, levels = c("Real", "Synthetic"))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$time, y = .data$value, colour = .data$source)) +
    ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggplot2::facet_grid(source ~ series, scales = "free_y") +
    ggplot2::scale_colour_manual(values = .sp_cols, drop = TRUE) +
    ggplot2::labs(x = "time", y = NULL, title = "Time series") +
    .theme_sp()
}

#' @rdname plot_ts
#' @export
plot_acf <- function(real, synthetic, lag_max = 24L) {
  acf_df <- function(x, src) {
    M <- as.matrix(x)
    nm <- colnames(M) %||% if (ncol(M) == 1L) "Series" else paste0("Series", seq_len(ncol(M)))
    do.call(rbind, lapply(seq_len(ncol(M)), function(j) {
      a <- stats::acf(M[, j], lag.max = lag_max, plot = FALSE, na.action = stats::na.pass)
      data.frame(lag = as.vector(a$lag)[-1L], acf = as.vector(a$acf)[-1L],
                 series = nm[j], source = src, stringsAsFactors = FALSE)
    }))
  }
  long <- rbind(acf_df(real, "Real"), acf_df(synthetic, "Synthetic"))
  long$source <- factor(long$source, levels = c("Real", "Synthetic"))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$lag, y = .data$acf, fill = .data$source)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#c3c2b7") +
    ggplot2::geom_col(position = ggplot2::position_dodge2(padding = 0.1), width = 0.8) +
    ggplot2::facet_wrap(~series, scales = "free_x") +
    ggplot2::scale_fill_manual(values = .sp_cols) +
    ggplot2::labs(x = "lag", y = "autocorrelation", title = "Autocorrelation functions") +
    .theme_sp()
}

#' Plot distances to the closest record
#'
#' Compares the distribution of distances from synthetic records to their
#' closest real record with the real-to-real baseline (see [dcr()]).
#'
#' @param x The result of [dcr()] or [compare_synthetic()].
#' @return A `ggplot` object.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' plot_dcr(dcr(iris, syn, seed = 1))
#' @export
plot_dcr <- function(x) {
  if (inherits(x, "sp_comparison")) x <- x$privacy
  if (is.null(x$synthetic) || is.null(x$real)) {
    stop("'x' must be the result of dcr() or compare_synthetic().", call. = FALSE)
  }
  long <- data.frame(distance = c(x$real, x$synthetic),
                     source = factor(rep(c("Real", "Synthetic"),
                                         c(length(x$real), length(x$synthetic))),
                                     levels = c("Real", "Synthetic")))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$distance, colour = .data$source,
                                     fill = .data$source)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.7, bounds = c(0, Inf)) +
    ggplot2::scale_colour_manual(values = .sp_cols) +
    ggplot2::scale_fill_manual(values = .sp_cols) +
    ggplot2::labs(x = "Gower distance to closest real record", y = "density",
                  title = "Distance to closest record",
                  subtitle = "Real: distance to the nearest *other* real record") +
    .theme_sp()
}

#' Plot a comparison of real and synthetic data
#'
#' @param object,x An `sp_comparison` from [compare_synthetic()].
#' @param type What to plot: `"marginals"`, `"association"`, `"qq"`,
#'   `"pairs"` or `"dcr"`.
#' @param ... Passed to the underlying plot function.
#' @return `autoplot()` returns a `ggplot` (or `sp_plot_list`); `plot()`
#'   prints it and returns it invisibly.
#' @examples
#' syn <- synthesize(iris, seed = 1)
#' cmp <- compare_synthetic(iris, syn, seed = 1)
#' plot(cmp, type = "qq")
#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.sp_comparison <- function(object,
                                   type = c("marginals", "association", "qq", "pairs", "dcr"),
                                   ...) {
  type <- match.arg(type)
  switch(type,
    marginals = plot_marginals(object$real, object$synthetic, ...),
    association = plot_association(object$real, object$synthetic, ...),
    qq = plot_qq(object$real, object$synthetic, ...),
    pairs = plot_pairs(object$real, object$synthetic, ...),
    dcr = {
      if (is.null(object$privacy)) stop("No privacy metrics in this comparison.", call. = FALSE)
      plot_dcr(object$privacy)
    }
  )
}

#' @rdname autoplot.sp_comparison
#' @export
plot.sp_comparison <- function(x, type = c("marginals", "association", "qq", "pairs", "dcr"),
                               ...) {
  p <- autoplot.sp_comparison(x, type = type, ...)
  print(p)
  invisible(p)
}

#' @export
ggplot2::autoplot
