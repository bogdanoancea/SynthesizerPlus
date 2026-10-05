# Marginal models --------------------------------------------------------------
#
# Each column of the training data is described by an `sp_marginal` object:
#   type = "continuous"  interpolated empirical quantile function (linear /
#                        monotone spline between order statistics)
#          "discrete"    finite numeric support with probabilities
#          "categorical" finite set of labels with probabilities
#          "empty"       column is entirely missing
# Values are mapped to the unit interval (pseudo-observations) for fitting the
# copula, and uniforms are mapped back through the quantile function when
# generating.

.fit_marginal <- function(x, name, interpolation = "linear",
                          discrete_threshold = 20L, knots = NULL) {
  cls <- .value_class(x)
  if (is.na(cls)) {
    stop(sprintf("Column '%s' has unsupported class '%s'.", name,
                 paste(class(x), collapse = "/")), call. = FALSE)
  }
  miss <- is.na(x)
  out <- list(
    name = name,
    class = cls,
    na_rate = mean(miss),
    n = length(x),
    tz = if (cls == "POSIXct") (attr(x, "tzone") %||% "") else NULL
  )
  xo <- x[!miss]

  if (length(xo) == 0L) {
    out$type <- "empty"
    return(structure(out, class = "sp_marginal"))
  }

  if (.is_categorical_class(cls)) {
    lev <- switch(cls,
      factor = , ordered = levels(x),
      logical = c("FALSE", "TRUE"),
      character = sort(unique(xo))
    )
    counts <- tabulate(match(as.character(xo), lev), nbins = length(lev))
    out$type <- "categorical"
    out$values <- lev
    out$probs <- counts / sum(counts)
    return(structure(out, class = "sp_marginal"))
  }

  num <- as.numeric(xo)
  uv <- sort(unique(num))
  out$integer_valued <- all(num == round(num))

  repeats <- length(uv) <= length(num) / 2
  if (length(uv) <= discrete_threshold && (out$integer_valued || repeats)) {
    counts <- tabulate(match(num, uv), nbins = length(uv))
    out$type <- "discrete"
    out$values <- uv
    out$probs <- counts / sum(counts)
    return(structure(out, class = "sp_marginal"))
  }

  ys <- sort(num)
  if (!is.null(knots) && length(ys) > knots) {
    # Compress the empirical quantile function to `knots` quantiles. This keeps
    # the model small and avoids storing every original value.
    pk <- (seq_len(knots) - 0.5) / knots
    ys <- as.numeric(stats::quantile(ys, probs = pk, type = 8, names = FALSE))
    ys[1L] <- min(num)
    ys[knots] <- max(num)
  }
  out$type <- "continuous"
  out$quantiles <- ys
  out$interpolation <- interpolation
  structure(out, class = "sp_marginal")
}

# Map observed values to pseudo-observations on (0, 1).
# Continuous: randomised ranks / (n + 1). Discrete / categorical: the
# distributional transform F(x-) + V * P(X = x), V ~ U(0, 1), which makes the
# pseudo-observations uniform and lets one Gaussian copula handle mixed data.
.pseudo_obs <- function(m, x) {
  u <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  if (!any(ok) || m$type == "empty") return(u)
  if (m$type == "continuous") {
    r <- rank(as.numeric(x[ok]), ties.method = "random")
    u[ok] <- r / (sum(ok) + 1)
    return(u)
  }
  key <- if (m$type == "categorical") as.character(x[ok]) else as.numeric(x[ok])
  idx <- match(key, m$values)
  cum <- c(0, cumsum(m$probs))
  lo <- cum[idx]
  width <- m$probs[idx]
  u[ok] <- lo + stats::runif(length(idx)) * width
  # keep strictly inside (0, 1)
  u[ok] <- pmin(pmax(u[ok], 1e-10), 1 - 1e-10)
  u
}

# Map uniforms to values (numeric representation or labels)
.marginal_quantile <- function(m, u) {
  n <- length(u)
  if (m$type == "empty") return(rep(NA, n))
  if (m$type == "continuous") {
    ys <- m$quantiles
    k <- length(ys)
    if (k == 1L) return(rep(ys, n))
    p <- (seq_len(k) - 0.5) / k
    pp <- p[1L] + u * (p[k] - p[1L])
    v <- if (identical(m$interpolation, "spline")) {
      stats::splinefun(p, ys, method = "hyman")(pp)
    } else {
      stats::approx(p, ys, xout = pp, ties = "ordered")$y
    }
    v <- pmin(pmax(v, ys[1L]), ys[k])
    if (isTRUE(m$integer_valued)) v <- round(v)
    return(v)
  }
  cum <- cumsum(m$probs)
  cum[length(cum)] <- 1
  idx <- findInterval(u, cum, left.open = TRUE) + 1L
  idx <- pmin(idx, length(m$values))
  m$values[idx]
}

# Convert generated values back to the class of the template column
.restore_class <- function(v, m, template_col) {
  cls <- m$class
  if (m$type == "empty") {
    return(template_col[rep(NA_integer_, length(v))])
  }
  switch(cls,
    factor = factor(v, levels = levels(template_col)),
    ordered = factor(v, levels = levels(template_col), ordered = TRUE),
    character = as.character(v),
    logical = as.logical(v),
    integer = as.integer(round(as.numeric(v))),
    numeric = as.numeric(v),
    numeric_classed = .restore_numeric_classed(v, template_col),
    Date = as.Date(round(as.numeric(v)), origin = "1970-01-01"),
    POSIXct = as.POSIXct(as.numeric(v), origin = "1970-01-01", tz = m$tz %||% ""),
    v
  )
}

# Indicator marginal for missingness of a column (1 = observed, 2 = missing)
.missing_marginal <- function(m) {
  structure(list(
    name = paste0(".missing_", m$name),
    class = "logical",
    na_rate = 0,
    type = "categorical",
    values = c("observed", "missing"),
    probs = c(1 - m$na_rate, m$na_rate)
  ), class = "sp_marginal")
}

#' @export
print.sp_marginal <- function(x, ...) {
  cat(sprintf("<sp_marginal> %s [%s, %s]", x$name, x$class, x$type))
  if (x$na_rate > 0) cat(sprintf(", %s missing", .fmt_pct(x$na_rate)))
  cat("\n")
  invisible(x)
}
