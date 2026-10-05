# Data augmentation for machine learning -------------------------------------

#' Augment a small training set with synthetic data
#'
#' @description
#' Fits a synthesizer to `data` and returns the real records together with
#' synthetic ones, ready to train a machine-learning model.
#'
#' * For a data frame (or matrix), `n` synthetic rows are appended to the
#'   real rows. An indicator column can be added with `label`. A vector is
#'   extended with `n` synthetic values.
#' * For a time series, `nsim` synthetic series of length `n` are drawn
#'   (see [fit_synthesizer.ts()]) and returned in a list, preceded by the real
#'   series. Use [ts_windows()] to turn the list into a table of lagged
#'   features.
#'
#' Synthetic data cannot add information that is not in the real data. What
#' they add is a smoothed version of the training distribution, which acts as
#' a regulariser and mainly helps flexible, high-variance learners (trees,
#' nearest neighbours, neural networks) trained on small samples. Always
#' **fit the synthesizer on the training data only** and measure the effect
#' on held-out *real* data; see `vignette("augmentation")`.
#'
#' @param data Training data: a data frame, matrix or `ts` object.
#' @param n For tabular data: number of synthetic rows (default: 10 times the
#'   training size). For time series: length of each synthetic series
#'   (default: the length of `data`).
#' @param nsim For time series: number of synthetic series.
#' @param keep_real Include the real data in the result?
#' @param label Optional name of a logical column marking synthetic rows
#'   (tabular data only).
#' @param closeness,perturb How close the synthetic rows are to the real
#'   ones (tabular data only); see [generate()]. Lower closeness gives more
#'   diverse synthetic rows.
#' @param seed Optional random seed.
#' @param ... Passed to [fit_synthesizer()] (e.g. `by`, `copula`, or for
#'   time series `method`, `order`).
#' @return A data frame (tabular data) or a list of `ts` objects (time series).
#' @seealso [ts_windows()], [fit_synthesizer()]
#' @examples
#' train <- mtcars[1:20, ]
#' big <- augment_data(train, n = 500, label = "synthetic", seed = 1)
#' table(big$synthetic)
#'
#' series <- augment_data(ldeaths, nsim = 5, order = 12, seed = 1)
#' length(series)
#' head(ts_windows(series, lags = 12))
#' @export
augment_data <- function(data, ...) {
  UseMethod("augment_data")
}

#' @rdname augment_data
#' @export
augment_data.data.frame <- function(data, n = 10L * nrow(data), keep_real = TRUE,
                                    label = NULL, closeness = 1, perturb = NULL,
                                    seed = NULL, ...) {
  .check_flag(keep_real, "keep_real")
  if (!is.null(label) && (!is.character(label) || length(label) != 1L || label %in% names(data))) {
    stop("'label' must be a single new column name.", call. = FALSE)
  }
  syn <- .with_seed(seed, generate(fit_synthesizer(data, ...), n = n,
                                   closeness = closeness, perturb = perturb))
  out <- if (keep_real) rbind(.as_plain_df(data), syn) else syn
  if (!is.null(label)) {
    out[[label]] <- c(rep(FALSE, if (keep_real) nrow(data) else 0L), rep(TRUE, nrow(syn)))
  }
  attr(out, "row.names") <- .set_row_names(nrow(out))
  out
}

#' @rdname augment_data
#' @export
augment_data.default <- function(data, ...) {
  if (!is.atomic(data) || !is.null(dim(data)) || is.na(.value_class(data))) {
    stop("Cannot augment an object of class '", paste(class(data), collapse = "/"), "'.",
         call. = FALSE)
  }
  augment_data.data.frame(data.frame(x = data, stringsAsFactors = FALSE), ...)$x
}

#' @rdname augment_data
#' @export
augment_data.matrix <- function(data, ...) {
  as.matrix(augment_data.data.frame(as.data.frame(data, stringsAsFactors = FALSE), ...))
}

#' @rdname augment_data
#' @export
augment_data.ts <- function(data, nsim = 20L, n = NROW(data), keep_real = TRUE,
                            seed = NULL, ...) {
  .check_flag(keep_real, "keep_real")
  dots <- list(...)
  if (any(c("closeness", "perturb", "label") %in% names(dots))) {
    warning("'closeness', 'perturb' and 'label' apply to tabular data only and were ignored.",
            call. = FALSE)
    dots <- dots[setdiff(names(dots), c("closeness", "perturb", "label"))]
  }
  nsim <- .check_count(nsim, "nsim", allow_zero = FALSE)
  sims <- .with_seed(seed, {
    fit <- do.call(fit_synthesizer, c(list(data), dots))
    lapply(seq_len(nsim), function(i) generate(fit, n = n))
  })
  if (keep_real) c(list(data), sims) else sims
}

#' Build a lagged-feature table from one or more time series
#'
#' Converts time series into a data frame suitable for training forecasting
#' models: each row holds the target `y` (the value `horizon` steps ahead) and
#' the `lags` preceding values `lag1` (most recent) ... `lagp`. With a list of
#' series (e.g. from [augment_data()]), the windows of all series are stacked
#' and a `series` column identifies their origin; windows never cross series
#' boundaries.
#'
#' For multivariate series, lags of every component are created
#' (`<name>_lag1`, ...) and the target is the column named by `target`.
#'
#' @param x A `ts` object, numeric vector/matrix, or a list of them.
#' @param lags Number of lags.
#' @param horizon Forecast horizon (1 = one step ahead).
#' @param target For multivariate series: name or index of the target column.
#' @return A data frame with columns `y`, the lag columns and (for lists)
#'   `series`.
#' @examples
#' w <- ts_windows(ldeaths, lags = 3)
#' head(w)
#' fit <- lm(y ~ ., data = w)
#' @export
ts_windows <- function(x, lags, horizon = 1L, target = 1L) {
  lags <- .check_count(lags, "lags", allow_zero = FALSE)
  horizon <- .check_count(horizon, "horizon", allow_zero = FALSE)
  if (is.list(x) && !is.data.frame(x)) {
    parts <- lapply(seq_along(x), function(i) {
      w <- ts_windows(x[[i]], lags = lags, horizon = horizon, target = target)
      w$series <- rep(i, nrow(w))
      w
    })
    out <- do.call(rbind, parts)
    attr(out, "row.names") <- .set_row_names(nrow(out))
    return(out)
  }
  M <- as.matrix(x)
  if (!is.numeric(M)) stop("Series must be numeric.", call. = FALSE)
  N <- nrow(M)
  if (N < lags + horizon) {
    stop(sprintf("A series of length %d is too short for %d lags and horizon %d.",
                 N, lags, horizon), call. = FALSE)
  }
  if (is.character(target)) {
    if (!target %in% colnames(M)) stop("Unknown target column: ", target, call. = FALSE)
    target <- match(target, colnames(M))
  }
  if (target < 1L || target > ncol(M)) stop("Invalid 'target'.", call. = FALSE)
  rows <- (lags + horizon):N
  out <- list(y = M[rows, target])
  nm <- colnames(M)
  for (j in seq_len(ncol(M))) {
    for (k in seq_len(lags)) {
      col <- if (ncol(M) == 1L) paste0("lag", k) else paste0(nm[j] %||% paste0("V", j), "_lag", k)
      out[[col]] <- M[rows - horizon - k + 1L, j]
    }
  }
  out <- as.data.frame(out, optional = TRUE)
  attr(out, "row.names") <- .set_row_names(length(rows))
  out
}
