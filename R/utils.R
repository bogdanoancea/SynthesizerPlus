# Internal utilities ---------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Evaluate `code` with a temporary RNG seed, restoring the caller's RNG state
# afterwards (including the case where no `.Random.seed` existed before).
.with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be NULL or a single number.", call. = FALSE)
  }
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.check_count <- function(n, name = "n", allow_zero = TRUE) {
  ok <- is.numeric(n) && length(n) == 1L && !is.na(n) && is.finite(n) &&
    n == round(n) && (n > 0 || (allow_zero && n == 0))
  if (!ok) {
    stop(sprintf("'%s' must be a single %s integer.", name,
                 if (allow_zero) "non-negative" else "positive"), call. = FALSE)
  }
  as.integer(n)
}

# Warn about arguments passed through `...` that nothing uses (typos)
.warn_unused_dots <- function(dots, fun, hint = NULL) {
  if (!length(dots)) return(invisible(NULL))
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  nm[nm == ""] <- "<unnamed>"
  warning(sprintf("Unused argument(s) in %s(): %s.%s", fun, paste(nm, collapse = ", "),
                  if (is.null(hint)) "" else paste0(" ", hint)), call. = FALSE)
  invisible(NULL)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop(sprintf("'%s' must be a single number in [0, 1].", name), call. = FALSE)
  }
  x
}

# Classify the storage class of a column into the classes this package can
# reproduce. Returns NA for unsupported classes.
.value_class <- function(x) {
  if (is.ordered(x)) return("ordered")
  if (is.factor(x)) return("factor")
  if (inherits(x, "POSIXct")) return("POSIXct")
  if (inherits(x, "Date")) return("Date")
  if (inherits(x, "integer64")) return(NA_character_)   # 64-bit integers: not supported
  if ((is.double(x) || is.integer(x)) && !is.null(oldClass(x))) {
    return("numeric_classed")                           # e.g. difftime, hms, units
  }
  if (is.logical(x)) return("logical")
  if (is.character(x)) return("character")
  if (is.integer(x)) return("integer")
  if (is.double(x)) return("numeric")
  NA_character_
}

.is_categorical_class <- function(cls) {
  cls %in% c("factor", "ordered", "character", "logical")
}

# Numeric columns (incl. dates) that are treated as continuous in metrics/plots
.is_numeric_like <- function(x) {
  !is.factor(x) && !is.logical(x) &&
    (is.numeric(x) || is.double(x) || is.integer(x))
}

# Restore a classed numeric vector (difftime, hms, ...) from a template column
.restore_numeric_classed <- function(v, template) {
  v <- as.numeric(v)
  a <- attributes(template)
  a <- a[setdiff(names(a), c("names", "dim", "dimnames"))]
  attributes(v) <- a
  v
}

.as_plain_df <- function(x) {
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
    attr(x, "row.names") <- .set_row_names(nrow(x))
    return(x)
  }
  stop("Expected a data frame.", call. = FALSE)
}

.set_row_names <- function(n) {
  if (n == 0L) integer(0) else c(NA_integer_, -as.integer(n))
}

.require <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required %s. Install it with install.packages(\"%s\").",
                 pkg, why, pkg), call. = FALSE)
  }
  invisible(TRUE)
}

# Subsample rows of a data frame to at most `max_rows`
.subsample <- function(x, max_rows) {
  if (is.null(max_rows) || nrow(x) <= max_rows) return(x)
  x[sort(sample.int(nrow(x), max_rows)), , drop = FALSE]
}

# Empty data frame with the same columns/classes as `template`
.prototype <- function(template) {
  template[0L, , drop = FALSE]
}

.fmt_pct <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}
