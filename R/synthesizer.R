# Synthesizer: fitting and generation ----------------------------------------

#' Fit a synthetic-data model
#'
#' @description
#' `fit_synthesizer()` learns a generative model from a data set. The model
#' combines
#'
#' * a **marginal model per column**: an empirical quantile function for
#'   continuous variables (linear or monotone-spline interpolation), and an
#'   empirical probability table for discrete, categorical, ordinal and
#'   logical variables;
#' * a **copula** (Gaussian, Student-*t* or independence) that couples the
#'   columns. Mixed continuous / discrete data are handled with the
#'   distributional transform, so categorical variables take part in the
#'   dependence structure;
#' * an optional **missing-data model**: for each column with missing values
#'   a missingness indicator is added to the copula, so that the rate of
#'   missingness *and* its association with other variables are reproduced;
#' * optional **stratification** (`by`), fitting separate models within
#'   groups, with a pooled fallback for small groups.
#'
#' Use [generate()] to draw synthetic data from the fitted model, or
#' [synthesize()] to fit and generate in one step.
#'
#' @param data A data frame, matrix, atomic vector (numeric, integer,
#'   logical, character, factor, `Date`, `POSIXct`) or `ts` object. Time
#'   series are handled by [fit_synthesizer.ts()].
#' @param copula Dependence model: `"gaussian"` (default), `"t"` (heavier
#'   joint tails; degrees of freedom estimated by profile likelihood) or
#'   `"independence"`.
#' @param by Optional character vector of column names used to stratify the
#'   model. Each stratum gets its own marginals and copula; strata are sampled
#'   in proportion to their training size.
#' @param min_stratum_size Strata with fewer rows than this use a pooled model
#'   fitted on all rows.
#' @param missing How missing values are treated: `"joint"` (default; model
#'   missingness inside the copula), `"independent"` (reproduce the missing
#'   rate of each column independently) or `"drop"` (never generate missing
#'   values).
#' @param interpolation Interpolation of the empirical quantile function for
#'   continuous variables: `"linear"` (default) or `"spline"` (monotone Hyman
#'   spline, smoother densities).
#' @param discrete_threshold Numeric columns with at most this many distinct
#'   values are treated as discrete (sampled from their observed support),
#'   provided their values are whole numbers or actually repeat (fewer
#'   distinct values than half the observations). A small sample of a
#'   continuous variable is therefore still modelled as continuous.
#' @param knots Optional maximum number of quantile knots stored per
#'   continuous column. `NULL` (default) keeps every observation; a smaller
#'   value shrinks the model and avoids storing the original values.
#' @param id_cols Optional character vector of identifier columns. They are
#'   not modelled; synthetic data get fresh sequential identifiers.
#' @param pd_method How to repair a correlation matrix that is not positive
#'   definite: `"auto"`/`"higham"` (nearest correlation matrix), `"eigen"`
#'   (eigenvalue flooring) or `"none"` (fail if not positive definite).
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `sp_synthesizer`.
#'
#' @seealso [generate()], [synthesize()], [compare_synthetic()]
#'
#' @references
#' Sklar, A. (1959). Fonctions de répartition à n dimensions et leurs marges.
#' *Publications de l'Institut de Statistique de l'Université de Paris*, 8,
#' 229--231.
#'
#' Rüschendorf, L. (2009). On the distributional transform, Sklar's theorem,
#' and the empirical copula process. *Journal of Statistical Planning and
#' Inference*, 139(11), 3921--3927. \doi{10.1016/j.jspi.2009.05.030}
#'
#' Higham, N. J. (2002). Computing the nearest correlation matrix -- a
#' problem from finance. *IMA Journal of Numerical Analysis*, 22(3),
#' 329--343. \doi{10.1093/imanum/22.3.329}
#'
#' @examples
#' fit <- fit_synthesizer(iris)
#' fit
#' syn <- generate(fit, n = 300, seed = 1)
#' head(syn)
#'
#' # Stratified model and missing-data model
#' fit_aq <- fit_synthesizer(airquality, by = "Month")
#' colSums(is.na(generate(fit_aq, seed = 2)))
#' @export
fit_synthesizer <- function(data, ...) {
  UseMethod("fit_synthesizer")
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.data.frame <- function(data,
                                       copula = c("gaussian", "t", "independence"),
                                       by = NULL,
                                       min_stratum_size = 20L,
                                       missing = c("joint", "independent", "drop"),
                                       interpolation = c("linear", "spline"),
                                       discrete_threshold = 20L,
                                       knots = NULL,
                                       id_cols = NULL,
                                       pd_method = c("auto", "higham", "eigen", "none"),
                                       ...) {
  copula <- match.arg(copula)
  missing <- match.arg(missing)
  interpolation <- match.arg(interpolation)
  pd_method <- match.arg(pd_method)
  .warn_unused_dots(list(...), "fit_synthesizer",
                    if (any(c("closeness", "perturb", "dependence") %in% names(list(...)))) {
                      "'closeness', 'perturb' and 'dependence' are arguments of generate()."
                    })
  data <- .as_plain_df(data)
  if (nrow(data) < 2L) stop("'data' must have at least 2 rows.", call. = FALSE)
  if (ncol(data) < 1L) stop("'data' must have at least 1 column.", call. = FALSE)
  if (anyDuplicated(names(data)) || any(is.na(names(data)) | names(data) == "")) {
    stop("'data' must have unique, non-empty column names.", call. = FALSE)
  }
  discrete_threshold <- .check_count(discrete_threshold, "discrete_threshold")
  min_stratum_size <- .check_count(min_stratum_size, "min_stratum_size")
  if (!is.null(knots)) knots <- .check_count(knots, "knots", allow_zero = FALSE)
  if (!is.null(knots) && knots < 2L) stop("'knots' must be at least 2.", call. = FALSE)

  .check_cols(by, names(data), "by")
  .check_cols(id_cols, names(data), "id_cols")
  if (length(intersect(by, id_cols))) {
    stop("A column cannot be both in 'by' and 'id_cols'.", call. = FALSE)
  }
  model_cols <- setdiff(names(data), c(by, id_cols))
  if (!length(model_cols)) stop("No columns left to model.", call. = FALSE)
  bad <- model_cols[is.na(vapply(data[model_cols], .value_class, character(1)))]
  if (length(bad)) {
    stop("Unsupported column class in: ", paste(bad, collapse = ", "),
         ". Supported: numeric, integer, logical, character, factor, ordered, Date, POSIXct ",
         "and numeric classes such as difftime.",
         call. = FALSE)
  }

  settings <- list(copula = copula, missing = missing,
                   interpolation = interpolation,
                   discrete_threshold = discrete_threshold, knots = knots,
                   pd_method = pd_method)

  pooled <- .fit_stratum(data[model_cols], settings)

  strata <- NULL
  if (length(by)) {
    key <- .stratum_key(data[by])
    tab <- table(key, useNA = "no")
    keys <- names(tab)
    first <- match(keys, key)
    values <- data[first, by, drop = FALSE]
    attr(values, "row.names") <- .set_row_names(length(keys))
    models <- lapply(keys, function(k) {
      rows <- which(key == k)
      if (length(rows) < min_stratum_size) return(NULL)
      .fit_stratum(data[rows, model_cols, drop = FALSE], settings)
    })
    names(models) <- keys
    strata <- list(keys = keys, values = values,
                   probs = as.numeric(tab) / sum(tab),
                   counts = as.integer(tab), models = models)
  }

  structure(list(
    prototype = .prototype(data),
    columns = names(data),
    model_cols = model_cols,
    by = by,
    id_cols = id_cols,
    pooled = pooled,
    strata = strata,
    settings = settings,
    n_train = nrow(data),
    input = "data.frame"
  ), class = "sp_synthesizer")
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.matrix <- function(data, ...) {
  if (is.null(colnames(data))) colnames(data) <- paste0("V", seq_len(ncol(data)))
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  obj <- fit_synthesizer.data.frame(df, ...)
  obj$input <- "matrix"
  obj
}

#' @rdname fit_synthesizer
#' @export
fit_synthesizer.default <- function(data, ...) {
  if (!is.atomic(data) || !is.null(dim(data)) || is.na(.value_class(data))) {
    stop("Cannot fit a synthesizer to an object of class '",
         paste(class(data), collapse = "/"), "'.", call. = FALSE)
  }
  obj <- fit_synthesizer.data.frame(data.frame(x = data, stringsAsFactors = FALSE), ...)
  obj$input <- "vector"
  obj
}

.check_cols <- function(cols, available, arg) {
  if (is.null(cols)) return(invisible(NULL))
  if (!is.character(cols) || anyNA(cols)) {
    stop(sprintf("'%s' must be a character vector of column names.", arg), call. = FALSE)
  }
  miss <- setdiff(cols, available)
  if (length(miss)) {
    stop(sprintf("'%s' columns not found in data: %s", arg,
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

.stratum_key <- function(df) {
  parts <- lapply(df, function(v) ifelse(is.na(v), "<NA>", as.character(v)))
  do.call(paste, c(parts, sep = "\r"))
}

# Fit marginals + copula on one block of data
.fit_stratum <- function(df, settings) {
  marg <- lapply(names(df), function(nm) {
    .fit_marginal(df[[nm]], nm,
                  interpolation = settings$interpolation,
                  discrete_threshold = settings$discrete_threshold,
                  knots = settings$knots)
  })
  names(marg) <- names(df)

  copula_marg <- marg[vapply(marg, function(m) m$type != "empty", logical(1))]
  U <- lapply(names(copula_marg), function(nm) .pseudo_obs(copula_marg[[nm]], df[[nm]]))
  names(U) <- names(copula_marg)

  miss_marg <- list()
  if (settings$missing == "joint") {
    has_na <- vapply(marg, function(m) m$na_rate > 0 && m$na_rate < 1, logical(1))
    for (nm in names(marg)[has_na]) {
      mm <- .missing_marginal(marg[[nm]])
      miss_marg[[nm]] <- mm
      ind <- ifelse(is.na(df[[nm]]), "missing", "observed")
      U[[mm$name]] <- .pseudo_obs(mm, ind)
    }
  }

  Umat <- if (length(U)) do.call(cbind, U) else matrix(numeric(0), nrow(df), 0L)
  if (length(U)) colnames(Umat) <- names(U)
  cop <- .fit_copula(Umat, family = settings$copula,
                     pd_method = settings$pd_method)
  list(marginals = marg, missing_marginals = miss_marg, copula = cop,
       n = nrow(df))
}

# Generate n rows from one fitted block
.generate_stratum <- function(model, n, pert, prototype, missing) {
  marg <- model$marginals
  cop <- model$copula
  vars <- colnames(cop$R)
  U <- .rcopula_fitted(n, cop, dependence = pert$dependence)
  colnames(U) <- vars
  out <- vector("list", length(marg))
  names(out) <- names(marg)
  for (nm in names(marg)) {
    m <- marg[[nm]]
    u <- if (nm %in% vars) U[, nm] else stats::runif(n)
    v <- .restore_class(.perturbed_values(m, u, pert), m, prototype[[nm]])
    if (m$type != "empty" && m$na_rate > 0) {
      if (missing == "joint" && !is.null(model$missing_marginals[[nm]])) {
        mm <- model$missing_marginals[[nm]]
        is_na <- .marginal_quantile(mm, U[, mm$name]) == "missing"
      } else if (missing == "independent") {
        is_na <- stats::runif(n) < m$na_rate
      } else {
        is_na <- rep(FALSE, n)
      }
      v[is_na] <- NA
    }
    out[[nm]] <- v
  }
  attr(out, "row.names") <- .set_row_names(n)
  class(out) <- "data.frame"
  out
}

#' Generate synthetic data from a fitted model
#'
#' @description
#' `generate()` draws synthetic records from a fitted model such as an
#' `sp_synthesizer` (see [fit_synthesizer()]) or an `sp_ts_synthesizer` (see
#' [fit_synthesizer.ts()]).
#'
#' @param object A fitted model.
#' @param n Number of records (rows) to generate. Defaults to the number of
#'   rows in the training data.
#' @param seed Optional random seed. The global random number generator state
#'   is restored afterwards.
#' @param dependence Number in \[0, 1\] scaling the strength of the dependence
#'   between variables: 1 (default) reproduces the fitted copula, 0 generates
#'   independent columns. Intermediate values shrink the copula correlation
#'   matrix towards the identity; the marginal distributions are unaffected.
#' @param closeness Number in \[0, 1\] controlling how close the synthetic
#'   data are to the real data: 1 (default) reproduces the fitted model as
#'   closely as possible; smaller values move records away from the real
#'   ones (noise on continuous variables, random re-drawing of categories) and
#'   loosen the distributions (flatter category probabilities, weaker
#'   dependence). With `closeness = 0` the variables are independent and
#'   strongly smoothed. See [perturbation()] for the exact mapping and
#'   [calibrate_closeness()] to choose a value that reaches a target distance.
#' @param perturb Optional [perturbation()] for fine control: record-level
#'   noise and swapping, and distribution-level shifts, spread, category
#'   temperature and dependence strength. Combined with `closeness` and
#'   `dependence`.
#' @param stratum_sizes For stratified models: `"multinomial"` (default;
#'   random stratum sizes) or `"proportional"` (deterministic sizes
#'   proportional to the training data).
#' @param ... Further arguments passed to methods.
#'
#' @return Synthetic data of the same type as the training data: a data
#'   frame, matrix, vector or `ts` object.
#'
#' @seealso [fit_synthesizer()], [synthesize()], [generate_to_file()]
#' @examples
#' fit <- fit_synthesizer(mtcars)
#' generate(fit, n = 5, seed = 42)
#'
#' # weaker dependence
#' syn <- generate(fit, n = 500, dependence = 0.5, seed = 1)
#' @export
generate <- function(object, ...) {
  UseMethod("generate")
}

#' @rdname generate
#' @export
generate.sp_synthesizer <- function(object, n = object$n_train, seed = NULL,
                                    dependence = 1, closeness = 1, perturb = NULL,
                                    stratum_sizes = c("multinomial", "proportional"),
                                    ...) {
  .warn_unused_dots(list(...), "generate")
  n <- .check_count(n)
  dependence <- .check_prob(dependence, "dependence")
  closeness <- .check_prob(closeness, "closeness")
  stratum_sizes <- match.arg(stratum_sizes)
  pert <- .combine_perturbation(.closeness_perturbation(closeness), perturb)
  pert$dependence <- pert$dependence * dependence
  .check_pert_names(pert, object$model_cols, object$by, object$id_cols)
  .with_seed(seed, .generate_impl(object, n, pert, stratum_sizes))
}

.generate_impl <- function(object, n, pert, stratum_sizes) {
  proto <- object$prototype
  missing <- object$settings$missing

  if (is.null(object$strata)) {
    body <- .generate_stratum(object$pooled, n, pert, proto, missing)
  } else {
    st <- object$strata
    sizes <- if (stratum_sizes == "multinomial") {
      as.integer(stats::rmultinom(1L, n, st$probs))
    } else {
      .largest_remainder(n, st$probs)
    }
    parts <- vector("list", length(st$keys))
    for (i in seq_along(st$keys)) {
      if (sizes[i] == 0L) next
      mod <- st$models[[i]] %||% object$pooled
      blk <- .generate_stratum(mod, sizes[i], pert, proto, missing)
      for (b in object$by) blk[[b]] <- rep(st$values[[b]][i], sizes[i])
      parts[[i]] <- blk
    }
    parts <- parts[!vapply(parts, is.null, logical(1))]
    body <- if (length(parts)) do.call(rbind, parts) else NULL
    if (!is.null(body) && nrow(body) > 1L) {
      body <- body[sample.int(nrow(body)), , drop = FALSE]
    }
  }

  out <- proto[rep(NA_integer_, n), , drop = FALSE]
  if (!is.null(body)) {
    for (nm in setdiff(object$columns, object$id_cols)) out[[nm]] <- body[[nm]]
  }
  for (nm in object$id_cols) out[[nm]] <- .make_ids(n, proto[[nm]])
  attr(out, "row.names") <- .set_row_names(n)

  switch(object$input,
    vector = out[[1L]],
    matrix = as.matrix(out),
    out
  )
}

.largest_remainder <- function(n, p) {
  raw <- n * p
  base <- floor(raw)
  rem <- n - sum(base)
  if (rem > 0) {
    idx <- order(raw - base, decreasing = TRUE)[seq_len(rem)]
    base[idx] <- base[idx] + 1
  }
  as.integer(base)
}

.make_ids <- function(n, template) {
  ids <- seq_len(n)
  if (is.factor(template)) return(factor(sprintf("syn%0*d", nchar(n), ids)))
  if (is.character(template)) return(sprintf("syn%0*d", nchar(n), ids))
  if (is.integer(template)) return(ids)
  if (is.numeric(template)) return(as.numeric(ids))
  ids
}

#' Fit a synthesizer and generate synthetic data in one step
#'
#' @description
#' `synthesize()` is a convenience wrapper around [fit_synthesizer()] and
#' [generate()]. Optionally the result is written to a file in any format
#' supported by [write_data()].
#'
#' @param data Training data (data frame, matrix, vector or `ts`).
#' @param n Number of synthetic records. Defaults to the size of `data`.
#' @param seed Optional random seed (applies to fitting and generation).
#' @param dependence Dependence strength in \[0, 1\]; see [generate()].
#' @param closeness,perturb Closeness to the real data and optional
#'   perturbation; see [generate()].
#' @param file Optional output file path. The format is inferred from the
#'   extension (see [supported_formats()]).
#' @param format Optional explicit file format, overriding the extension.
#' @param ... Arguments passed to [fit_synthesizer()].
#'
#' @return The synthetic data (invisibly if `file` is given).
#'
#' @examples
#' syn <- synthesize(iris, n = 200, seed = 1)
#' summary(syn)
#'
#' out <- tempfile(fileext = ".csv")
#' synthesize(mtcars, file = out, seed = 1)
#' head(read_data(out))
#' @export
synthesize <- function(data, n = NULL, seed = NULL, dependence = 1,
                       closeness = 1, perturb = NULL,
                       file = NULL, format = NULL, ...) {
  res <- .with_seed(seed, {
    fit <- fit_synthesizer(data, ...)
    if (inherits(fit, "sp_ts_synthesizer")) {
      if (dependence != 1 || closeness != 1 || !is.null(perturb)) {
        warning(.ts_closeness_msg, call. = FALSE)
      }
      generate(fit, n = n %||% fit$n_train)
    } else {
      generate(fit, n = n %||% fit$n_train, dependence = dependence,
               closeness = closeness, perturb = perturb)
    }
  })
  if (!is.null(file)) {
    write_data(res, file, format = format)
    return(invisible(res))
  }
  res
}

#' Simulate several synthetic data sets
#'
#' Draws `nsim` independent synthetic data sets from a fitted model, for
#' example to account for synthesis uncertainty with multiple synthetic data
#' sets (Raghunathan, Reiter and Rubin, 2003).
#'
#' @param object A fitted `sp_synthesizer` or `sp_ts_synthesizer`.
#' @param nsim Number of synthetic data sets.
#' @param seed Optional random seed.
#' @param ... Passed to [generate()], e.g. `n`.
#' @return A list of `nsim` synthetic data sets.
#' @references Raghunathan, T. E., Reiter, J. P. and Rubin, D. B. (2003).
#'   Multiple imputation for statistical disclosure limitation. *Journal of
#'   Official Statistics*, 19(1), 1--16.
#' @examples
#' fit <- fit_synthesizer(mtcars)
#' sims <- simulate(fit, nsim = 3, seed = 1)
#' vapply(sims, function(d) mean(d$mpg), numeric(1))
#' @importFrom stats simulate
#' @export
simulate.sp_synthesizer <- function(object, nsim = 1, seed = NULL, ...) {
  nsim <- .check_count(nsim, "nsim", allow_zero = FALSE)
  .with_seed(seed, lapply(seq_len(nsim), function(i) generate(object, ...)))
}

#' @rdname simulate.sp_synthesizer
#' @export
simulate.sp_ts_synthesizer <- simulate.sp_synthesizer

#' Inspect a fitted synthesizer
#'
#' `copula_correlation()` returns the (latent) copula correlation matrix of a
#' fitted model; `summary()` returns a table describing how each column is
#' modelled.
#'
#' @param object A fitted `sp_synthesizer`.
#' @param stratum Optional stratum label for stratified models (as shown by
#'   `summary()`); defaults to the pooled model.
#' @param include_missing Include the missingness indicators in the matrix?
#' @param ... Unused.
#' @return `copula_correlation()`: a correlation matrix. `summary()`: an
#'   object of class `summary.sp_synthesizer` (a list with a `columns` data
#'   frame).
#' @examples
#' fit <- fit_synthesizer(airquality)
#' round(copula_correlation(fit), 2)
#' summary(fit)
#' @export
copula_correlation <- function(object, stratum = NULL, include_missing = FALSE) {
  if (!inherits(object, "sp_synthesizer")) {
    stop("'object' must be an sp_synthesizer.", call. = FALSE)
  }
  mod <- object$pooled
  if (!is.null(stratum)) {
    if (is.null(object$strata)) stop("The model is not stratified.", call. = FALSE)
    i <- match(stratum, .stratum_label(object$strata$keys))
    if (is.na(i)) stop("Unknown stratum: ", stratum, call. = FALSE)
    mod <- object$strata$models[[i]] %||% object$pooled
  }
  R <- mod$copula$R
  if (!include_missing) {
    keep <- !startsWith(colnames(R), ".missing_")
    R <- R[keep, keep, drop = FALSE]
  }
  R
}

.stratum_label <- function(keys) gsub("\r", " / ", keys, fixed = TRUE)

#' @rdname copula_correlation
#' @export
summary.sp_synthesizer <- function(object, ...) {
  m <- object$pooled$marginals
  cols <- data.frame(
    column = names(m),
    class = vapply(m, `[[`, character(1), "class"),
    model = vapply(m, `[[`, character(1), "type"),
    support = vapply(m, function(z) {
      switch(z$type,
        continuous = sprintf("[%s, %s]", format(signif(min(z$quantiles), 4)),
                             format(signif(max(z$quantiles), 4))),
        discrete = , categorical = sprintf("%d values", length(z$values)),
        empty = "-")
    }, character(1)),
    missing = vapply(m, function(z) z$na_rate, numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  strata <- NULL
  if (!is.null(object$strata)) {
    st <- object$strata
    strata <- data.frame(stratum = .stratum_label(st$keys), n = st$counts,
                         model = ifelse(vapply(st$models, is.null, logical(1)),
                                        "pooled", "own"),
                         stringsAsFactors = FALSE)
  }
  structure(list(columns = cols, strata = strata, settings = object$settings,
                 n_train = object$n_train,
                 df = object$pooled$copula$df),
            class = "summary.sp_synthesizer")
}

#' @export
print.summary.sp_synthesizer <- function(x, ...) {
  cat("Synthesizer trained on", x$n_train, "rows\n")
  cat("Copula:", x$settings$copula)
  if (is.finite(x$df)) cat(sprintf(" (df = %g)", x$df))
  cat(" | missing data:", x$settings$missing, "\n\n")
  cols <- x$columns
  cols$missing <- .fmt_pct(cols$missing)
  print(cols, row.names = FALSE)
  if (!is.null(x$strata)) {
    cat("\nStrata:\n")
    print(x$strata, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.sp_synthesizer <- function(x, ...) {
  m <- x$pooled$marginals
  types <- table(vapply(m, `[[`, character(1), "type"))
  cat("<sp_synthesizer>\n")
  cat(sprintf("  trained on : %d rows, %d modelled column(s)\n",
              x$n_train, length(x$model_cols)))
  cat("  marginals  :", paste(sprintf("%d %s", as.integer(types), names(types)),
                              collapse = ", "), "\n")
  cop <- x$pooled$copula
  cat("  copula     :", cop$family,
      if (is.finite(cop$df)) sprintf("(df = %g)", cop$df) else "", "\n")
  cat("  missing    :", x$settings$missing, "\n")
  if (!is.null(x$strata)) {
    own <- sum(!vapply(x$strata$models, is.null, logical(1)))
    cat(sprintf("  strata     : %d by %s (%d with own model)\n",
                length(x$strata$keys), paste(x$by, collapse = " x "), own))
  }
  if (length(x$id_cols)) cat("  id columns :", paste(x$id_cols, collapse = ", "), "\n")
  invisible(x)
}
