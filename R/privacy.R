# Disclosure risk for synthetic microdata -------------------------------------

#' Disclosure risk of synthetic microdata
#'
#' @description
#' Measures the two classical disclosure risks of a synthetic version of
#' confidential microdata (census, survey, register data), from the point of
#' view of an intruder who knows the *key variables* (quasi-identifiers such
#' as age, sex and region) of a person in the real data.
#'
#' * **Identity disclosure -- replicated uniques.** Records whose combination
#'   of key variables is unique in the real data are the most exposed. The
#'   share of these real uniques whose key combination also occurs exactly
#'   once in the synthetic data (`replicated_uniques`, Nowok, Raab and
#'   Dibben, 2016) indicates how many synthetic records could be taken for a
#'   real, identifiable person.
#' * **Attribute disclosure -- correct attribution probability (CAP).** For
#'   every real record, the intruder looks up the synthetic records with the
#'   same key values and predicts each `target` (sensitive) variable: the
#'   most frequent category, or the median for numeric targets. `cap` is the
#'   share of real records for which the prediction is correct (numeric
#'   targets: within `tolerance` relative error), among those with at least
#'   one matching synthetic record. `cap_baseline` is the same share when the
#'   intruder ignores the synthetic data and predicts the overall mode or
#'   median; `cap_uniques` restricts the CAP to real uniques. A CAP close to
#'   the baseline means that the synthetic data reveal little beyond what
#'   is known from population-level statistics (Taub et al., 2018).
#'
#' Numeric key variables are matched exactly; coarsen them first (e.g. age
#' in years or 5-year bands) if they are continuous.
#'
#' @param real,synthetic Data frames with the key and target variables.
#' @param keys Character vector of key variables (quasi-identifiers).
#' @param target Character vector of sensitive variables for attribute
#'   disclosure (optional).
#' @param ignore Optional columns to leave out entirely, typically direct
#'   identifiers (which differ between real and synthetic files and would
#'   otherwise hide exact copies).
#' @param tolerance Relative tolerance for a correct prediction of numeric
#'   targets (default 0.1, i.e. 10 percent).
#' @return An object of class `sp_disclosure` with elements
#'   `replicated_uniques`, `real_unique_rate` (share of real records that are
#'   unique on the keys), `exact_copies` (share of synthetic records that are
#'   identical to some real record on all common variables except `ignore`)
#'   and a data frame
#'   `attribute` with, per target, `match_rate`, `cap`, `cap_baseline` and
#'   `cap_uniques`.
#' @references
#' Nowok, B., Raab, G. M. and Dibben, C. (2016). synthpop: Bespoke creation of
#' synthetic data in R. *Journal of Statistical Software*, 74(11), 1--26.
#' \doi{10.18637/jss.v074.i11}
#'
#' Taub, J., Elliot, M., Pampaka, M. and Smith, D. (2018). Differential
#' correct attribution probability for synthetic data: an exploration. In
#' *Privacy in Statistical Databases*, LNCS 11126, 122--137. Springer.
#' \doi{10.1007/978-3-319-99771-1_9}
#' @examples
#' syn <- synthesize(mtcars, seed = 1)
#' disclosure_risk(mtcars, syn, keys = c("cyl", "gear", "am"), target = "mpg")
#'
#' # census-like microdata with a direct identifier
#' fit <- fit_synthesizer(census_sim, id_cols = "person_id", by = "employment")
#' syn <- generate(fit, seed = 1)
#' disclosure_risk(census_sim, syn, keys = c("age", "sex", "region"),
#'                 target = c("income", "education"), ignore = "person_id")
#' @export
disclosure_risk <- function(real, synthetic, keys, target = NULL, ignore = NULL,
                            tolerance = 0.1) {
  p <- .check_pair(real, synthetic)
  real <- p$real
  synthetic <- p$synthetic
  if (missing(keys) || !length(keys)) stop("'keys' must name at least one variable.", call. = FALSE)
  .check_cols(keys, names(real), "keys")
  .check_cols(target, names(real), "target")
  if (length(intersect(keys, target))) stop("A variable cannot be both a key and a target.", call. = FALSE)
  .check_cols(ignore, names(real), "ignore")
  if (length(intersect(ignore, c(keys, target)))) {
    stop("'ignore' cannot contain key or target variables.", call. = FALSE)
  }
  keep <- setdiff(names(real), ignore)
  real <- real[keep]
  synthetic <- synthetic[keep]
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0) {
    stop("'tolerance' must be a non-negative number.", call. = FALSE)
  }

  kk <- .row_keys(real[keys], synthetic[keys])
  kr <- kk[[1L]]
  ks <- kk[[2L]]
  freq_r <- table(kr)
  freq_s <- table(ks)
  uniq_r <- names(freq_r)[freq_r == 1L]
  rep_uniq <- if (length(uniq_r)) {
    mean(uniq_r %in% names(freq_s)[freq_s == 1L])
  } else {
    0
  }

  full <- .row_keys(real, synthetic)
  full_r <- full[[1L]]
  full_s <- full[[2L]]
  exact <- mean(full_s %in% full_r)

  attr_tab <- NULL
  if (length(target)) {
    is_uniq <- kr %in% uniq_r
    rows <- lapply(target, function(tg) {
      yr <- real[[tg]]
      ys <- synthetic[[tg]]
      num <- .is_numeric_like(yr)
      pred <- .predict_by_key(ks, ys, num)
      hat <- pred[kr]
      matched <- !is.na(names(pred)[match(kr, names(pred))]) & !is.na(yr)
      correct <- .correct(hat, yr, num, tolerance)
      base_val <- if (num) {
        stats::median(as.numeric(yr), na.rm = TRUE)
      } else {
        names(which.max(table(as.character(yr))))
      }
      base_correct <- .correct(rep(base_val, length(yr)), yr, num, tolerance)
      data.frame(
        target = tg,
        match_rate = mean(matched),
        cap = if (any(matched)) mean(correct[matched]) else NA_real_,
        cap_baseline = mean(base_correct[!is.na(yr)]),
        cap_uniques = if (any(matched & is_uniq)) mean(correct[matched & is_uniq]) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    attr_tab <- do.call(rbind, rows)
  }
  structure(list(keys = keys, replicated_uniques = rep_uniq,
                 real_unique_rate = mean(kr %in% uniq_r),
                 n_real_uniques = length(uniq_r), exact_copies = exact,
                 attribute = attr_tab),
            class = "sp_disclosure")
}

# Prediction of a target per key combination: mode or median
.predict_by_key <- function(key, y, numeric) {
  ok <- !is.na(y)
  if (!any(ok)) return(stats::setNames(character(0), character(0)))
  if (numeric) {
    out <- tapply(as.numeric(y[ok]), key[ok], stats::median)
  } else {
    out <- tapply(as.character(y[ok]), key[ok], function(v) {
      tb <- table(v)
      names(tb)[which.max(tb)]
    })
  }
  stats::setNames(as.vector(out), names(out))
}

.correct <- function(hat, y, numeric, tolerance) {
  res <- rep(FALSE, length(y))
  ok <- !is.na(hat) & !is.na(y)
  if (numeric) {
    yy <- as.numeric(y[ok])
    hh <- as.numeric(hat[ok])
    res[ok] <- abs(hh - yy) <= tolerance * pmax(abs(yy), .Machine$double.eps)
  } else {
    res[ok] <- as.character(hat[ok]) == as.character(y[ok])
  }
  res
}

#' @export
print.sp_disclosure <- function(x, digits = 3, ...) {
  cat("<sp_disclosure> keys:", paste(x$keys, collapse = ", "), "\n")
  cat(sprintf("  real records unique on the keys : %s (%d combinations)\n",
              .fmt_pct(x$real_unique_rate), x$n_real_uniques))
  cat(sprintf("  replicated uniques              : %s of real uniques\n",
              .fmt_pct(x$replicated_uniques)))
  cat(sprintf("  synthetic records = real record : %s\n", .fmt_pct(x$exact_copies, 2)))
  if (!is.null(x$attribute)) {
    cat("\nAttribute disclosure (correct attribution probability):\n")
    a <- x$attribute
    num <- vapply(a, is.numeric, logical(1))
    a[num] <- lapply(a[num], round, digits)
    print(a, row.names = FALSE)
  }
  invisible(x)
}
