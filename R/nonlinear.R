# Optional hybrid copula / conditional synthesis (CART or GAM) -----------------
#
# The copula remains responsible for the predictors and missingness. For ONE
# specified, continuously valued response we replace the generated response by
# a fitted conditional prediction (CART or penalised GAM) plus a residual draw.
# The result is optionally marginally quantile-calibrated to the observed response. This
# creates a semi-parametric conditional synthesizer; its final joint
# distribution is NOT an elliptical copula, even when the baseline is Gaussian.
#
# The tree and its residual pools contain information from the training data.
# Fitted models must NOT be disclosed as if they were privacy-safe artifacts.

.hybrid_encode <- function(df, specs = NULL) {
  train <- is.null(specs)
  if (train) specs <- vector("list", ncol(df))
  if (length(specs) != ncol(df)) stop("Hybrid predictor specification mismatch.", call. = FALSE)
  out <- list()
  for (j in seq_along(df)) {
    x <- df[[j]]
    stem <- paste0("x", j)
    if (train) {
      is_cat <- is.factor(x) || is.character(x) || is.logical(x)
      if (is_cat) {
        lv <- sort(unique(as.character(x[!is.na(x)])), method = "radix")
        if (length(lv) > 30L) {
          stop("Hybrid CART currently supports at most 30 observed levels per categorical predictor.",
               call. = FALSE)
        }
        specs[[j]] <- list(kind = "categorical", values = lv)
      } else {
        xx <- suppressWarnings(as.numeric(x))
        if (length(xx) != length(x)) stop("Unsupported nonlinear predictor.", call. = FALSE)
        med <- suppressWarnings(stats::median(xx, na.rm = TRUE))
        if (!is.finite(med)) med <- 0
        specs[[j]] <- list(kind = "numeric", median = med)
      }
    }
    sp <- specs[[j]]
    if (identical(sp$kind, "categorical")) {
      chars <- as.character(x)
      code <- match(chars, sp$values)
      if (any(!is.na(chars) & is.na(code))) {
        stop("Hybrid generation encountered a previously unseen categorical predictor.",
             call. = FALSE)
      }
      labels <- paste0("L", seq_along(sp$values))
      tokens <- ifelse(is.na(chars), "MISSING", paste0("L", code))
      out[[stem]] <- factor(tokens, levels = c(labels, "MISSING"))
    } else if (identical(sp$kind, "numeric")) {
      z <- suppressWarnings(as.numeric(x))
      if (any(!is.na(z) & !is.finite(z))) {
        stop("Hybrid predictor contains non-finite numeric values.", call. = FALSE)
      }
      out[[stem]] <- ifelse(is.na(z), sp$median, z)
      out[[paste0(stem, "_missing")]] <- as.integer(is.na(z))
    } else {
      stop("Unsupported hybrid predictor specification.", call. = FALSE)
    }
  }
  list(data = as.data.frame(out, stringsAsFactors = FALSE), specs = specs)
}

.fit_nonlinear <- function(df, settings) {
  engine <- settings$nonlinear_engine %||% "cart"
  if (identical(engine, "auto")) {
    selection <- .select_nonlinear_engine(df, settings)
    chosen <- settings
    chosen$nonlinear_engine <- selection$engine
    chosen$nonlinear_gam_interactions <- selection$interactions
    model <- .fit_nonlinear(df, chosen)
    if (!is.null(model)) model$selection <- selection
    return(model)
  }
  if (identical(engine, "gam")) {
    return(.fit_nonlinear_gam(df, settings))
  }
  target <- settings$nonlinear_target
  if (is.null(target)) return(NULL)
  y <- df[[target]]
  ok <- is.finite(y) & !is.na(y)
  if (sum(ok) < max(10L, 2L * settings$nonlinear_minbucket)) {
    # A sparsely observed stratum falls back to the fitted copula.
    return(NULL)
  }
  predictors <- setdiff(names(df), target)
  if (!length(predictors)) return(NULL)
  # Encode on ALL rows, so that categories occurring only where the target is
  # missing are known to the encoder (the copula can generate them).
  encoded <- .hybrid_encode(df[predictors])
  encoded$data <- encoded$data[ok, , drop = FALSE]
  seen <- .seen_tokens(encoded$data)
  train <- data.frame(.response = as.numeric(y[ok]), encoded$data, check.names = FALSE)
  ctrl <- rpart::rpart.control(
    cp = settings$nonlinear_cp,
    minbucket = settings$nonlinear_minbucket,
    minsplit = max(2L * settings$nonlinear_minbucket, 4L),
    maxdepth = settings$nonlinear_maxdepth,
    xval = 0L
  )
  tree <- rpart::rpart(.response ~ ., data = train, method = "anova",
                       control = ctrl, model = FALSE, x = FALSE, y = FALSE)
  fitted <- as.numeric(stats::predict(tree, newdata = encoded$data))
  if (length(fitted) != sum(ok) || any(!is.finite(fitted))) {
    stop("Nonlinear tree fit returned invalid predictions.", call. = FALSE)
  }
  centers <- sort(unique(fitted))
  ix <- match(fitted, centers)
  residuals <- as.numeric(y[ok]) - fitted
  pools <- split(residuals, factor(ix, levels = seq_along(centers)))
  # Center each residual pool. A leaf's fitted mean should already be its
  # sample mean, but explicit centering reduces rounding effects.
  pools <- lapply(pools, function(z) z - mean(z))
  model <- list(target = target, predictors = predictors, encoder = encoded$specs,
                tree = tree, centers = centers, pools = pools, seen = seen,
                residual_mode = settings$nonlinear_residual %||% "leaf",
                calibration_mode = settings$nonlinear_calibration %||% "quantile")
  if (!isTRUE(settings$.cv_fit)) model$reference <- .calibration_reference(model, df)
  model
}

# Categories (tokens) of each encoded factor column that occur in the rows
# used to fit the conditional model
.seen_tokens <- function(encoded) {
  out <- list()
  for (nm in names(encoded)) {
    if (is.factor(encoded[[nm]])) out[[nm]] <- unique(as.character(encoded[[nm]]))
  }
  out
}

# Categories without any observed response cannot inform the tree: they are
# passed to rpart as unknown (NA), which it handles through surrogate splits,
# instead of being sent down an arbitrary branch.
.mask_unseen <- function(encoded, seen) {
  if (is.null(seen)) return(encoded)
  for (nm in names(seen)) {
    x <- encoded[[nm]]
    if (is.null(x)) next
    x[!(as.character(x) %in% seen[[nm]])] <- NA
    encoded[[nm]] <- x
  }
  encoded
}

# Reference distribution for the marginal calibration of the hybrid target:
# uncalibrated conditional draws for a large resample of the training
# predictor rows. Calibrating each generated value against this fixed
# reference (rather than against the other values of the same batch) makes
# the result independent of the number of rows generated per call, stratum
# or chunk. The reference is drawn from its own fixed random stream, so it is
# a deterministic function of the fitted model and does not consume or alter
# the caller's random number stream.
.calibration_reference <- function(model, df, size = 20000L) {
  if (identical(model$calibration_mode, "none")) return(NULL)
  .with_seed(20260101L, {
    rows <- sample.int(nrow(df), size, replace = TRUE)
    tmp <- model
    tmp$calibration_mode <- "none"
    sort(.generate_nonlinear(tmp, df[rows, , drop = FALSE], NULL))
  })
}

# Map raw draws to (0, 1) by the mid-point empirical CDF of the reference
.calibrate_u <- function(raw, reference, residual_mode) {
  if (is.null(reference) || !length(reference)) {
    # models fitted before the reference was stored: within-batch ranks
    return(rank(raw, ties.method = if (residual_mode == "none") "average" else "random") /
             (length(raw) + 1))
  }
  N <- length(reference)
  u <- (findInterval(raw, reference, left.open = TRUE) + findInterval(raw, reference)) / (2 * N)
  pmin(pmax(u, 0.5 / N), 1 - 0.5 / N)
}

.generate_nonlinear <- function(hybrid, df, marginal) {
  if (identical(hybrid$engine %||% "cart", "gam")) {
    return(.generate_nonlinear_gam(hybrid, df, marginal))
  }
  enc <- .hybrid_encode(df[hybrid$predictors], hybrid$encoder)
  pred <- as.numeric(stats::predict(hybrid$tree, newdata = .mask_unseen(enc$data, hybrid$seen)))
  n <- nrow(df)
  if (length(pred) != n || any(!is.finite(pred))) {
    stop("Nonlinear tree prediction returned invalid values.", call. = FALSE)
  }
  # Prediction values are terminal-node means; mapping to the set of
  # observed terminal predictions identifies the residual pool. A nearest
  # center fallback handles tiny floating-point differences.
  ix <- match(pred, hybrid$centers)
  if (anyNA(ix)) {
    for (i in which(is.na(ix))) ix[i] <- which.min(abs(hybrid$centers - pred[i]))
  }
  residual_mode <- hybrid$residual_mode %||% "leaf"
  calibration_mode <- hybrid$calibration_mode %||% "quantile"
  # Three ablations on exactly the same fitted conditional mean and predictor
  # distribution. The original implementation is residual_mode="leaf".
  raw <- pred
  if (identical(residual_mode, "leaf")) {
    for (k in unique(ix)) {
      rows <- which(ix == k)
      pool <- hybrid$pools[[k]]
      raw[rows] <- pred[rows] + pool[sample.int(length(pool), length(rows), replace = TRUE)]
    }
  } else if (identical(residual_mode, "global")) {
    pool <- unlist(hybrid$pools, use.names = FALSE)
    if (!length(pool)) stop("Hybrid residual pool is empty.", call. = FALSE)
    raw <- pred + pool[sample.int(length(pool), n, replace = TRUE)]
  } else if (!identical(residual_mode, "none")) {
    stop("Unknown hybrid residual mode in fitted model.", call. = FALSE)
  }
  # 'none' returns raw conditional draws; these can exceed the training
  # response's observed range. 'quantile' is the pre-existing default.
  if (identical(calibration_mode, "none")) return(raw)
  if (!identical(calibration_mode, "quantile")) {
    stop("Unknown hybrid calibration mode in fitted model.", call. = FALSE)
  }
  # Calibrate against the stored reference distribution. With residual mode
  # "none" tied terminal means receive the mid-point of their tie block, so
  # no artificial within-leaf noise is introduced.
  .marginal_quantile(marginal, .calibrate_u(raw, hybrid$reference, residual_mode))
}

# GAM hybrid implementation ---------------------------------------------------
# GAM is an alternative conditional mean model, not an alternative copula.
# Factor predictors enter parametrically. Numeric columns with enough unique
# values use penalised cubic regression splines. Columns with two or three
# distinct values enter linearly. Predictors constant in the observed-response
# sample are omitted. Interaction terms are *opt in*.
#
# For GAM, the pre-existing "leaf" residual setting is interpreted as
# nearest-neighbour resampling in the one-dimensional fitted-mean space;
# unlike a CART leaf this cannot capture heteroskedasticity orthogonal to E[Y|X].

.gam_predictor_frame <- function(encoded, observed_levels = NULL) {
  frame <- encoded
  if (is.null(observed_levels)) observed_levels <- list()
  reference <- list()
  for (name in names(frame)) {
    x <- frame[[name]]
    if (!is.factor(x)) next
    if (is.null(observed_levels[[name]])) {
      levels_here <- levels(droplevels(x))
      if (!length(levels_here)) next
      observed_levels[[name]] <- levels_here
    }
    lv <- observed_levels[[name]]
    values <- as.character(x)
    values[is.na(values) | !(values %in% lv)] <- lv[1L]
    frame[[name]] <- factor(values, levels = lv)
  }
  list(data = frame, observed_levels = observed_levels)
}

.fit_nonlinear_gam <- function(df, settings) {
  target <- settings$nonlinear_target
  y <- df[[target]]
  ok <- !is.na(y) & is.finite(y)
  if (sum(ok) < max(10L, 2L * settings$nonlinear_minbucket)) return(NULL)
  predictors <- setdiff(names(df), target)
  if (!length(predictors)) return(NULL)
  # Build encoder on all predictor rows to preserve categories occurring only
  # among training rows whose target is missing.
  enc <- .hybrid_encode(df[predictors])
  trained <- .gam_predictor_frame(enc$data[ok, , drop = FALSE])
  train <- trained$data
  train$.response <- as.numeric(y[ok])
  terms <- character(0L)
  smooth_terms <- character(0L)
  k_values <- integer(0L)
  for (name in setdiff(names(train), ".response")) {
    x <- train[[name]]
    if (is.factor(x)) {
      if (nlevels(x) > 1L) terms <- c(terms, name)
    } else if (is.numeric(x)) {
      unique_n <- length(unique(x[is.finite(x)]))
      if (unique_n >= 4L && !grepl("_missing$", name)) {
        kk <- min(as.integer(settings$nonlinear_gam_k), unique_n,
                  max(4L, floor(sum(ok) / 4)))
        terms <- c(terms, sprintf("s(%s, bs='cr', k=%d)", name, kk))
        smooth_terms <- c(smooth_terms, name)
        k_values[name] <- kk
      } else if (unique_n > 1L) {
        terms <- c(terms, name)
      }
    }
  }
  if (isTRUE(settings$nonlinear_gam_interactions) && length(smooth_terms) > 1L) {
    if (length(smooth_terms) > 3L) {
      stop("GAM tensor interactions require at most three eligible numeric ",
           "smooth predictors; disable 'nonlinear_gam_interactions' or ",
           "reduce the predictor set.", call. = FALSE)
    }
    pair_terms <- utils::combn(smooth_terms, 2L, simplify = FALSE)
    for (pair in pair_terms) {
      kk <- pmin(as.integer(k_values[pair]), 5L)
      terms <- c(terms, sprintf("ti(%s, %s, bs=c('cr','cr'), k=c(%d,%d))",
                                pair[1L], pair[2L], kk[1L], kk[2L]))
    }
  }
  # Set the formula environment to mgcv's namespace: s() and ti() are found
  # without attaching mgcv, and no user-controlled column names are pasted.
  f <- stats::as.formula(paste(".response ~", if (length(terms))
    paste(terms, collapse = " + ") else "1"), env = asNamespace("mgcv"))
  model <- mgcv::gam(f, data = train, method = "REML", select = TRUE,
                     drop.unused.levels = FALSE)
  fitted <- as.numeric(stats::predict(model, newdata = train))
  if (length(fitted) != sum(ok) || any(!is.finite(fitted))) {
    stop("Nonlinear GAM fit returned invalid predictions.", call. = FALSE)
  }
  residual <- as.numeric(y[ok]) - fitted
  residual <- residual - mean(residual)
  order_fit <- order(fitted)
  out <- list(target = target, predictors = predictors, encoder = enc$specs,
              engine = "gam", gam = model,
              observed_levels = trained$observed_levels,
              fitted_sorted = fitted[order_fit], residual_sorted = residual[order_fit],
              residual_mode = settings$nonlinear_residual %||% "leaf",
              calibration_mode = settings$nonlinear_calibration %||% "quantile",
              local_k = settings$nonlinear_local_k %||% 30L)
  if (!isTRUE(settings$.cv_fit)) out$reference <- .calibration_reference(out, df)
  out
}

.generate_nonlinear_gam <- function(hybrid, df, marginal) {
  n <- nrow(df)
  if (!n) return(numeric(0L))
  enc <- .hybrid_encode(df[hybrid$predictors], hybrid$encoder)
  prepared <- .gam_predictor_frame(enc$data, hybrid$observed_levels)
  pred <- as.numeric(stats::predict(hybrid$gam, newdata = prepared$data))
  if (length(pred) != n || any(!is.finite(pred))) {
    stop("Nonlinear GAM prediction returned invalid values.", call. = FALSE)
  }
  residual_mode <- hybrid$residual_mode %||% "leaf"
  calibration_mode <- hybrid$calibration_mode %||% "quantile"
  raw <- pred
  if (identical(residual_mode, "global")) {
    res <- hybrid$residual_sorted
    raw <- pred + res[sample.int(length(res), n, replace = TRUE)]
  } else if (identical(residual_mode, "leaf")) {
    # Binary search for each fitted mean, then a draw from the window of the
    # k nearest training fitted values; each window's residuals are centred
    # so the GAM conditional mean is kept. Window means come from cumulative
    # sums, so the whole step is vectorised.
    xp <- hybrid$fitted_sorted
    rp <- hybrid$residual_sorted
    available <- length(xp)
    k <- min(as.integer(hybrid$local_k), available)
    center <- findInterval(pred, xp)
    lo <- pmin(pmax(1L, center - k %/% 2L), available - k + 1L)
    cs <- c(0, cumsum(rp))
    window_mean <- (cs[lo + k] - cs[lo]) / k
    pick <- lo + sample.int(k, n, replace = TRUE) - 1L
    raw <- pred + rp[pick] - window_mean
  } else if (!identical(residual_mode, "none")) {
    stop("Unknown GAM residual mode in fitted model.", call. = FALSE)
  }
  if (identical(calibration_mode, "none")) return(raw)
  if (!identical(calibration_mode, "quantile")) {
    stop("Unknown GAM marginal calibration mode in fitted model.", call. = FALSE)
  }
  .marginal_quantile(marginal, .calibrate_u(raw, hybrid$reference, residual_mode))
}


# Training-only cross-validation for selecting the conditional engine.
# This deliberately does NOT use the fitted copula, generated synthetic data,
# downstream model, or the benchmark's independent holdout. By default each
# candidate is scored by the continuous ranked probability score (CRPS) of
# its conditional predictive distribution -- the conditional mean PLUS the
# residual resampling actually used for generation -- on held-out real rows.
# CRPS rewards a correct conditional spread as well as a correct conditional
# mean; the mean-squared error of the conditional mean is reported alongside
# and can be used for selection with nonlinear_cv_score = "mse".
# Selection does not directly optimise disclosure risk.
#
# Each fold's encoder/model is trained exclusively on the other folds.
# Unknown nominal levels in a validation fold are treated as missing for
# prediction, rather than leaking the held-out level vocabulary into fitting.
.nonlinear_cv_newdata <- function(df, model) {
  newdata <- df[model$predictors]
  for (j in seq_along(model$predictors)) {
    nm <- model$predictors[j]
    spec <- model$encoder[[j]]
    if (!identical(spec$kind, "categorical")) next
    v <- as.character(newdata[[nm]])
    v[!is.na(v) & !(v %in% spec$values)] <- NA_character_
    newdata[[nm]] <- v
  }
  newdata
}

.predict_nonlinear_mean <- function(model, df) {
  predictors <- .nonlinear_cv_newdata(df, model)
  enc <- .hybrid_encode(predictors, model$encoder)
  if (identical(model$engine %||% "cart", "gam")) {
    prepared <- .gam_predictor_frame(enc$data, model$observed_levels)
    pred <- stats::predict(model$gam, newdata = prepared$data)
  } else {
    pred <- stats::predict(model$tree, newdata = .mask_unseen(enc$data, model$seen))
  }
  as.numeric(pred)
}

# Sample CRPS of y under the empirical distribution of m conditional draws
# per held-out row: mean|X - y| - 0.5 * mean|X - X'|.
.crps_draws <- function(model, df, y, m = 50L) {
  nd <- df
  nd[model$predictors] <- .nonlinear_cv_newdata(df, model)
  rep_rows <- rep(seq_len(nrow(nd)), each = m)
  tmp <- model
  tmp$calibration_mode <- "none"
  draws <- .generate_nonlinear(tmp, nd[rep_rows, , drop = FALSE], NULL)
  if (length(draws) != length(rep_rows) || any(!is.finite(draws))) {
    stop("non-finite or incomplete conditional draws")
  }
  X <- matrix(draws, nrow = m)                    # one column per held-out row
  term1 <- colMeans(abs(X - rep(y, each = m)))
  Xs <- apply(X, 2L, sort)
  w <- (2 * seq_len(m) - m - 1) / m^2             # E|X - X'| from order statistics
  term2 <- colSums(Xs * w)
  mean(term1 - term2)
}

.select_nonlinear_engine <- function(df, settings) {
  target <- settings$nonlinear_target
  valid <- which(!is.na(df[[target]]) & is.finite(df[[target]]))
  folds_requested <- settings$nonlinear_cv_folds %||% 5L
  rule <- settings$nonlinear_cv_rule %||% "one_se"
  criterion <- settings$nonlinear_cv_score %||% "crps"
  # Candidate order determines the complexity preference used by 1-SE.
  candidates <- data.frame(name = c("cart", "gam", "gam_ti"),
                           engine = c("cart", "gam", "gam"),
                           interactions = c(FALSE, FALSE, TRUE),
                           complexity = 1:3, stringsAsFactors = FALSE)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    candidates <- candidates[1L, , drop = FALSE]
  }
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Automatic nonlinear engine selection requires suggested package 'rpart'.",
         call. = FALSE)
  }
  k <- min(as.integer(folds_requested), length(valid))
  min_n <- max(10L, 2L * settings$nonlinear_minbucket)
  # Sparse strata can support the full-data CART model but not K-fold fitting.
  # In that case expose an explicit, reproducible fallback to CART.
  if (k < 2L || length(valid) - ceiling(length(valid) / max(k, 1L)) < min_n) {
    return(list(engine = "cart", interactions = FALSE, selected = "cart",
                reason = "insufficient_observed_responses_for_cv",
                rule = rule, score = criterion, folds = 0L, scores = data.frame(),
                fold_scores = matrix(numeric(), nrow = 0L, ncol = 0L),
                fold_mse = matrix(numeric(), nrow = 0L, ncol = 0L)))
  }
  fold_ids <- integer(nrow(df))
  fold_ids[valid] <- sample(rep(seq_len(k), length.out = length(valid)))
  scores <- matrix(NA_real_, nrow = k, ncol = nrow(candidates),
                   dimnames = list(paste0("fold", seq_len(k)), candidates$name))
  mse <- scores
  issues <- rep("", nrow(candidates))
  for (j in seq_len(nrow(candidates))) {
    cfg <- settings
    cfg$.cv_fit <- TRUE
    cfg$nonlinear_engine <- candidates$engine[j]
    cfg$nonlinear_gam_interactions <- candidates$interactions[j]
    for (f in seq_len(k)) {
      test_idx <- which(fold_ids == f)
      fit_df <- df[fold_ids != f, , drop = FALSE]
      # Leave original rows with missing target in the fitting fold. They can
      # help encode predictor categories but do not supply response outcomes.
      res <- tryCatch({
        model <- .fit_nonlinear(fit_df, cfg)
        if (is.null(model)) stop("conditional model unavailable in fold")
        test_df <- df[test_idx, , drop = FALSE]
        y_test <- as.numeric(df[[target]][test_idx])
        pred <- .predict_nonlinear_mean(model, test_df)
        if (length(pred) != length(test_idx) || any(!is.finite(pred))) {
          stop("non-finite or incomplete conditional prediction")
        }
        c(mse = mean((y_test - pred)^2),
          crps = if (identical(criterion, "crps")) .crps_draws(model, test_df, y_test) else NA_real_)
      }, error = identity)
      if (inherits(res, "error")) {
        issues[j] <- conditionMessage(res)
        break
      }
      mse[f, j] <- res[["mse"]]
      scores[f, j] <- if (identical(criterion, "crps")) res[["crps"]] else res[["mse"]]
    }
  }
  # Reject a candidate when any fold failed; no silent complete-case averaging.
  complete <- colSums(is.finite(scores)) == k
  mu <- rep(NA_real_, ncol(scores))
  se <- rep(NA_real_, ncol(scores))
  mu[complete] <- colMeans(scores[, complete, drop = FALSE])
  se[complete] <- apply(scores[, complete, drop = FALSE], 2L,
                        stats::sd) / sqrt(k)
  mu_mse <- se_mse <- rep(NA_real_, ncol(mse))
  mu_mse[complete] <- colMeans(mse[, complete, drop = FALSE])
  se_mse[complete] <- apply(mse[, complete, drop = FALSE], 2L, stats::sd) / sqrt(k)
  score_table <- data.frame(candidate = candidates$name,
                            engine = candidates$engine,
                            interactions = candidates$interactions,
                            complexity = candidates$complexity,
                            n_valid = colSums(is.finite(scores)),
                            mean_score = mu, se_score = se,
                            mean_mse = mu_mse, se_mse = se_mse,
                            issue = issues, stringsAsFactors = FALSE)
  available <- which(complete)
  if (!length(available)) {
    # Avoid training data leakage and remain usable on challenging strata.
    # If the full CART refit also fails, the caller reports that failure.
    return(list(engine = "cart", interactions = FALSE, selected = "cart",
                reason = "all_cv_candidates_failed", rule = rule, score = criterion,
                folds = k, scores = score_table, fold_scores = scores, fold_mse = mse))
  }
  best <- available[which.min(mu[available])]
  chosen <- if (identical(rule, "one_se")) {
    eligible <- available[mu[available] <= mu[best] + se[best] + 1e-12]
    eligible[which.min(candidates$complexity[eligible])]
  } else best
  list(engine = candidates$engine[chosen],
       interactions = candidates$interactions[chosen],
       selected = candidates$name[chosen], reason = "cv",
       rule = rule, score = criterion, folds = k, scores = score_table,
       fold_scores = scores, fold_mse = mse)
}
