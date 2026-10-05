# SynthesizerPlus: augmenting a small training set for machine learning
# Run with:
#   source(system.file("examples", "augmentation.R", package = "SynthesizerPlus"))
#
# The rules that make this honest:
#   1. split the REAL data into training and test sets first;
#   2. fit the synthesizer on the TRAINING data only;
#   3. evaluate every model on the held-out REAL test data;
#   4. repeat over several random splits - one split proves nothing.
#
# Requires the recommended packages 'rpart' and 'MASS' (shipped with R).

library(SynthesizerPlus)
library(rpart)

rmse <- function(pred, obs) sqrt(mean((pred - obs)^2))
tree <- function(formula, data) {
  rpart(formula, data, control = rpart.control(cp = 0.005, minsplit = 5))
}

# =============================================================================
# Part 1. Tabular data: only 60 observations to learn house prices
# =============================================================================
data(Boston, package = "MASS")

one_split <- function(seed) {
  set.seed(seed)
  idx   <- sample(nrow(Boston), 60)
  train <- Boston[idx, ]                     # the small data set we "have"
  test  <- Boston[-idx, ]                    # real data never seen in training

  # real training data + 2000 synthetic rows, synthesizer fitted on train only
  big <- augment_data(train, n = 2000, seed = seed)

  c(tree_real      = rmse(predict(tree(medv ~ ., train), test), test$medv),
    tree_augmented = rmse(predict(tree(medv ~ ., big), test), test$medv),
    lm_real        = rmse(predict(lm(medv ~ ., train), test), test$medv),
    lm_augmented   = rmse(predict(lm(medv ~ ., big), test), test$medv))
}

res_tab <- t(sapply(1:10, one_split))
cat("\nTest RMSE over 10 random splits (60 training rows):\n")
print(round(colMeans(res_tab), 3))
cat("Share of splits where augmentation helped the tree:",
    mean(res_tab[, "tree_augmented"] < res_tab[, "tree_real"]), "\n")

# Check the synthetic training data look like the real training data
set.seed(1)
train <- Boston[sample(nrow(Boston), 60), ]
syn   <- generate(fit_synthesizer(train), n = 2000, seed = 1)
print(plot_marginals(train, syn, vars = c("medv", "lstat", "rm", "crim")))

# =============================================================================
# Part 2. Time series: forecast with a short history
# =============================================================================
# Only the first 4 years (48 months) of ldeaths are used for training;
# the remaining 24 months are the real test period.
lags    <- 12
n_train <- 48
train_ts <- window(ldeaths, end = time(ldeaths)[n_train])

# windows on the full series; the test windows are those whose target lies
# in the test period
all_w  <- ts_windows(ldeaths, lags = lags)
test_w <- all_w[(n_train - lags + 1):nrow(all_w), ]

forecast_rmse <- function(seed) {
  # 1 real + 50 synthetic training series. For monthly data with a short
  # history, fix the VAR order at the seasonal period (AIC picks too few lags).
  series  <- augment_data(train_ts, nsim = 50, order = 12, seed = seed)
  real_w  <- ts_windows(series[1], lags = lags)
  aug_w   <- ts_windows(series, lags = lags)
  real_w$series <- NULL
  aug_w$series  <- NULL
  c(tree_real      = rmse(predict(tree(y ~ ., real_w), test_w), test_w$y),
    tree_augmented = rmse(predict(tree(y ~ ., aug_w), test_w), test_w$y),
    lm_real        = rmse(predict(lm(y ~ ., real_w), test_w), test_w$y))
}

res_ts <- t(sapply(1:10, forecast_rmse))
cat("\nOne-step forecast RMSE on the 24 test months (10 repetitions):\n")
print(round(colMeans(res_ts), 1))

# Look at a few synthetic training series next to the real one
fit_ts <- fit_synthesizer(train_ts, order = 12)
print(plot_ts(train_ts, generate(fit_ts, seed = 3)))
print(plot_acf(train_ts, generate(fit_ts, n = 2000, seed = 3), lag_max = 24))

# =============================================================================
# What to take away
# =============================================================================
# * Augmentation helps flexible learners (here: regression trees) trained on
#   small samples, because synthetic data smooth the training distribution.
# * It does not help a linear model, which is already smooth, and a
#   well-specified simple model on the real data can still be the best choice.
# * The synthesizer must reproduce the structure the model relies on: for the
#   seasonal series the VAR order had to cover the seasonal period.
# * Synthetic data add no new information: never evaluate on synthetic data,
#   and never fit the synthesizer on data that include the test set.
