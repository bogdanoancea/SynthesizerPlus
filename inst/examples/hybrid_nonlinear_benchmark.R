# Reproducible exploratory benchmark for the optional hybrid method.
# This script is NOT a formal proof of superiority; it prints replicate-level
# results and records settings. Run after installing this source package.
#
# Rscript inst/examples/hybrid_nonlinear_benchmark.R
suppressPackageStartupMessages(library(SynthesizerPlus))
if (!requireNamespace("rpart", quietly = TRUE)) stop("Install rpart")

make_data <- function(n, seed) {
  set.seed(seed)
  x1 <- runif(n, -2, 2)
  x2 <- rnorm(n)
  group <- factor(sample(c("A", "B", "C"), n, TRUE))
  mu <- sin(3.5 * x1) + 1.5 * (x2 > 0) +
    0.6 * x1 * x2 + 1 * (group == "C")
  data.frame(x1, x2, group, y = mu + rnorm(n, sd = 0.35))
}

rmse <- function(a, b) sqrt(mean((a - b)^2))
# The evaluation model is a tree fitted to synthetic data and tested on an
# independent population sample. It is not the same fitted tree as used by
# the hybrid generator.
model_error <- function(train, holdout) {
  tr <- rpart::rpart(y ~ x1 + x2 + group, data = train,
                     control = rpart::rpart.control(cp = 0.001,
                                                   minbucket = 15, xval = 0))
  rmse(holdout$y, as.numeric(predict(tr, newdata = holdout)))
}

results <- list()
for (rep in seq_len(20L)) {
  train <- make_data(1000, 1000L + rep)
  holdout <- make_data(3000, 3000L + rep)
  baseline <- fit_synthesizer(train, seed = 5000L + rep)
  hybrid <- fit_synthesizer(train, nonlinear_target = "y",
                            nonlinear_minbucket = 20L, nonlinear_cp = 0.001,
                            seed = 5000L + rep)
  synth_base <- generate(baseline, n = 1000L, seed = 6000L + rep)
  synth_hybrid <- generate(hybrid, n = 1000L, seed = 6000L + rep)
  results[[rep]] <- data.frame(
    rep = rep,
    method = c("original", "copula", "copula_CART_hybrid"),
    holdout_rmse = c(model_error(train, holdout),
                     model_error(synth_base, holdout),
                     model_error(synth_hybrid, holdout))
  )
}
res <- do.call(rbind, results)
print(aggregate(holdout_rmse ~ method, res,
                function(x) c(mean = mean(x), sd = sd(x))))
# Paired Monte Carlo differences, negative = lower downstream error for hybrid.
wide <- reshape(res, idvar = "rep", timevar = "method", direction = "wide")
diff <- wide$holdout_rmse.copula_CART_hybrid - wide$holdout_rmse.copula
se <- stats::sd(diff) / sqrt(length(diff))
interval <- mean(diff) + stats::qt(c(0.025, 0.975), df = length(diff) - 1) * se
cat(sprintf("Hybrid minus baseline RMSE: %.4f (95%% paired MC CI %.4f to %.4f)\n",
            mean(diff), interval[1], interval[2]))
cat(sprintf("Hybrid better in %d/%d Monte Carlo replicates\n",
            sum(diff < 0), length(diff)))
write.csv(res, "hybrid_nonlinear_benchmark.csv", row.names = FALSE)
cat("Replicate-level data: hybrid_nonlinear_benchmark.csv\n")
cat("Software environment:\n")
print(sessionInfo())
