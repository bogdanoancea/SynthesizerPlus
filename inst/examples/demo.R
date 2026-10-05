# SynthesizerPlus: a short tour
# Run with: source(system.file("examples", "demo.R", package = "SynthesizerPlus"))

library(SynthesizerPlus)

# 1. Fit a synthesizer -------------------------------------------------------
# airquality has numeric columns and missing values; stratify by month
fit <- fit_synthesizer(airquality, by = "Month")
fit
summary(fit)

# 2. Generate synthetic data -------------------------------------------------
syn <- generate(fit, n = 1000, seed = 123)
head(syn)
colMeans(is.na(syn))          # missingness is reproduced
colMeans(is.na(airquality))

# one-step alternative, any data frame with mixed column types
syn_iris <- synthesize(iris, n = 500, by = "Species", seed = 1)

# weaker dependence between variables (marginals unchanged)
syn_weak <- generate(fit, n = 1000, dependence = 0.5, seed = 1)

# how close to the real data: 1 = as close as possible, 0 = far
syn_far <- generate(fit, n = 1000, closeness = 0.6, seed = 1)

# fine control: record-level noise/swap, distribution-level shift/scale/...
syn_scen <- generate(fit, n = 1000, seed = 1,
                     perturb = perturbation(noise = 0.3, shift = c(Temp = 1)))

# find the closeness that gives a target distance (here: DCR ratio = 1.3)
cal <- calibrate_closeness(fit, airquality, target = 1.3, metric = "dcr",
                           grid = seq(0, 1, 0.25), refine = 2)
cal

# 3. Evaluate quality and disclosure risk ------------------------------------
cmp <- compare_synthetic(airquality, syn, seed = 1)
cmp
summary(cmp)

# a regression gives the same answers on real and synthetic data?
ci_overlap(airquality, syn, Ozone ~ Temp + Wind)

# 4. Visualise ---------------------------------------------------------------
print(plot_marginals(airquality, syn))
print(plot_association(airquality, syn, vars = c("Ozone", "Solar.R", "Wind", "Temp")))
print(plot_pairs(iris, syn_iris, vars = c("Sepal.Length", "Petal.Length", "Petal.Width")))
print(plot_dcr(cmp))

# 5. Time series ---------------------------------------------------------------
fit_ts <- fit_synthesizer(ldeaths)        # copula-VAR by default
syn_ts <- generate(fit_ts, n = 120, seed = 1)
print(plot_ts(ldeaths, syn_ts))
print(plot_acf(ldeaths, syn_ts))

# 6. Multivariate distributions --------------------------------------------------
x <- r_mvnorm(1000, mean = c(0, 0, 0), sigma = make_corr(3, "ar1", 0.7), seed = 1)
round(cor(x), 2)

u <- r_copula(1000, "clayton", dim = 3, theta = 2, seed = 1)
print(plot_pairs(as.data.frame(u)))

sim <- r_mvdist(
  1000,
  margins = list(
    income = list(dist = "lnorm", meanlog = 10, sdlog = 0.5),
    age    = list(dist = "unif", min = 18, max = 80),
    sector = margin_categorical(c("agri", "industry", "services"), c(0.1, 0.3, 0.6))
  ),
  copula = "gaussian", corr = make_corr(3, "exchangeable", 0.4), seed = 1
)
str(sim)

# 7. Read and write many file formats ------------------------------------------
supported_formats(available_only = TRUE)

out <- tempfile("sp_demo_")
dir.create(out)
write_data(syn, file.path(out, "airquality_syn.csv"))
write_data(syn, file.path(out, "airquality_syn.rds"))
if (requireNamespace("arrow", quietly = TRUE) ||
    requireNamespace("nanoparquet", quietly = TRUE)) {
  write_data(syn, file.path(out, "airquality_syn.parquet"))
}
if (requireNamespace("hdf5r", quietly = TRUE)) {
  write_data(syn, file.path(out, "airquality_syn.h5"))
}
list.files(out)

# CSV does not store column types: restore them from a template
back <- read_data(file.path(out, "airquality_syn.csv"), template = airquality)
str(back)

# 8. Large data sets, written in chunks -----------------------------------------
res <- generate_to_file(fit, n = 100000, path = file.path(out, "big.csv.gz"),
                        chunk_size = 25000, seed = 1)
nrow(read_data(res$path))

unlink(out, recursive = TRUE)
