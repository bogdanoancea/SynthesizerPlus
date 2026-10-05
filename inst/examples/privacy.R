# SynthesizerPlus: releasing confidential microdata as synthetic data
# Run with:
#   source(system.file("examples", "privacy.R", package = "SynthesizerPlus"), echo = TRUE)
#
# census_sim is a SIMULATED census-like file (no real persons).

library(SynthesizerPlus)

keys      <- c("age", "sex", "region", "household_size")  # known to an intruder
sensitive <- c("income", "education")                     # confidential
model     <- income ~ age + sex + education + employment  # a typical analysis

# 1. Risk of releasing the real microdata ------------------------------------
disclosure_risk(census_sim, census_sim, keys = keys, target = sensitive,
                ignore = "person_id")

# 2. Synthetic file -------------------------------------------------------------
fit <- fit_synthesizer(census_sim, id_cols = "person_id",
                       by = c("employment", "sex"), knots = 500)
syn <- generate(fit, seed = 1)
disclosure_risk(census_sim, syn, keys = keys, target = sensitive,
                ignore = "person_id")
dcr(census_sim[-1], syn[-1], max_rows = 1500, seed = 1)$close_share  # ~0.05 = chance

# 3. How much noise? Noise level = 1 - closeness ------------------------------
for (noise in c(0, 0.1, 0.2, 0.4)) {
  s <- generate(fit, closeness = 1 - noise, seed = 1)
  cat(sprintf("noise %.1f: near-copies %.3f | income CAP %.3f | utility %.2f\n",
              noise,
              dcr(census_sim[-1], s[-1], max_rows = 1500, seed = 1)$close_share,
              disclosure_risk(census_sim, s, keys, "income")$attribute$cap,
              mean(pmax(ci_overlap(census_sim, s, model)$overlap, 0))))
}

# 4. Noise targeted at the variables that need it --------------------------------
s_t <- generate(fit, seed = 1,
                perturb = perturbation(noise = c(income = 0.5, age = 0.3),
                                       swap = c(region = 0.2)))
disclosure_risk(census_sim, s_t, keys = keys, target = sensitive,
                ignore = "person_id")

# 5. Choose the noise from a risk rule: at most 4.5% near-copies -----------------
near_copies <- function(real, synthetic) {
  dcr(real[-1], synthetic[-1], max_rows = 1500, seed = 1)$close_share
}
cal <- calibrate_closeness(fit, census_sim, target = 0.045, metric = near_copies,
                           grid = c(0.6, 0.8, 1), refine = 2)
cal
cat("noise level to report:", round(1 - cal$closeness, 3), "\n")

release <- generate(fit, closeness = cal$closeness, seed = 2026)
ci_overlap(census_sim, release, model)[, c("term", "estimate_real",
                                           "estimate_synthetic", "overlap")]
print(plot_marginals(census_sim[-1], release[-1],
                     vars = c("age", "income", "household_size")))

# Release `release` - never the fitted model `fit`.
