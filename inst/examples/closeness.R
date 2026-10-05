# SynthesizerPlus: controlling how close synthetic data are to the real data
# Run with:
#   source(system.file("examples", "closeness.R", package = "SynthesizerPlus"), echo = TRUE)

library(SynthesizerPlus)

fit <- fit_synthesizer(mtcars)

# 1. One knob: closeness (1 = as close as possible, 0 = far) -------------------
for (cl in c(1, 0.75, 0.5, 0.25, 0)) {
  s <- generate(fit, closeness = cl, seed = 1)
  cat(sprintf("closeness %.2f | discriminator AUC %.2f | DCR ratio %.2f | near-copies %.3f\n",
              cl, discriminator_auc(mtcars, s, seed = 1),
              dcr(mtcars, s, seed = 1)$ratio, dcr(mtcars, s, seed = 1)$close_share))
}

# 2. Record level: further from real records, same distributions ---------------
#    noise: smoothing of numeric variables (mean and sd preserved)
#    swap : random re-drawing of categorical / discrete values
p_rec <- perturbation(noise = 0.5, swap = 0.2)
p_rec
s_rec <- generate(fit, n = 5000, perturb = p_rec, seed = 1)
round(rbind(real = sapply(mtcars[c("mpg", "hp", "wt")], sd),
            perturbed = sapply(s_rec[c("mpg", "hp", "wt")], sd)), 2)

# 3. Distribution level: scenarios --------------------------------------------
#    shift (in sd), scale (spread), temperature (category shares), dependence
scenario <- perturbation(shift = c(hp = 1), scale = c(wt = 1.5),
                         temperature = c(cyl = 3), dependence = 1.5)
s_scen <- generate(fit, n = 5000, perturb = scenario, seed = 1)
c(real_hp = mean(mtcars$hp), scenario_hp = mean(s_scen$hp))
round(rbind(real = prop.table(table(mtcars$cyl)),
            scenario = prop.table(table(s_scen$cyl))), 2)
print(plot_marginals(mtcars, s_scen, vars = c("hp", "wt", "cyl")))

# keep values within the observed range, or allow any value
s_obs <- generate(fit, perturb = perturbation(noise = 1, bounds = "observed"), seed = 1)
range(s_obs$mpg); range(mtcars$mpg)

# 4. Knob and perturbation together --------------------------------------------
s_both <- generate(fit, closeness = 0.8, perturb = perturbation(shift = c(mpg = -1)), seed = 1)

# 5. Hit a target instead of guessing --------------------------------------------
cal <- calibrate_closeness(fit, mtcars, target = 0.85, metric = "auc",
                           grid = seq(0, 1, by = 0.25), refine = 3)
cal
print(plot(cal))
s_cal <- generate(fit, closeness = cal$closeness, seed = 7)

# other built-in targets: "dcr" (distance to closest record), "pmse", "marginal"
cal2 <- calibrate_closeness(fit, mtcars, target = 4, metric = "dcr",
                            grid = seq(0, 1, by = 0.25), refine = 2)
cal2

# ... or any criterion of your own, e.g. the mean correlation error
cor_error <- function(real, synthetic) {
  mean(abs(cor(synthetic[c("mpg", "hp", "wt")]) - cor(real[c("mpg", "hp", "wt")])))
}
cal3 <- calibrate_closeness(fit, mtcars, target = 0.2, metric = cor_error,
                            grid = seq(0, 1, by = 0.25), refine = 2)
cal3

# 6. The same controls in the other entry points --------------------------------
synthesize(mtcars, n = 10, closeness = 0.6, seed = 1)
augment_data(mtcars, n = 100, closeness = 0.8, seed = 1)        # diverse training data
sims <- simulate(fit, nsim = 3, closeness = 0.7, seed = 1)       # several noisy copies
out <- tempfile(fileext = ".csv")
generate_to_file(fit, n = 1000, path = out, closeness = 0.7, seed = 1)
