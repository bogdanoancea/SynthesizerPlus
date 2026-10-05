# Simulated census-like microdata used in examples and vignettes.
# The data are generated, not real: no confidential information is involved.
library(SynthesizerPlus)

regions <- c("Nord-Est", "Sud-Est", "Sud-Muntenia", "Sud-Vest Oltenia", "Vest",
             "Nord-Vest", "Centru", "Bucuresti-Ilfov")
reg_p <- c(0.17, 0.13, 0.15, 0.10, 0.09, 0.13, 0.12, 0.11)

R <- matrix(c(
  # age   sex  region  educ  empl  hhsize income
  1.00,  0.05, 0.00, -0.25,  0.55, -0.30,  0.15,
  0.05,  1.00, 0.00,  0.05,  0.15,  0.00, -0.20,
  0.00,  0.00, 1.00,  0.20, -0.05, -0.10,  0.25,
 -0.25,  0.05, 0.20,  1.00, -0.30, -0.15,  0.55,
  0.55,  0.15,-0.05, -0.30,  1.00,  0.05, -0.45,
 -0.30,  0.00,-0.10, -0.15,  0.05,  1.00, -0.05,
  0.15, -0.20, 0.25,  0.55, -0.45, -0.05,  1.00), 7, 7)

census_sim <- r_mvdist(
  5000,
  margins = list(
    age = function(p) as.integer(round(qbeta(p, 1.6, 1.9) * 72 + 18)),
    sex = margin_categorical(c("male", "female"), c(0.48, 0.52)),
    region = margin_categorical(regions, reg_p),
    education = margin_categorical(c("primary", "lower secondary", "upper secondary", "tertiary"),
                                   c(0.15, 0.25, 0.40, 0.20), ordered = TRUE),
    employment = margin_categorical(c("employed", "unemployed", "inactive", "retired"),
                                    c(0.50, 0.05, 0.17, 0.28)),
    household_size = function(p) as.integer(qpois(p, 1.6) + 1L),
    income = function(p) round(qlnorm(p, meanlog = 8.1, sdlog = 0.6), -1)
  ),
  copula = "gaussian", corr = R, seed = 2026
)
census_sim$income[census_sim$employment %in% c("unemployed", "inactive")] <-
  round(census_sim$income[census_sim$employment %in% c("unemployed", "inactive")] * 0.25, -1)
census_sim$employment[census_sim$age < 20 & census_sim$employment == "retired"] <- "inactive"
census_sim$person_id <- sprintf("P%05d", seq_len(nrow(census_sim)))
census_sim <- census_sim[c("person_id", setdiff(names(census_sim), "person_id"))]

save(census_sim, file = "data/census_sim.rda", compress = "xz")
