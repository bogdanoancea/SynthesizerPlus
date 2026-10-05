#' Simulated census-like microdata
#'
#' A **simulated** data set of 5000 persons with the structure of census or
#' household-survey microdata, used to illustrate synthetic data for
#' statistical disclosure control. The records were generated with
#' [r_mvdist()] (Gaussian copula with plausible correlations, see
#' `data-raw/census_sim.R` in the package sources); they do not describe real
#' people. Region names are the eight development regions of Romania; their
#' shares only roughly follow the population distribution.
#'
#' @format A data frame with 5000 rows and 8 columns:
#' \describe{
#'   \item{person_id}{Direct identifier (character).}
#'   \item{age}{Age in years, 18--90 (integer).}
#'   \item{sex}{Factor: male, female.}
#'   \item{region}{Factor: development region (8 levels).}
#'   \item{education}{Ordered factor: primary < lower secondary < upper
#'     secondary < tertiary.}
#'   \item{employment}{Factor: employed, unemployed, inactive, retired.}
#'   \item{household_size}{Number of persons in the household (integer).}
#'   \item{income}{Monthly income in national currency (numeric).}
#' }
#' @source Simulated with SynthesizerPlus.
#' @examples
#' str(census_sim)
#' table(census_sim$employment, census_sim$sex)
"census_sim"
