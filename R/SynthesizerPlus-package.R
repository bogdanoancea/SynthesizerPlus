#' @keywords internal
#' @section Main functions:
#' * Fitting and generation: [fit_synthesizer()], [generate()],
#'   [synthesize()], [simulate.sp_synthesizer()], [generate_to_file()].
#' * Disclosure control: [disclosure_risk()], [dcr()], [census_sim] data.
#' * Closeness to the real data: `closeness` in [generate()],
#'   [perturbation()], [calibrate_closeness()].
#' * Time series: [fit_synthesizer.ts()].
#' * Machine-learning augmentation: [augment_data()], [ts_windows()].
#' * Multivariate distributions: [r_mvnorm()], [r_mvt()], [r_mvskewnorm()],
#'   [r_dirichlet()], [r_mvmixture()], [r_copula()], [r_mvdist()],
#'   [make_corr()].
#' * Evaluation: [compare_synthetic()], [marginal_metrics()],
#'   [association_matrix()], [pmse()], [discriminator_auc()],
#'   [ci_overlap()], [dcr()].
#' * Visualisation: [plot_marginals()], [plot_association()], [plot_pairs()],
#'   [plot_qq()], [plot_ts()], [plot_acf()], [plot_dcr()].
#' * Input/output: [read_data()], [write_data()], [supported_formats()],
#'   [match_types()].
#' @importFrom ggplot2 .data
"_PACKAGE"
