

relevant_mixture_functions <- list(
  "em R" = list(name_fonction = emnmix_multivariate, list_params = list()),
  "Rmixmod" = list(name_fonction = RGMMBench::em_Rmixmod_multivariate, list_params = list()),
  "mixtools" = list(name_fonction = em_mixtools_multivariate, list_params = list()),
  "bgmm" = list(name_fonction = em_bgmm_multivariate, list_params = list()),
  "mclust" = list(name_fonction = em_mclust_multivariate, list_params = list(prior = NULL)),
  "EMCluster" = list(name_fonction = em_EMCluster_multivariate, list_params = list()),
  "GMKMcharlie" = list(name_fonction = em_GMKMcharlie_multivariate, list_params = list()),
  "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list()),
  "DCEM" = list(name_fonction = em_DCEM_multivariate, list_params = list()))


corr_sequence <- seq(-0.8, 0.8, 0.2)
sigma_values <- list()
for (corr_1 in corr_sequence) {
  for (corr_2 in corr_sequence) {
    sigma_values[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
      array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
  }
}
# simulated_distribution <- readRDS("./errors/initialisation_failures/init_algo_rebmix_error_1.rds")
# simulated_distribution <- readRDS("./errors/estimation_failures/package_em R_init_algo_kmeans_step_16_100_observations.rds")$simulated_distribution
# initial_estimates <- readRDS("./errors/estimation_failures/package_em R_init_algo_kmeans_step_16_100_observations.rds")$initial_estimates
# RNGkind("L'Ecuyer-CMRG")
# set.seed(20)
# test_parameters_distribution <- benchmark_multivariate_GMM_estimation(
#   mixture_functions = relevant_mixture_functions[1],
#   initialisation_algorithms = c("kmeans"),
#   sigma_values = sigma_values[1:2],
#   mean_values = list("small OVL"=matrix(c(0, 20, 20, 0), nrow = 2, ncol = 2)),
#   proportions = list("balanced"=c(0.5, 0.5), "unbalanced"=c(0.9, 0.1)),
#   Nbootstrap = 6, nobservations = c(200))

# saveRDS(test_parameters_distribution,
#         file.path("/home/bncl_cb/rstudio/working/mixture_models/results", "multivariate_test_distribution.rds"))


sigma_values_relevant <- list()
for (corr_1 in c(-0.8, 0.8)) {
  for (corr_2 in c(-0.8, 0.8)) {
    sigma_values_relevant[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
      array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
  }
}
sigma_values_relevant[["well-separated"]] <- array(rep(c(1, 0, 0, 1), 2), dim = c(2, 2, 2))
time_computation <- RGMMBench::compute_microbenchmark_multivariate(mixture_functions = relevant_mixture_functions[1:2],
                                                                   initialisation_algorithms = c("kmeans"),
                                                                   sigma_values = sigma_values_relevant[1:2],
                                                                   mean_values = list("high OVL"=matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
                                                                   proportions = list("unbalanced"=c(0.9, 0.1)),
                                                                   Nbootstrap = 100, nobservations = c(100, 200))




# RNGkind("L'Ecuyer-CMRG"); set.seed(20)
# simulated_distribution <- readRDS("./errors/initialisation_failures/init_algo_rebmix_step_100_500_observations.rds")$simulated_distribution
# running_times_multivariate <- compute_microbenchmark_multivariate(
#   mixture_functions = relevant_mixture_functions[1:2],
#   initialisation_algorithms = c("kmeans", "hc", "random"),
#   sigma_values = sigma_values[1],
#   mean_values = list("high OVL"=matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
#   proportions = list("unbalanced"=c(0.9, 0.1)),
#   Nbootstrap = 4, nobservations = c(50, 100))


#
# density_plots <- plot_density_distribution(
#   sigma_values = list("null OVL" = c(0.3, 0.3), "small OVL" = c(1, 1), "high OVL" = c(2, 2)),
#   mean_values = list(c(0, 4)), skewness_values = list("null skewness" = c(0, 0), "high skewness" = c(20, 20)),
#   proportions = list("balanced" = c(0.5, 0.5), "little unbalanced" = c(0.65, 0.35), "highly unbalanced" = c(0.9, 0.1))
# )
