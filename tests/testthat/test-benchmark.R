RNGkind("L'Ecuyer-CMRG")
set.seed(20)

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

# simulated_distribution <- readRDS("./errors/package_em R_init_algo_rebmix.rds")$simulated_distribution

test_parameters_distribution <- benchmark_multivariate_GMM_estimation_parallel(
  mixture_functions = relevant_mixture_functions,
  initialisation_algorithms = c("kmeans", "hc"),
  sigma_values = sigma_values[1],
  mean_values = list("high OVL"=matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
  proportions = list("balanced"=c(0.5, 0.5)),
  Nbootstrap = 10, nobservations = c(100))

# saveRDS(test_parameters_distribution,
#         file.path("/home/bncl_cb/rstudio/working/mixture_models/results", "multivariate_test_distribution.rds"))

multi_parameters_distribution <- benchmark_multivariate_GMM_estimation(
  mixture_functions = relevant_mixture_functions,
  initialisation_algorithms = c("kmeans", "random", "hc", "small em", "rebmix"),
  sigma_values = sigma_values,
  mean_values = list("small OVL"=matrix(c(20, 40, 40, 20), nrow = 2, ncol = 2)),
  # mean_values = list("high OVL"=matrix(c(20, 22, 22, 20), nrow = 2, ncol = 2),
  #                    "small OVL"=matrix(c(20, 40, 40, 20), nrow = 2, ncol = 2)),
  proportions = list("balanced" = c(0.5, 0.5)),
  Nbootstrap = 100, nobservations = c(1000))




#
# capture.output(running_times <- compute_microbenchmark(
#   mixture_functions = relevant_mixture_functions, sigma_values = list("null OVL" = rep(0.3, 4)),
#   mean_values = list(c(0, 4)), proportions = list("balanced" = c(0.5, 0.5)), prop_outliers = c(0),
#   Nbootstrap = 2, nobservations = c(50, 100), initialisation_algorithms = c("quantiles")
# )) %>% invisible()


#
# density_plots <- plot_density_distribution(
#   sigma_values = list("null OVL" = c(0.3, 0.3), "small OVL" = c(1, 1), "high OVL" = c(2, 2)),
#   mean_values = list(c(0, 4)), skewness_values = list("null skewness" = c(0, 0), "high skewness" = c(20, 20)),
#   proportions = list("balanced" = c(0.5, 0.5), "little unbalanced" = c(0.65, 0.35), "highly unbalanced" = c(0.9, 0.1))
# )
