test_that("univariate benchmark", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20)
  relevant_mixture_functions <- list(
    "RGMMBench" = list(name_fonction = emnmix_univariate, list_params = list()),
    "Rmixmod" = list(name_fonction = RGMMBench::em_Rmixmod_univariate, list_params = list()),
    "mixtools" = list(name_fonction = em_mixtools_univariate, list_params = list()),
    "bgmm" = list(name_fonction = em_bgmm_univariate, list_params = list()),
    "mclust" = list(name_fonction = em_mclust_univariate, list_params = list(prior = NULL)),
    "EMCluster" = list(name_fonction = em_EMCluster_univariate, list_params = list()),
    "GMKMcharlie" = list(name_fonction = em_GMKMcharlie_univariate, list_params = list()),
    "flexmix" = list(name_fonction = em_flexmix_univariate, list_params = list()))

  test_parameters_distribution <- benchmark_univariate_GMM_estimation(mixture_functions=relevant_mixture_functions[2],
                                                                      initialisation_algorithms = c("kmeans"),
                                                                      sigma_values=list("low OVL"= rep(0.3, 2)),
                                                                      mean_values=list(c(0, 4)),
                                                                      proportions = list("small imbalanced"=c(0.8, 0.2)),
                                                                      prop_outliers = c(0),
                                                                      Nbootstrap=5,  nobservations=c(100))

  saveRDS(test_parameters_distribution,
          file.path("./results", "univariate_test_distribution.rds"))
})



test_that("multivariate benchmark", {
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
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list()))


  corr_sequence <- seq(-0.8, 0.8, 0.2)
  sigma_values <- list()
  for (corr_1 in corr_sequence) {
    for (corr_2 in corr_sequence) {
      sigma_values[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
        array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
    }
  }

  test_parameters_distribution <- benchmark_multivariate_GMM_estimation(
    mixture_functions = relevant_mixture_functions[2:3],
    initialisation_algorithms = c("kmeans"),
    sigma_values = sigma_values[1:2],
    mean_values = list("small OVL"=matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
    proportions = list("balanced"=c(0.5, 0.5)),
    Nbootstrap = 4, nobservations = c(100))

  saveRDS(test_parameters_distribution,
          file.path("./results", "multivariate_test_distribution.rds"))
})




test_that("computation time in multivariate", {
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
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list()))


  corr_sequence <- seq(-0.8, 0.8, 0.2)
  sigma_values <- list()
  for (corr_1 in corr_sequence) {
    for (corr_2 in corr_sequence) {
      sigma_values[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
        array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
    }
  }

  test_parameters_distribution <- compute_microbenchmark_multivariate(
    mixture_functions = relevant_mixture_functions[2:3],
    initialisation_algorithms = c("kmeans"),
    sigma_values = sigma_values[1:2],
    mean_values = list("small OVL"=matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
    proportions = list("balanced"=c(0.5, 0.5)),
    Nbootstrap = 4, nobservations = c(100, 200, 500))

  saveRDS(test_parameters_distribution,
          file.path("./results", "multivariate_test_time_computations.rds"))
})

test_that("univariate time computation", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20)
  relevant_mixture_functions <- list(
    "RGMMBench" = list(name_fonction = emnmix_univariate, list_params = list()),
    "Rmixmod" = list(name_fonction = RGMMBench::em_Rmixmod_univariate, list_params = list()),
    "mixtools" = list(name_fonction = em_mixtools_univariate, list_params = list()),
    "bgmm" = list(name_fonction = em_bgmm_univariate, list_params = list()),
    "mclust" = list(name_fonction = em_mclust_univariate, list_params = list(prior = NULL)),
    "EMCluster" = list(name_fonction = em_EMCluster_univariate, list_params = list()),
    "GMKMcharlie" = list(name_fonction = em_GMKMcharlie_univariate, list_params = list()),
    "flexmix" = list(name_fonction = em_flexmix_univariate, list_params = list()))

  test_parameters_distribution <- compute_microbenchmark_univariate(mixture_functions=relevant_mixture_functions[2:3],
                                                                      initialisation_algorithms = c("kmeans", "quantiles"),
                                                                      sigma_values=list("low OVL"= rep(0.3, 2)),
                                                                      mean_values=list(c(0, 4)),
                                                                      proportions = list( "balanced"=c(0.5, 0.5),
                                                                                          "small imbalanced"=c(0.8, 0.2)),
                                                                      prop_outliers = c(0),
                                                                      Nbootstrap=5,  nobservations=c(100, 200))

  saveRDS(test_parameters_distribution,
          file.path("./results", "univariate_test_time_computations.rds"))
})



