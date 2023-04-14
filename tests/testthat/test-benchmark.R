chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
if (nzchar(chk) && chk == "TRUE") {
  # use 2 cores in CRAN/Travis/AppVeyor
  num_cores <- 2L
} else {
  # use all cores in devtools::test()
  num_cores <- parallel::detectCores()
}
library(dplyr)



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
    "flexmix" = list(name_fonction = em_flexmix_univariate, list_params = list())
  )

  univariate_distribution_parameters <- benchmark_univariate_GMM_estimation(
    mixture_functions = relevant_mixture_functions[2],
    initialisation_algorithms = c("kmeans"),
    sigma_values = list("low OVL" = rep(0.3, 2)),
    mean_values = list(c(0, 4)),
    proportions = list("small imbalanced" = c(0.8, 0.2)),
    prop_outliers = c(0), cores = num_cores,
    Nbootstrap = 2, nobservations = c(100)
  )

  # saveRDS(univariate_distribution_parameters,
  #         test_path("results", "univariate_test_distribution.rds"))
  expect_equal(univariate_distribution_parameters, readRDS(test_path("results", "univariate_test_distribution.rds")))
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
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list())
  )


  corr_sequence <- seq(-0.8, 0.8, 0.2)
  sigma_values <- list()
  for (corr_1 in corr_sequence) {
    for (corr_2 in corr_sequence) {
      sigma_values[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
        array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
    }
  }

  multivariate_distribution_parameters <- benchmark_multivariate_GMM_estimation(
    mixture_functions = relevant_mixture_functions[2:3],
    initialisation_algorithms = c("kmeans"), cores = num_cores,
    sigma_values = sigma_values[1:2],
    mean_values = list("small OVL" = matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
    proportions = list("balanced" = c(0.5, 0.5)),
    Nbootstrap = 4, nobservations = c(100)
  )

  # saveRDS(multivariate_distribution_parameters,
  #         test_path("results", "multivariate_test_distribution.rds"))
  original_dist <- readRDS(test_path("results", "multivariate_test_distribution.rds"))
  expect_equal(multivariate_distribution_parameters, original_dist)
})




test_that("computation time in multivariate", {
  skip_on_cran()
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
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list())
  )


  corr_sequence <- seq(-0.8, 0.8, 0.2)
  sigma_values <- list()
  for (corr_1 in corr_sequence) {
    for (corr_2 in corr_sequence) {
      sigma_values[[glue::glue("comp_1_corr_{corr_1}_comp_2_{corr_2}")]] <-
        array(c(1, corr_1, corr_1, 1, 1, corr_2, corr_2, 1), dim = c(2, 2, 2))
    }
  }

  multivariate_time_computations <- compute_microbenchmark_multivariate(
    mixture_functions = relevant_mixture_functions[2:3],
    initialisation_algorithms = c("kmeans"),
    sigma_values = sigma_values[1:2],
    mean_values = list("small OVL" = matrix(c(0, 2, 2, 0), nrow = 2, ncol = 2)),
    proportions = list("balanced" = c(0.5, 0.5)),
    Nbootstrap = 4, nobservations = c(100, 200, 500), cores = num_cores
  )

  # saveRDS(multivariate_time_computations,
  #         test_path("results", "multivariate_test_time.rds"))
  # expect_equal(multivariate_time_computations, readRDS(test_path("results", "multivariate_test_time.rds")))
})

test_that("univariate time computation", {
  skip_on_cran()
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
    "flexmix" = list(name_fonction = em_flexmix_univariate, list_params = list())
  )

  univariate_time_computations <- compute_microbenchmark_univariate(
    mixture_functions = relevant_mixture_functions[2:3],
    initialisation_algorithms = c("kmeans"),
    sigma_values = list("low OVL" = rep(0.3, 2)),
    mean_values = list(c(0, 4)),
    proportions = list(
      "balanced" = c(0.5, 0.5),
      "small imbalanced" = c(0.8, 0.2)
    ),
    prop_outliers = c(0), cores = num_cores,
    Nbootstrap = 4, nobservations = c(100, 200)
  )

  # saveRDS(univariate_time_computations,
  #         test_path("results", "univariate_test_time.rds"))
  expect_equal(univariate_time_computations, readRDS(test_path("results", "univariate_test_time.rds")))
})



test_that("HD multivariate benchmark", {

  k <- 2
  relevant_mixture_functions <- list(
    "EMMIXmfa" = list(name_fonction = em_EMMIXmfa_multivariate, list_params = list()),
    "HDclassif" = list(name_fonction = em_HDclassif_multivariate, list_params = list()),
    "em R" = list(name_fonction = emnmix_multivariate, list_params = list()),
    "Rmixmod" = list(name_fonction = RGMMBench::em_Rmixmod_multivariate, list_params = list()),
    "mixtools" = list(name_fonction = em_mixtools_multivariate, list_params = list()),
    "bgmm" = list(name_fonction = em_bgmm_multivariate, list_params = list()),
    "mclust" = list(name_fonction = em_mclust_multivariate, list_params = list(prior = NULL)),
    "EMCluster" = list(name_fonction = em_EMCluster_multivariate, list_params = list()),
    "GMKMcharlie" = list(name_fonction = em_GMKMcharlie_multivariate, list_params = list()),
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list()))

  #################################################################
  ##                         low overlap                         ##
  #################################################################
  RNGkind("L'Ecuyer-CMRG"); set.seed(20)
  theta_low_OVL_balanced_circular <- MixSim::MixSim(BarOmega = 10^-4,
                                                    K=2, p=10, sph = TRUE, hom = FALSE,
                                                    ecc = 0.9, PiLow = 1, int = c(0.0, 2.0))
  theta_low_OVL_balanced_circular_formatted <- list(p=theta_low_OVL_balanced_circular$Pi,
                                                    mu=t(theta_low_OVL_balanced_circular$Mu),
                                                    sigma=theta_low_OVL_balanced_circular$S)


  HD_low_OVL_balanced_circular_distribution_parameters <- benchmark_multivariate_GMM_estimation(
    id_scenario = 1,
    initialisation_algorithms = c("random"),
    mixture_functions = list("em R" = list(name_fonction = emnmix_multivariate, list_params = list())),
    sigma_values = list(theta_low_OVL_balanced_circular_formatted$sigma),
    mean_values = list(theta_low_OVL_balanced_circular_formatted$mu),
    proportions = list(theta_low_OVL_balanced_circular_formatted$p),
    Nbootstrap = 1, nobservations = c(2000))


  multivariate_time_computations <- compute_microbenchmark_multivariate (
    mixture_functions = list("HDclassif" = list(name_fonction = em_HDclassif_multivariate, list_params = list()),
                             "em R" = list(name_fonction = emnmix_multivariate, list_params = list())),
    sigma_values = list(theta_low_OVL_balanced_eccentric_formatted$sigma),
    mean_values = list(theta_low_OVL_balanced_eccentric_formatted$mu),
    proportions = list(theta_low_OVL_balanced_eccentric_formatted$p, c(0.9, 0.1)),
    nobservations = c(100, 200, 500),
    Nbootstrap = 5, epsilon = 10^-4, itmax = 5,
    initialisation_algorithms = c("random", "rebmix"))






  RNGkind("L'Ecuyer-CMRG"); set.seed(20)
  relevant_mixture_functions <- list(
    "HDclassif" = list(name_fonction = em_HDclassif_multivariate, list_params = list()),
    "mixtools" = list(name_fonction = em_mixtools_multivariate, list_params = list()),
    "bgmm" = list(name_fonction = em_bgmm_multivariate, list_params = list()),
    "EMCluster" = list(name_fonction = em_EMCluster_multivariate, list_params = list()),
    "flexmix" = list(name_fonction = em_flexmix_multivariate, list_params = list()))

  theta_high_OVL_unbalanced_circular <- readRDS("../mixture_models/theta_high_OVL_unbalanced_circular.rds")[[1]]

  HD_high_OVL_unbalanced_circular_distribution_parameters <- benchmark_multivariate_GMM_estimation(
    id_scenario = 5, initialisation_algorithms = c("kmeans", "rebmix", "hc"),
    mixture_functions = relevant_mixture_functions,
    sigma_values = list(theta_high_OVL_unbalanced_circular$sigma),
    mean_values = list(theta_high_OVL_unbalanced_circular$mu),
    proportions = list(theta_high_OVL_unbalanced_circular$p),
    Nbootstrap = 30, nobservations = c(2000)) %>% magrittr::extract2("distributions")



  plot_boxplots_parameters(HD_high_OVL_unbalanced_circular_distribution_parameters %>% filter(ID == "5a"),
                           num_col = 2, true_theta = enforce_identifiability(theta_high_OVL_unbalanced_circular),
                           match_symbol = "^p[[:digit:]]+$|mu_var1_|sd_var1_var1|sd_var2_var1|sd_var2_var2")

  plot_correlation_Heatmap(HD_high_OVL_unbalanced_circular_distribution_parameters %>% filter(ID == "5a"))$kmeans




})









