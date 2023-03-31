# grep -rl 'benchmark_distribution_parameters' *.R | xargs -i@ sed -i 's/benchmark_distribution_parameters/benchmark_univariate_GMM_estimation/g' @
test_that("simulation of univariate GMM", {
  set.seed(20)
  # define parameters of a two component multivariate GMM
  true_theta <- list(p = c(0.5, 0.5), mu = c(0, 2), sigma = c(1, 1))
  expect_true(check_parameters_validity_univariate(true_theta, k = 2))

  univariate_simulation <- simulate_univariate_GMM(theta = true_theta, n = 100)

  observed_estimates <- estimate_supervised_univariate_GMM(
    x = univariate_simulation$x,
    s = univariate_simulation$s
  )

  #################################################################
  ##          test initialisation in univariate context          ##
  #################################################################
  expect_warning(initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "hc"
  ))

  initial_estimates_kmeans <- initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "kmeans"
  )

  initial_estimates_quantiles <- initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "quantiles"
  )

  initial_estimates_rebmix <- initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "rebmix"
  )

  initial_estimates_small_em <- initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "small em"
  )

  initial_estimates_random <- initialize_em_univariate(
    x = univariate_simulation$x, k = 2,
    initialisation_algorithm = "random"
  )




  ##################################################################
  ##               estimation in univariate context               ##
  ##################################################################

  univariate_em_estimates <- emnmix_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_mclust <- em_mclust_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_bgmm <- em_bgmm_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_Rmixmod <- em_Rmixmod_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_mixtools <- em_mixtools_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_flexmix <- em_flexmix_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_EMCluster <- em_EMCluster_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  univariate_em_GMKMcharlie <- em_GMKMcharlie_univariate(univariate_simulation$x,
    k = 2,
    start = true_theta, iter = 20
  )

  expect_equal(univariate_em_mclust, univariate_em_Rmixmod, tolerance = 10^-3)
})


test_that("GMM estimation in multivariate case", {
  set.seed(20)

  # multivariate scenario
  true_theta <- list(
    p = c(0.5, 0.5),
    mu = matrix(c(20, 40, 40, 20), nrow = 2),
    sigma = array(c(1, 0.2, 0.2, 1, 1, -0.2, -0.2, 1), dim = c(2, 2, 2))
  )
  multivariate_simulation <- simulate_multivariate_GMM(theta = true_theta, n = 1000)
  inititial_kmeans_estimates <- initialize_em_multivariate(multivariate_simulation$x,
    k = 2,
    initialisation_algorithm = "kmeans"
  )

  # test all packages implementing GMM mixture
  own_EM_implementation <- emnmix_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6,
    start = inititial_kmeans_estimates
  )

  Rmixmod_multi_estimates <- em_Rmixmod_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1000, epsilon = 10^-12,
    start = inititial_kmeans_estimates
  )

  EMCluster_multi_estimates <- em_EMCluster_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6,
    start = inititial_kmeans_estimates
  )

  bgmm_multi_estimates <- em_bgmm_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6,
    start = inititial_kmeans_estimates
  )


  flexmix_multi_estimates <- em_flexmix_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6, minprior = 0,
    start = inititial_kmeans_estimates
  )

  mixtools_multi_estimates <- em_mixtools_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6,
    start = inititial_kmeans_estimates
  )

  mclust_multi_estimates <- em_mclust_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6,
    start = inititial_kmeans_estimates
  )

  GMKMcharlie_multi_estimates <- em_GMKMcharlie_multivariate(
    x = multivariate_simulation$x, k = 2,
    itmax = 1, epsilon = 10^-6, embedNoise = 0,
    start = inititial_kmeans_estimates, parallel = F
  )
  expect_equal(own_EM_implementation, flexmix_multi_estimates, tolerance = 10^-2) # GMKMCharlie as well
  expect_equal(own_EM_implementation, Rmixmod_multi_estimates)
})



test_that("GMM estimation in supervised case", {
  set.seed(20)
  true_theta <- list(
    p = c(0.2, 0.8),
    mu = matrix(c(20, 20, 20, 40, 40, 40), nrow = 3),
    sigma = array(rep(c(1, 0, 0.1, 0, 1, -0.1, 0.1, -0.1, 1), 2), dim = c(3, 3, 2))
  )

  multivariate_simulation <- simulate_multivariate_GMM(theta = true_theta, n = 2000)
  estimated_theta <- estimate_supervised_multivariate_GMM(multivariate_simulation$x, multivariate_simulation$s)
  # saveRDS(estimated_theta, test_path("fixtures", "two_component_3D_supervised_estimation.rds"))
  old_theta <- readRDS(test_path("fixtures", "two_component_3D_supervised_estimation.rds"))
  expect_equal(estimated_theta, old_theta)


  simulation_test_small_overlap <- MixSim::MixSim(BarOmega = 10^-4, MaxOmega = NULL, K=2, p=10, sph = FALSE, hom = FALSE,
                            ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), resN = 1000)

  simulation_test_strong_overlap <- MixSim::MixSim(BarOmega = 0.2, MaxOmega = NULL, K=2, p=10, sph = FALSE, hom = FALSE,
                                                  ecc = 0.90, PiLow = 1.0, int = c(0, 20), resN = 1000)

})

# small test
