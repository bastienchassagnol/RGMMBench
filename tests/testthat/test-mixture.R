# grep -rl 'benchmark_distribution_parameters' *.R | xargs -i@ sed -i 's/benchmark_distribution_parameters/benchmark_univariate_GMM_estimation/g' @
dir.create("./errors", showWarnings = F, recursive = T)
test_that("simulation of univariate or multivariate GMM", {
  set.seed(20)
  # define parameters of a two component multivariate GMM
  true_theta <- list(p=c(0.2, 0.8),
                     mu=matrix(c(20, 22, 22, 20), nrow = 2),
                     sigma=array(c(1, 0.2, 0.2, 1, 1, -0.2, -0.2, 1), dim=c(2, 2, 2)))

  expect_true(check_parameters_validity_multivariate(true_theta, k=2))

  multivariate_simulation <- simulate_multivariate_GMM (theta=true_theta, n=100)

  expect_equal(multivariate_simulation,
                   readRDS(test_path("fixtures", "two_component_multivariate_GMM.rds")))



  true_theta <- list(p=c(0.5, 0.5),
                     mu=matrix(c(20, 40, 40, 20), nrow = 2),
                     sigma=array(c(1, 0.2, 0.2, 1, 1, -0.2, -0.2, 1), dim=c(2, 2, 2)))
  multivariate_simulation <- simulate_multivariate_GMM (theta=true_theta, n=1000)

  # assign("failed_hc_counter", 0, envir=globalenv())
  inititial_hc_estimates <- initialize_em_multivariate(multivariate_simulation$x, k = 2,
                                                    initialisation_algorithm = "hc")

  inititial_small_EM_estimates <- initialize_em_multivariate(multivariate_simulation$x, k = 2,
                                                           initialisation_algorithm = "small em")

  Rmixmod_multi_estimates <- em_Rmixmod_multivariate (x = multivariate_simulation$x, k = 2,
                                                      itmax = 200, epsilon = 10^-6,
                                                      start = inititial_kmeans_estimates)

  expect_error(initialize_em_multivariate(multivariate_simulation$x, k = 2,
                                          initialisation_algorithm = "quantiles"))

})


test_that("GMM estimation in unsupervised case", {
  set.seed(20)

  # multivariate scenario
  true_theta <- list(p=c(0.5, 0.5),
                     mu=matrix(c(20, 40, 40, 20), nrow = 2),
                     sigma=array(c(1, 0.2, 0.2, 1, 1, -0.2, -0.2, 1), dim=c(2, 2, 2)))
  multivariate_simulation <- simulate_multivariate_GMM (theta=true_theta, n=1000)
  inititial_kmeans_estimates <- initialize_em_multivariate(multivariate_simulation$x, k = 2,
                                                             initialisation_algorithm = "kmeans")

  # test all packages implementing GMM mixture
  own_EM_implementation <- emnmix_multivariate (x = multivariate_simulation$x, k = 2,
                                                    itmax = 1, epsilon = 10^-6,
                                                    start = inititial_kmeans_estimates)

  Rmixmod_multi_estimates <- em_Rmixmod_multivariate (x = multivariate_simulation$x, k = 2,
                                                      itmax = 1000, epsilon = 10^-12,
                                                      start = inititial_kmeans_estimates)

  EMCluster_multi_estimates <- em_EMCluster_multivariate (x = multivariate_simulation$x, k = 2,
                                                      itmax = 1, epsilon = 10^-6,
                                                      start = inititial_kmeans_estimates)

  bgmm_multi_estimates <- em_bgmm_multivariate (x = multivariate_simulation$x, k = 2,
                                                          itmax = 1, epsilon = 10^-6,
                                                          start = inititial_kmeans_estimates)

  flexmix_multi_estimates <- em_flexmix_multivariate (x = multivariate_simulation$x, k = 2,
                                                itmax = 1, epsilon = 10^-6, minprior = 0,
                                                start = inititial_kmeans_estimates)
  #

  mixtools_multi_estimates <- em_mixtools_multivariate (x = multivariate_simulation$x, k = 2,
                                                      itmax = 1, epsilon = 10^-6,
                                                      start = inititial_kmeans_estimates)

  mclust_multi_estimates <- em_mclust_multivariate (x = multivariate_simulation$x, k = 2,
                                                itmax = 1, epsilon = 10^-6,
                                                start = inititial_kmeans_estimates)

  DCEM_multi_estimates <- em_DCEM_multivariate (x = multivariate_simulation$x, k = 2,
                                                itmax = 1, epsilon = 10^-6,
                                                start = inititial_kmeans_estimates)

  GMKMcharlie_multi_estimates <- em_GMKMcharlie_multivariate (x = multivariate_simulation$x, k = 2,
                                                        itmax = 1, epsilon = 10^-6, embedNoise = 0,
                                                        start = inititial_kmeans_estimates, parallel = F)



  expect_equal(own_EM_implementation, flexmix_multi_estimates) # GMKMCharlie as well



})

test_that("specific focus on my own developped function", {

  distribution <- readRDS("./errors/estimation_failures/package_em R_init_algo_kmeans_step_16_100_observations.rds")

  # initial_rebmix_estimates <- initialize_em_multivariate (x = distribution$simulated_distribution$x,
  #                                                  k = 2,initialisation_algorithm = "rebmix")
  #
  # GMKMcharlie_multi_estimates <- em_GMKMcharlie_multivariate (x = multivariate_simulation$x, k = 2,
  #                                                             itmax = 1, start = start, parallel = F)

  own_EM_implementation <- em_mixtools_multivariate (x = distribution$simulated_distribution$x, k = 2,
                                                itmax = 1000, epsilon = 10^-12, start = distribution$initial_estimates)


})


test_that("GMM estimation in supervised case", {
  set.seed(20)
  true_theta <- list(p=c(0.2, 0.8),
                     mu=matrix(c(20, 20, 20, 40, 40, 40), nrow = 3),
                     sigma=array(rep(c(1, 0, 0.1, 0, 1, -0.1, 0.1, -0.1, 1), 2), dim=c(3, 3, 2)))

  multivariate_simulation <- simulate_multivariate_GMM (theta=true_theta, n=2000)
  estimated_theta <- estimate_supervised_multivariate_GMM(multivariate_simulation$x, multivariate_simulation$s)

  expect_equal(estimated_theta,
               readRDS(test_path("fixtures", "two_component_3D_supervised_estimation.rds")))
})







