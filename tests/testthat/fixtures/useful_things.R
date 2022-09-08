##################################################################
##                       test simulations                       ##
##################################################################
set.seed(20)
# define parameters of a two component multivariate GMM
true_theta <- list(p=c(0.2, 0.8),
                   mu=matrix(c(20, 22, 22, 20), nrow = 2),
                   sigma=array(c(1, 0.2, 0.2, 1, 1, -0.2, -0.2, 1), dim=c(2, 2, 2)))
obervations <- simulate_multivariate_GMM (theta=true_theta, n=100)
saveRDS(obervations, testthat::test_path("fixtures", "two_component_multivariate_GMM.rds"))


true_theta <- list(p=c(0.2, 0.8),
                   mu=matrix(c(20, 20, 20, 40, 40, 40), nrow = 3),
                   sigma=array(rep(c(1, 0, 0.1, 0, 1, -0.1, 0.1, -0.1, 1), 2), dim=c(3, 3, 2)))

multivariate_simulation <- simulate_multivariate_GMM (theta=true_theta, n=2000)
estimated_theta <- estimate_supervised_multivariate_GMM(multivariate_simulation$x, multivariate_simulation$s)

saveRDS(estimated_theta, testthat::test_path("fixtures", "two_component_3D_supervised_estimation.rds"))
