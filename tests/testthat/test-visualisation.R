library(dplyr)
test_that("univariate visualisation", {
  univariate_parameters_distribution <- readRDS(testthat::test_path("fixtures", "univariate_test_distribution.rds"))$distributions
  univariate_config <- readRDS(testthat::test_path("fixtures", "univariate_test_distribution.rds"))$config
  true_theta <- univariate_config %>%
    dplyr::filter(ID == 1) %>%
    pull(true_parameters) %>%
    magrittr::extract2(1)

  univariate_boxplot <- plot_boxplots_parameters(univariate_parameters_distribution %>% filter(ID == 1),
    num_col = 2, true_theta = true_theta
  )
  univariate_hellinger <- plot_Hellinger(univariate_parameters_distribution %>% filter(ID == 1), true_theta = true_theta)
  univariate_density_distribution <- plot_univariate_normal_density_distribution(true_theta)

  univariate_time_data <- readRDS(testthat::test_path("fixtures", "univariate_test_time_computations.rds"))$time_data
  univariate_time_computation <- plot_time_computations(univariate_time_data %>% dplyr::filter(ID == 1))
  vdiffr::expect_doppelganger("univariate time computation", univariate_time_computation)


  # initialisation_time_data <- readRDS("../mixture_models/results/univariate/univariate_initialisation_time_computation.rds")
  # inialisation_plot <- plot_initialisation_time_computations(initialisation_time_data, size=4, x = 2.6)
  # ggsave("./figs/univariate_initialisation_time_computations.png", inialisation_plot, dpi = 300)
})

test_that("bivariate visualisation", {
  bivariate_parameters_distribution <- readRDS(testthat::test_path("fixtures", "bivariate_test_distribution.rds"))$distributions
  bivariate_config <- readRDS(testthat::test_path("fixtures", "bivariate_test_distribution.rds"))$config
  true_theta <- bivariate_config %>%
    dplyr::filter(ID == 1) %>%
    pull(true_parameters) %>%
    magrittr::extract2(1)


  bivariate_boxplot <- plot_boxplots_parameters(bivariate_parameters_distribution %>% filter(ID == 1),
    num_col = 2, true_theta = true_theta
  )
  bivariate_hellinger <- plot_Hellinger(bivariate_parameters_distribution %>% filter(ID == 1), true_theta = true_theta)
  bivariate_ellipses <- plot_ellipses_bivariate(bivariate_parameters_distribution %>% filter(ID == 1), true_theta)
  bivariate_density_distribution <- plot_bivariate_normal_density_distribution(true_theta)
  correlation_scores <- plot_correlation_Heatmap(bivariate_parameters_distribution %>% filter(ID == 1))

  bivariate_time_data <- readRDS(testthat::test_path("fixtures", "bivariate_test_time_computations.rds"))$time_data
  bivariate_time_computation <- plot_time_computations(bivariate_time_data %>% filter(ID == 1))
  vdiffr::expect_doppelganger("bivariate time computation", bivariate_time_computation)
})


# plot the 14 parametrisations available
test_that("bivariate parametrisations", {

  vec_pos_correlated <- eigen(matrix(c(2, 0.8,0.8, 1), nrow = 2))$vectors
  vec_neg_correlated <- eigen(matrix(c(2, -0.8, -0.8, 1), nrow = 2))$vectors



  list_models <- list ("EII"= list(p=c(1/3, 1/3, 1/3),
                                  mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                  sigma=array(rep(c(diag(c(1,1))), 3), dim = c(2, 2, 3))),
                      "VII" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(diag(c(1,1)), diag(c(0.5,0.5)), diag(c(2,2))), dim = c(2, 2, 3))),


                      "EEI" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(rep(c(diag(c(2, 1))), 3), dim = c(2, 2, 3))),
                      "VEI" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(1 * diag(c(sqrt(2)/2, sqrt(2))), 0.5 * diag(c(sqrt(2)/2, sqrt(2))), 2 * diag(c(sqrt(2)/2, sqrt(2)))), dim = c(2, 2, 3))),
                      "EVI" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(diag(c(2,1)), diag(c(1, 2)), diag(c(0.5, 4))), dim = c(2, 2, 3))),
                      "VVI" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(diag(c(2.5,1)), diag(c(1, 1)), diag(c(0.5, 1))), dim = c(2, 2, 3))),


                      "EEE" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(rep(c(1, 0.8,0.8, 2), 3), dim = c(2, 2, 3))),
                      "EVE" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(vec_pos_correlated %*% diag(c(2/3, 3)) %*% t(vec_pos_correlated),
                                                 vec_pos_correlated %*%  diag(c(1, 2)) %*% t(vec_pos_correlated),
                                                 vec_pos_correlated %*% diag(c(0.5, 4)) %*% t(vec_pos_correlated)), dim = c(2, 2, 3))),
                      "VEE" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(1 * vec_pos_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_pos_correlated),
                                                 0.5 * vec_pos_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_pos_correlated),
                                                 2 * vec_pos_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_pos_correlated)), dim = c(2, 2, 3))),
                      "VVE" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(1 * vec_pos_correlated %*% diag(c(1, 2))%*% t(vec_pos_correlated),
                                                 0.5 * vec_pos_correlated %*% diag(c(2/3, 3)) %*% t(vec_pos_correlated),
                                                 2 * vec_pos_correlated %*% diag(c(0.5, 4)) %*% t(vec_pos_correlated)), dim = c(2, 2, 3))),
                      "EEV" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(vec_pos_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_pos_correlated),
                                                 vec_neg_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_neg_correlated),
                                                 diag(c(sqrt(2)/2, sqrt(2)))), dim = c(2, 2, 3))),
                      "VEV" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(2 * vec_pos_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_pos_correlated),
                                                 0.5 * vec_neg_correlated %*% diag(c(sqrt(2)/2, sqrt(2))) %*% t(vec_neg_correlated),
                                                 diag(c(sqrt(2)/2, sqrt(2)))), dim = c(2, 2, 3))),
                      "EVV" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(vec_pos_correlated %*% diag(c(2,1))%*% t(vec_pos_correlated),
                                                 vec_neg_correlated %*% diag(c(2, 1)) %*% t(vec_neg_correlated),
                                                 diag(c(0.5, 4))), dim = c(2, 2, 3))),
                      "VVV" = list(p=c(1/3, 1/3, 1/3),
                                   mu=matrix(c(0,4, 4, 4, 2, 0), nrow = 2),
                                   sigma=array(c(2, 0.8,0.8, 1, 1, -0.8, -0.8, 2, 1, 0, 0, 1), dim = c(2, 2, 3))))


  # dir.create ("./figs/bivariate_parametrisations", showWarnings = F, recursive = T)
  # list_figures <- purrr::imap(list_models, function(.x, .y) {
  #   filename <- file.path(".", "figs", "bivariate_parametrisations", paste0(.y, ".png"))
  #   plot_param <- plot_multi_parametrisation(.x)
  #   ggsave(filename, plot_param, dpi = 150)
  #   return(plot_param)
  # })
})



test_that("HD visualisation", {
  library(dplyr)
  HD_parameters_distribution <- readRDS(testthat::test_path("fixtures", "HD_high_OVL_balanced_circular.rds"))$distributions
  HD_config <- readRDS(testthat::test_path("fixtures", "HD_high_OVL_balanced_circular.rds"))$config
  true_theta <- HD_config %>%
    dplyr::filter(ID == 1) %>%
    pull(true_parameters) %>%
    magrittr::extract2(1)


  HD_boxplot <- plot_boxplots_parameters(HD_parameters_distribution %>% filter(ID == 1),
                                         num_col = 2, true_theta = true_theta,
                                         match_symbol = "^p[[:digit:]]+$|mu_var1_|sd_var1_var1|sd_var2_var1|sd_var2_var2")


  HD_hellinger <- plot_Hellinger(HD_parameters_distribution %>% filter(ID == 1), true_theta = true_theta)

  correlation_scores <- plot_correlation_Heatmap(HD_parameters_distribution %>% filter(ID == 1))


  HD_density_distribution <- plot_HD_density_distribution(true_theta, nobservations = 200)
  ggsave(testthat::test_path("results", "HD_second_representation_test.pdf"), HD_density_distribution,
         width = 20, height = 16,dpi = 600)



  # HD_parallel <- MixSim::pdplot(Pi=true_theta$p, Mu=t(true_theta$mu), S=true_theta$sigma,
  #                               MaxInt = 0.8, marg = c(2,1,1,1))

  # HD_time_data <- readRDS(testthat::test_path("fixtures", "HD_test_time_computations.rds"))$time_data
  # HD_time_computation <- plot_time_computations(HD_time_data %>% filter(ID == 1))
  # vdiffr::expect_doppelganger("HD time computation", HD_time_computation)
})





