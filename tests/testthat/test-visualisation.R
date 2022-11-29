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

test_that("multivariate visualisation", {
  multivariate_parameters_distribution <- readRDS(testthat::test_path("fixtures", "multivariate_test_distribution.rds"))$distributions
  multivariate_config <- readRDS(testthat::test_path("fixtures", "multivariate_test_distribution.rds"))$config
  true_theta <- multivariate_config %>%
    dplyr::filter(ID == 1) %>%
    pull(true_parameters) %>%
    magrittr::extract2(1)


  multivariate_boxplot <- plot_boxplots_parameters(multivariate_parameters_distribution %>% filter(ID == 1),
    num_col = 2, true_theta = true_theta
  )
  multivariate_hellinger <- plot_Hellinger(multivariate_parameters_distribution %>% filter(ID == 1), true_theta = true_theta)
  multivariate_ellipses <- plot_ellipses_bivariate(multivariate_parameters_distribution %>% filter(ID == 1), true_theta)
  multivariate_density_distribution <- plot_bivariate_normal_density_distribution(true_theta)
  correlation_scores <- plot_correlation_Heatmap(multivariate_parameters_distribution %>% filter(ID == 1))

  multivariate_time_data <- readRDS(testthat::test_path("fixtures", "multivariate_test_time_computations.rds"))$time_data
  multivariate_time_computation <- plot_time_computations(multivariate_time_data %>% filter(ID == 1))
  vdiffr::expect_doppelganger("multivariate time computation", multivariate_time_computation)
})


# plot the 14 parametrisations available
test_that("multivariate parametrisations", {

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

  vdiffr::expect_doppelganger("parametrisation of gaussian mixtures", plot_multi_parametrisation(list_models$"EII"))
  # dir.create ("./figs/multivariate_parametrisations", showWarnings = F, recursive = T)
  # list_figures <- purrr::imap(list_models, function(.x, .y) {
  #   filename <- file.path(".", "figs", "multivariate_parametrisations", paste0(.y, ".png"))
  #   plot_param <- plot_multi_parametrisation(.x)
  #   ggsave(filename, plot_param, dpi = 150)
  #   return(plot_param)
  # })
})


# #### univariate distributions
# univariate_distribution_parameters <- readRDS("../mixture_models/results/univariate/univariate_distributions.rds") %>%
#   dplyr::group_by(ID, package, initialisation_method) %>%
#   dplyr::mutate(N.bootstrap=dplyr::row_number()) %>% dplyr::ungroup()
# saveRDS(univariate_distribution_parameters, "../mixture_models/results/univariate/univariate_distributions.rds")
# univariate_configuration <- readRDS("../mixture_models/results/univariate/univariate_configuration_scenario.rds")
# ## test representations
# true_theta <- univariate_configuration %>% filter(ID=="U9") %>%
#   pull(true_parameters) %>% magrittr::extract2(1)
# RGMMBench::plot_boxplots_parameters(univariate_distribution_parameters %>% filter(ID=="U9"), true_theta = true_theta)
# RGMMBench::plot_correlation_Heatmap(univariate_distribution_parameters %>% filter(ID=="U9"))
# RGMMBench::plot_Hellinger(univariate_distribution_parameters %>% filter(ID=="U9"), true_theta = true_theta)
# RGMMBench::plot_univariate_normal_density_distribution(true_theta = true_theta, nobservations = 2000)
#
#


# multivariate_distribution_parameters <- readRDS("../mixture_models/results/multivariate/bivariate_distributions.rds") %>%
#   dplyr::group_by(ID, package, initialisation_method) %>%
#   dplyr::mutate(N.bootstrap=dplyr::row_number()) %>% dplyr::ungroup()

# multivariate_distribution_parameters <- multivariate_distribution_parameters %>%
#   dplyr::relocate(sd_var1_var1_comp2, .after=sd_var1_var1_comp1) %>%
#   dplyr::relocate(sd_var2_var1_comp2, .after=sd_var2_var1_comp1) %>%
#   dplyr::relocate(sd_var2_var2_comp2, .after=sd_var2_var2_comp1)
#
# saveRDS(multivariate_distribution_parameters, "../mixture_models/results/multivariate/bivariate_distributions.rds")

# multivariate_configuration <- readRDS("../mixture_models/results/multivariate/bivariate_configuration_scenario.rds")
#
# liste_multivariate_plots <- purrr::map(multivariate_configuration %>% pull(ID), function(ID_index) {
#   true_theta <- multivariate_configuration %>% filter(ID==ID_index) %>%
#     pull(true_parameters) %>% magrittr::extract2(1)
#
#   boxplot <- RGMMBench::plot_boxplots_parameters(multivariate_distribution_parameters %>% filter(ID==ID_index), true_theta = true_theta)
#   hellinger <- RGMMBench::plot_Hellinger(multivariate_distribution_parameters %>% filter(ID==ID_index), true_theta = true_theta)
#   density <- RGMMBench::plot_bivariate_normal_density_distribution(true_theta = true_theta, nobservations = 2000)
#   intervals <- RGMMBench::plot_ellipses_bivariate(multivariate_distribution_parameters %>% filter(ID==ID_index), true_theta = true_theta, npoints = 2000)
#   general_plot <- gridExtra::arrangeGrob(grobs=list(density, intervals, hellinger, boxplot), nrow = 2, ncol = 2, top = paste0("ID is ", ID_index))
#   return(general_plot)
# })
#
# ggsave("./images/general_plots_summary.pdf",
#        gridExtra::marrangeGrob(grobs = liste_multivariate_plots, ncol = 1, nrow = 1, top=NULL), width = 20, height = 30 )





# distributions_plots <- purrr::imap(split(distribution_parameters, distribution_parameters$true_parameters_factor),
#                                    function(dist_per_fact, title) {
#                                 return(list("boxplot_parameter"=plot_boxplots_parameters_bivariate(dist_per_fact, num_col = 4) + ggtitle(title),
#                                 "ellipse"=plot_ellipses_bivariate(dist_per_fact) + ggtitle(title),
#                                 "hellinger"=plot_Hellinger(dist_per_fact) + ggtitle(title))) })
#
# ellipse_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("ellipse"), ncol = 1, nrow = 1, top=NULL)
# ggsave("./images/ellipse_plots.pdf", ellipse_plots, width = 15, height = 8,dpi = 300)
#
# boxplot_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("boxplot_parameter"), ncol = 1, nrow = 1, top=NULL)
# ggsave("./images/boxplot_parameters.pdf", boxplot_plots, width = 15, height = 15,dpi = 300)
#
# hellinger_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("hellinger"), ncol = 1, nrow = 1, top=NULL)
# ggsave("./images/hellinger_plot.pdf", hellinger_plots, width = 15, height = 15,dpi = 300)
#
# hellinger_plot <- plot_Hellinger(splitted_distribution[[11]])
# ggsave("./images/hellinger_plot_specific.pdf", hellinger_plot, width = 15, height = 15,dpi = 300)
