library(dplyr)
test_that("univariate visualisation", {
  univariate_parameters_distribution <- readRDS(testthat::test_path("fixtures", "univariate_test_distribution.rds"))$distributions
  univariate_config <- readRDS(testthat::test_path("fixtures", "univariate_test_distribution.rds"))$config
  true_theta <- univariate_config %>% dplyr::filter(ID==1) %>% pull(true_parameters) %>% magrittr::extract2(1)

  univariate_boxplot <- plot_boxplots_parameters(univariate_parameters_distribution %>% filter(ID==1),
                                                 num_col = 2, true_theta = true_theta)
  univariate_hellinger <- plot_Hellinger(univariate_parameters_distribution %>% filter(ID==1), true_theta = true_theta)
  univariate_density_distribution <- plot_univariate_normal_density_distribution(true_theta)

  univariate_time_data <- readRDS(testthat::test_path("fixtures", "univariate_test_time_computations.rds"))$time_data
  univariate_time_computation <- plot_time_computations(univariate_time_data %>% dplyr::filter(ID==1))
  vdiffr::expect_doppelganger("univariate time computation", univariate_time_computation)
})

test_that("multivariate visualisation", {
  multivariate_parameters_distribution <- readRDS(testthat::test_path("fixtures", "multivariate_test_distribution.rds"))$distributions
  multivariate_config <- readRDS(testthat::test_path("fixtures", "multivariate_test_distribution.rds"))$config
  true_theta <- multivariate_config %>% dplyr::filter(ID==1) %>% pull(true_parameters) %>% magrittr::extract2(1)


  multivariate_boxplot <- plot_boxplots_parameters(multivariate_parameters_distribution %>% filter(ID==1),
                                                 num_col = 2, true_theta = true_theta)
  multivariate_hellinger <- plot_Hellinger(multivariate_parameters_distribution %>% filter(ID==1), true_theta = true_theta)
  multivariate_ellipses <- plot_ellipses_bivariate (multivariate_parameters_distribution %>% filter(ID==1), true_theta)
  multivariate_density_distribution <- plot_bivariate_normal_density_distribution(true_theta)
  correlation_scores <- plot_correlation_Heatmap(multivariate_parameters_distribution %>% filter(ID==1))

  multivariate_time_data <- readRDS(testthat::test_path("fixtures", "multivariate_test_time_computations.rds"))$time_data
  multivariate_time_computation <- plot_time_computations(multivariate_time_data %>% filter(ID==1))
  vdiffr::expect_doppelganger("multivariate time computation", multivariate_time_computation)
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











