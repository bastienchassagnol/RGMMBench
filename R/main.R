# load auxiliary function files
source("mixture.R")
source("benchmark.R")
source("visualisation.R")

# load useful libraries and packages
library(ggplot2)
import::from(magrittr, "%>%", .into = "operators")
import::from(rebmix, .except = c("AIC", "BIC", "split"))
library(mclust)
library(Rmixmod)


## test functions

# simulated_data <- rnmix_skewed_with_outliers (200, theta =list(p=c(0.5, 0.5), mu=c(0, 4),sigma=c(0.3, 0.3), skew=c(0, 0)))
# start_theta <- initialize_em(simulated_data$x, k = 2)
# estimated_theta <- emnmix(simulated_data$x, k=2)

relevant_mixture_functions <- list(
  "otrimle" = list(name_fonction = em_otrimle, list_params = list()),
  "mixsmsn" = list(name_fonction = em_mixsmsn, list_params = list()),
  "em R" = list(name_fonction = emnmix, list_params = list()),
  "Rmixmod" = list(name_fonction = em_Rmixmod, list_params = list()),
  "mixtools" = list(name_fonction = em_mixtools, list_params = list()),
  "bgmm" = list(name_fonction = em_bgmm, list_params = list()),
  "mclust" = list(name_fonction = em_mclust, list_params = list(prior = NULL)),
  "EMCluster" = list(name_fonction = em_EMCluster, list_params = list()),
  "GMKMcharlie" = list(name_fonction = em_GMKMcharlie, list_params = list()),
  "flexmix" = list(name_fonction = em_flexmix, list_params = list()),
  "DCEM" = list(name_fonction = em_DCEM, list_params = list())
)

##################################################################
##       Compare statistical performances of the packages       ##
##################################################################

# Main simulation with four components
four_components_statistical_performances <- benchmark_distribution_parameters(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("null OVL" = rep(0.3, 4), "high OVL" = rep(2, 4), "average OVL" = rep(1, 4)),
  mean_values = list(c(0, 4, 8, 12)),
  proportions = list("balanced" = rep(1 / 4, 4), "small unbalance" = c(0.2, 0.4, 0.2, 0.2), "highly unbalanced" = c(0.1, 0.7, 0.1, 0.1)),
  skewness_values = list("null skewness" = c(0, 0, 0, 0)), prop_outliers = c(0),
  Nbootstrap = 200, nobservations = c(2000)
)

# Main simulation with two components
two_components_statistical_performances <- benchmark_distribution_parameters(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("null OVL" = rep(0.3, 2), "average OVL" = rep(1, 2), "high OVL" = rep(2, 2)),
  mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.5, 0.5), "small unbalanced" = c(0.65, 1), "high unbalanced" = c(0.9, 0.1)),
  skewness_values = list("null skewness" = c(0, 0)), prop_outliers = c(0),
  Nbootstrap = 200, nobservations = c(2000)
)

# Test influence of the number of observations
rmse_distribution <- benchmark_distribution_parameters(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("high OVL" = c(2, 2)), mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.9, 0.1)), skewness_values = list("null skewness" = c(0, 0)),
  prop_outliers = c(0), Nbootstrap = 200, nobservations = c(100, 1000, 10000)
)

# Test influence of adding outliers in the distribution
outliers_distribution <- benchmark_distribution_parameters(
  mixture_functions = relevant_mixture_functions, initialisation_algorithms = c("hc"),
  sigma_values = list("null OVL" = c(0.3, 0.3)),
  mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.5, 0.5)),
  skewness_values = list("null skewness" = c(0, 0)), prop_outliers = c(0.02, 0.04),
  Nbootstrap = 200, nobservations = c(2000)
)

# Test influence of adding a scale factor on the quality of the estimates
skewed_distributions <- benchmark_distribution_parameters(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("null OVL" = c(0.3, 0.3)),
  mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.5, 0.5)),
  skewness_values = list("null skewness" = c(0, 0), "high skewness" = c(20, 20)),
  prop_outliers = c(0), Nbootstrap = 200, nobservations = c(2000)
)

##################################################################
##      Compare computational performances of the packages      ##
##################################################################
# Main simulation with four components
four_components_computational_performances <- compute_microbenchmark(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("null OVL" = rep(0.3, 4), "high OVL" = rep(2, 4), "average OVL" = rep(1, 4)),
  mean_values = list(c(0, 4, 8, 12)),
  proportions = list("balanced" = rep(1 / 4, 4), "small unbalance" = c(0.2, 0.4, 0.2, 0.2), "highly unbalanced" = c(0.1, 0.7, 0.1, 0.1)),
  skewness_values = list("null skewness" = c(0, 0, 0, 0)), prop_outliers = c(0),
  Nbootstrap = 5, nobservations = c(100, 200, 500, 1000, 2000, 10000)
)

# Main simulation with two components
two_components_computational_performances <- compute_microbenchmark(
  mixture_functions = relevant_mixture_functions,
  sigma_values = list("null OVL" = rep(0.3, 2), "average OVL" = rep(1, 2), "high OVL" = rep(2, 2)),
  mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.5, 0.5), "small unbalanced" = c(0.65, 1), "high unbalanced" = c(0.9, 0.1)),
  skewness_values = list("null skewness" = c(0, 0)), prop_outliers = c(0),
  Nbootstrap = 200, nobservations = c(100, 200, 500, 1000, 2000, 10000)
)


# Test influence of adding outliers in the distribution
outliers_computational_results <- compute_microbenchmark(
  mixture_functions = relevant_mixture_functions, initialisation_algorithms = "hc",
  sigma_values = list("small OVL" = c(0.3, 0.3)),
  mean_values = list(c(0, 4)), prop_outliers = c(0.02, 0.04),
  proportions = list("balanced" = c(0.5, 0.5)),
  skewness_values = list("null skewness" = c(0, 0)), Nbootstrap = 50, nobservations = c(100, 200, 500, 1000, 2000, 5000, 10000)
)

# Test influence of adding a scale factor on the quality of the estimates
skewed_computational_results <- compute_microbenchmark(
  mixture_functions = relevant_mixture_functions, initialisation_algorithms = "kmeans",
  sigma_values = list("small OVL" = c(0.3, 0.3)),
  mean_values = list(c(0, 4)),
  proportions = list("balanced" = c(0.5, 0.5)),
  skewness_values = list("high skewness" = c(20, 20)), Nbootstrap = 50, nobservations = c(100, 200, 500, 1000, 2000, 5000, 10000)
)


#################################################################
##  Save results (example with the four components simulation  ##
#################################################################

# save summary scores and distributions of the bootstrap simulations
openxlsx::write.xlsx(four_components_computational_performances$local_scores, file = "tables/four_components_local_scores.xlsx", asTable = T)
openxlsx::write.xlsx(four_components_computational_performances$global_scores, file = "tables/four_components_global_scores.xlsx", asTable = T)
openxlsx::write.xlsx(four_components_computational_performances$distributions, file = "tables/four_components_distributions.xlsx", asTable = T)

# save boxplots associated to the distribution of the estimates
boxplot_parameters <- gridExtra::marrangeGrob(four_components_computational_performances$plots, nrow = 1, ncol = 1, top = "")
ggsave("images/four_components_parameters_boxplots.pdf", boxplot_parameters,
  width = 15, height = 14, dpi = 600
)
# save correlation Heatmaps associated to the distribution of the estimates
correlation_scores_plots <- lapply(four_components_computational_performances$correlation_scores_plots, function(x) gridExtra::arrangeGrob(grobs = x, nrow = 3, ncol = 2, top = ""))
correlation_scores_plots <- gridExtra::marrangeGrob(correlation_scores_plots, nrow = 1, ncol = 1, top = "")
ggsave("images/four_components_parameters_correlation_heatmap.pdf", correlation_scores_plots,
  width = 15, height = 14, dpi = 600
)

# save quantiles of the computational time taken for the estimation of GMMs
openxlsx::write.xlsx(four_components_computational_performances$time_data, file = "tables/four_components_EM_time_computations.xlsx", asTable = T)
openxlsx::write.xlsx(four_components_computational_performances$init_time_data, file = "tables/four_components_initialization_time_computations.xlsx", asTable = T)


# save time estimates associated to the initialization step and the EM estimation of the GMMs
time_plots <- gridExtra::marrangeGrob(four_components_computational_performances$time_plots, nrow = 1, ncol = 1, top = "")
ggsave("./images/four_components_time_computations.pdf", time_plots, width = 7, height = 8, dpi = 600)
ggsave("./images/four_components_initialization_time_computations.pdf",
  four_components_computational_performances$init_time_plots,
  width = 7, height = 7, dpi = 600
)
