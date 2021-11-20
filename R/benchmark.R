#' Launch the benchmark to compare statistical performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param sigma_values,mean_values,proportions,skewness_values the true parameters to be retrieved
#' @param prop_outliers the proportion of outliers added in the simulation
#' @param nobservations the number of observations drawn to generate the random sample
#' @param Nbootstrap the number of bootstrap simulations and repetitions to perform
#' @param epsilon,itmax respectively criterion threshold and maximal number of iterations to reach it
#' @param nstart,short_iter,short_eps hyper-parameters to control the initialisation step
#' @param prior_prob add minimal uncertainty on the cluster assignment returned by hierarchical clustering method
#' @param initialisation_algorithms among 6 methods, which algorithms to be chosen for the initialisation phase
#'
#' @return a list with the simulated distributions of the estimates, some summary scores per parameter and aggregated measures
#' as well as boxplot and Heatmap correlation representations of the estimates
#'
#' @importFrom magrittr "%>%"
#' @import ggplot2
#'
#' @export

benchmark_distribution_parameters <- function(mixture_functions,
                                              sigma_values, mean_values, proportions, skewness_values,
                                              prop_outliers = 0, nobservations = c(100, 1000, 10000),
                                              Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                              initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################

  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble()
  global_scores <- tibble::tibble() # summary scores
  parameters_plots <- list()
  correlation_scores <- list()
  correlation_scores_plots <- list() # visual representations of the summary scores


  for (prop_out in prop_outliers) {
    for (skew in skewness_values) {
      for (p in proportions) {
        for (sigma in sigma_values) {
          for (mu in mean_values) {
            #################################################################
            ##               simulation scenario description               ##
            #################################################################
            true_theta <- list(p = p, mu = mu, sigma = sigma, skew = skew) # true parameters of the distribution
            k <- length(p) # number of components
            bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
            balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega # compute OVL
            pairwise_ovl <- MixSim::overlap(Pi = true_theta$p, Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega
            entropy_value <- compute_shannon_entropy(p) # compute entropy


            for (n in nobservations) {
              filename <- paste0(
                n, "_observations_entropy", signif(entropy_value, digits = 2), "_skewness_", skew[1],
                "_OVL_", signif(balanced_ovl, digits = 2), "_prop_outliers_", prop_out
              )
              global_scores_temp <- tibble::tibble()
              distribution_parameters_temp <- tibble::tibble() # store temp summary scores
              for (t in 1:Nbootstrap) {
                simulated_distribution <- rnmix_skewed_with_outliers(n = n, theta = true_theta, prop_outliers = prop_out, interval = 2) # simulation
                for (init_algo in initialisation_algorithms) {
                  ##################################################################
                  ##             estimation of the initial estimates             ##
                  ##################################################################
                  start <- initialize_em(
                    x = simulated_distribution$x, k = k, nstart = nstart,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  )

                  # compute the estimates data with complete observations
                  observed_estimated_theta <- compute_parameters_complete_observations(simulated_distribution)

                  ##################################################################
                  ##           estimation of GMMs with the EM algorithm           ##
                  ##################################################################
                  for (index in 1:length(mixture_functions)) {
                    # retrieve function values
                    mixture_function <- mixture_functions[[index]]$name_fonction
                    package_name <- names(mixture_functions)[index]
                    if (init_algo != "hc" & identical(em_otrimle, mixture_function)) {
                      next
                    } # other methods than hc are not implemented with otrimle, we go out from the loop


                    success_estimation <- TRUE # checking that the algorithm is not trapped in boundary space or simply failed
                    tryCatch(
                      {
                        missing_estimated_theta_list <- do.call(mixture_function, c(
                          x = list(simulated_distribution$x), k = simulated_distribution$k,
                          epsilon = epsilon, itmax = itmax,
                          start = list(start), mixture_functions[[index]]$list_params
                        ))
                      },
                      error = function(e) {
                        success_estimation <- FALSE
                        warning(paste("Error is:", e))
                      }
                    )
                    if (!check_parameters_validity(missing_estimated_theta_list, k = k) | !success_estimation) {
                      warning(paste("Estimates trapped in the boundary space or failure of the package", package_name))
                      next
                    }
                    # normalize returned estimates
                    missing_estimated_theta_list <- missing_estimated_theta_list[c("p", "mu", "sigma")]
                    missing_estimated_theta_list$p <- missing_estimated_theta_list$p / sum(missing_estimated_theta_list$p)
                    missing_estimated_theta <- stats::setNames(unlist(missing_estimated_theta_list), bootstrap_colnames)
                    logLikelihood <- EMCluster::logL(as.matrix(simulated_distribution$x),
                      emobj = list(
                        pi = missing_estimated_theta_list$p, Mu = as.matrix(missing_estimated_theta_list$mu),
                        LTSigma = as.matrix(missing_estimated_theta_list$sigma^2)
                      )
                    )
                    # store distribution of the estimates
                    distribution_parameters_temp <- distribution_parameters_temp %>%
                      dplyr::bind_rows(tibble::tibble(
                        package = package_name,
                        initialisation_method = init_algo,
                        tibble::as_tibble_row(missing_estimated_theta),
                        OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2),
                        entropy = signif(entropy_value, digits = 2), logLikelihood = logLikelihood,
                        skew = skew[1], nobservations = n, prop_outliers = prop_out, N.bootstrap = t
                      ))

                    global_scores_temp <- global_scores_temp %>% dplyr::bind_rows(
                      tibble::tibble(logLikelihood = logLikelihood, initialisation_method = init_algo, package = package_name)
                    )
                  } # estimation per package
                } # initialization algorithm
              } # bootstrap repetitions

              #################################################################
              ##     summarize results obtained per scenario configuration   ##
              #################################################################
              distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_temp)

              # summary scores per parameter
              local_scores_temp <- distribution_parameters_temp %>%
                dplyr::group_by(initialisation_method, package) %>%
                dplyr::summarize(get_local_scores(dplyr::cur_data() %>%
                  dplyr::select(dplyr::all_of(bootstrap_colnames)), true_theta[c("p", "mu", "sigma")])) %>%
                tibble::add_column(
                  OVL = signif(balanced_ovl, digits = 2), entropy = signif(entropy_value, digits = 2),
                  OVL_pairwise = signif(pairwise_ovl, digits = 2), skew = skew[1], nobservations = n, prop_outliers = prop_out
                )
              local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)

              # aggregate scores obtained per parameter
              global_bias <- local_scores_temp %>%
                dplyr::filter(scores == "bias") %>%
                dplyr::rowwise() %>%
                dplyr::transmute(global_bias = sum(abs(dplyr::c_across(dplyr::all_of(bootstrap_colnames))))) # global bias
              global_mse <- local_scores_temp %>%
                dplyr::filter(scores == "mse") %>%
                dplyr::rowwise() %>%
                dplyr::transmute(global_mse = sum(dplyr::c_across(dplyr::all_of(bootstrap_colnames)))) # global mse
              missed_cases <- distribution_parameters_temp %>%
                dplyr::count(initialisation_method, package, name = "N.missed") %>%
                dplyr::mutate(N.missed = Nbootstrap - N.missed) # count number of times the package failed in the estimation with the given number of clusters

              global_scores_temp <- global_scores_temp %>%
                dplyr::group_by(initialisation_method, package) %>%
                dplyr::summarise(
                  logLikelihood = mean(logLikelihood),
                  entropy = signif(entropy_value, digits = 2),
                  OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2), skew = skew[1], nobservations = n, prop_outliers = prop_out
                ) %>%
                dplyr::inner_join(missed_cases, by = c("initialisation_method", "package")) %>%
                dplyr::inner_join(global_bias, by = c("initialisation_method", "package")) %>%
                dplyr::inner_join(global_mse, by = c("initialisation_method", "package"))
              global_scores <- global_scores %>% dplyr::bind_rows(global_scores_temp)

              # correlation matrices between the estimates returned by the packages
              distribution_parameters_temp_long <- distribution_parameters_temp %>%
                dplyr::select(c(
                  dplyr::all_of(bootstrap_colnames), package, initialisation_method,
                  OVL, skew, nobservations, prop_outliers, N.bootstrap
                )) %>%
                tidyr::pivot_longer(dplyr::all_of(bootstrap_colnames), names_to = "name_parameter", values_to = "value_parameter") %>%
                dplyr::mutate(
                  name_parameter = factor(name_parameter, levels = bootstrap_colnames),
                  package = factor(package, levels = unique(distribution_parameters_temp$package))
                )

              correlation_scores[[filename]] <- lapply(
                split(distribution_parameters_temp_long, distribution_parameters_temp_long$initialisation_method),
                function(x) {
                  cor_data <- tidyr::pivot_wider(x, names_from = c("package"), values_from = "value_parameter") %>%
                    dplyr::select(unique(x %>% dplyr::pull(package))) %>%
                    stats::cor(use = "complete.obs")
                  return(cor_data)
                }
              )


              ### boxplots of the estimated parameters
              parameters_plots[[filename]] <- plot_boxplots_parameters(distribution_parameters_temp_long, p, mu, sigma)

              ### Heatmap of the estimated parameters
              correlation_scores_plots[[filename]] <- Map(function(cor_matrix, init_method) {
                complex_heatmap <- ComplexHeatmap::Heatmap(cor_matrix,
                  name = "mat", heatmap_legend_param = list(title = ""),
                  cluster_rows = TRUE, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix),
                  row_title_gp = grid::gpar(fontsize = 4), column_names_rot = 45,
                  cluster_columns = TRUE, column_names_gp = grid::gpar(fontsize = 8), column_labels = colnames(cor_matrix),
                  width = unit(8, "cm"), height = unit(8, "cm"),
                  column_title = paste(init_method),
                  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
                )
                return(grid::grid.grabExpr(ComplexHeatmap::draw(complex_heatmap)))
              }, correlation_scores[[filename]], names(correlation_scores[[filename]]))
            } # number of observations
          }
        }
      } ##### theta configuration
    }
  }

  return(list(
    "distributions" = distribution_parameters,
    "local_scores" = local_scores, "global_scores" = global_scores,
    "plots" = parameters_plots, "correlation_scores_plots" = correlation_scores_plots
  ))
}













#' Launch the benchmark to compare computational performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param sigma_values,mean_values,proportions,skewness_values the true parameters to be retrieved
#' @param prop_outliers the proportion of outliers added in the simulation
#' @param nobservations the number of observations drawn to generate the random sample
#' @param Nbootstrap the number of bootstrap simulations and repetitions to perform
#' @param epsilon,itmax respectively criterion threshold and maximal number of iterations to reach it
#' @param nstart,short_iter,short_eps hyper-parameters to control the initialisation step
#' @param prior_prob add minimal uncertainty on the cluster assignment returned by hierarchical clustering method
#' @param initialisation_algorithms among 6 methods, which algorithms to be chosen for the initialisation phase
#'
#' @return a list with the running time of the initialisation and the EM estimation itself, as well as corresponding time curve representations
#'
#' @export

# function used to compare time computations between packages
compute_microbenchmark <- function(mixture_functions,
                                   sigma_values, mean_values, proportions,
                                   skewness_values, prop_outliers = 0,
                                   nobservations = c(100, 1000, 10000),
                                   Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                   nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                   initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {

  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  time_data <- tibble::tibble()
  init_time_data <- tibble::tibble()
  for (prop_out in prop_outliers) {
    for (skew in skewness_values) {
      for (p in proportions) {
        for (sigma in sigma_values) {
          for (mu in mean_values) {
            #################################################################
            ##               simulation scenario description               ##
            #################################################################
            true_theta <- list(p = p, mu = mu, sigma = sigma, skew = skew) # true parameters of the distribution
            k <- length(p) # number of components
            bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
            balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega # compute OVL
            pairwise_ovl <- MixSim::overlap(Pi = true_theta$p, Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega
            entropy_value <- compute_shannon_entropy(p) # compute entropy

            for (n in nobservations) {
              filename <- paste0(
                n, "_observations_entropy", signif(entropy_value, digits = 2), "_skewness_", skew[1],
                "_OVL_", signif(balanced_ovl, digits = 2), "_prop_outliers_", prop_out
              )
              simulated_distribution <- rnmix_skewed_with_outliers(n = n, theta = true_theta, prop_outliers = prop_out, interval = 2) # simulation of the experience
              mbm_temp <- tibble::tibble()
              init_temp <- tibble::tibble()

              for (t in 1:Nbootstrap) {
                for (init_algo in initialisation_algorithms) {
                  ##################################################################
                  ##        estimation of time taken for the initialisation       ##
                  ##################################################################
                  start <- initialize_em(
                    x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  )
                  # and determine time taken for the initialisation step
                  init_temp <- init_temp %>% dplyr::bind_rows(
                    microbenchmark::microbenchmark(initialisation_method = initialize_em(
                      x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                      short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                    ), times = 1) %>%
                      tibble::as_tibble() %>% dplyr::rename(initialisation_method = expr) %>%
                      dplyr::mutate(time = microbenchmark:::convert_to_unit(time, "s"), initialisation_method = init_algo, N.bootstrap = t)
                  )

                  # perform micro-benchmark
                  for (index in 1:length(mixture_functions)) {
                    ##################################################################
                    ##        estimation of time taken by the EM algorithm          ##
                    ##################################################################

                    mixture_function <- mixture_functions[[index]]$name_fonction
                    package_name <- names(mixture_functions)[index]
                    if (init_algo != "hc" & identical(em_otrimle, mixture_function)) {
                      next
                    } # skip otrimle, when no initialization algorithm is relevant
                    mbm_temp <- mbm_temp %>% dplyr::bind_rows(
                      microbenchmark::microbenchmark(package = do.call(mixture_function, c(
                        x = list(simulated_distribution$x), k = simulated_distribution$k,
                        epsilon = epsilon, itmax = itmax, start = list(start), mixture_functions[[index]]$list_params
                      )), times = 1) %>%
                        tibble::as_tibble() %>% dplyr::rename(package = expr) %>%
                        dplyr::mutate(
                          time = microbenchmark:::convert_to_unit(time, "s"), initialisation_method = init_algo,
                          N.bootstrap = t, package = package_name
                        )
                    )
                  } # mixture package
                } # initialization algorithm
              } # bootstrap loop
              init_temp <- init_temp %>% dplyr::mutate(
                entropy = signif(entropy_value, digits = 2),
                OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2),
                skew = skew[1], nobservations = n, prop_outliers = prop_out
              )
              init_time_data <- init_time_data %>% dplyr::bind_rows(init_temp) # store time computation taken by the initalisation step

              mbm_temp <- mbm_temp %>% dplyr::mutate(
                entropy = signif(entropy_value, digits = 2),
                OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2),
                skew = skew[1], nobservations = n, prop_outliers = prop_out
              )
              time_data <- time_data %>% dplyr::bind_rows(mbm_temp) # store time computation taken by the EM algorithm
            } # number of observations
          } # theta configuration
        }
      }
    }
  }
  ##################################################################
  ##          visualize time computations' distributions          ##
  ##################################################################

  init_time_plots <- plot_initialisation_time_computations(init_time_data) # plot quantiles of time taken by the initialization step
  time_plots <- plot_time_computations(time_data) # plot quantiles of time taken by the EM algorithm whether the package
  return(list(
    time_data = time_data, time_plots = time_plots,
    init_time_data = init_time_data, init_time_plots = init_time_plots
  ))
}
