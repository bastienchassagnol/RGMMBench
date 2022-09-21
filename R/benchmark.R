#' Launch the benchmark to compare statistical performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param sigma_values,mean_values,proportions the true parameters to be retrieved
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

benchmark_univariate_GMM_estimation <- function(mixture_functions,
                                              sigma_values, mean_values, proportions,
                                              prop_outliers = 0, nobservations = c(100, 1000, 10000),
                                              Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                              initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################

  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble() # store bias and mse for each parameter


  for (prop_out in prop_outliers) {
    for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) {
          #################################################################
          ##               simulation scenario description               ##
          #################################################################
          true_theta <- list(p = p, mu = mu, sigma = sigma) # true parameters of the distribution
          k <- length(p) # number of components
          bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
          balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega %>%
            signif(digits = 2)# compute OVL
          pairwise_ovl <- MixSim::overlap(Pi = true_theta$p, Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega %>%
            signif(digits = 2)
          entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy


          for (n in nobservations) {
            filename <- paste0(n, "_observations_entropy", entropy_value, "_OVL_", balanced_ovl, "_prop_outliers_", prop_out)

            distribution_parameters_temp <- tibble::tibble() # store temp summary scores
            for (t in 1:Nbootstrap) {
              simulated_distribution <- simulate_univariate_GMM(n = n, theta = true_theta, prop_outliers = prop_out, interval = 2) # simulation
              for (init_algo in initialisation_algorithms) {
                ##################################################################
                ##             estimation of the initial estimates             ##
                ##################################################################
                initial_estimates <- initialize_em_univariate(
                  x = simulated_distribution$x, k = k, nstart = nstart,
                  short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                )
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
                        start = list(initial_estimates), mixture_functions[[index]]$list_params
                      ))
                    },
                    error = function(e) {
                      success_estimation <- FALSE
                      warning(paste("Error is:", e))
                    }
                  )
                  if (!check_parameters_validity_univariate(missing_estimated_theta_list, k = k) | !success_estimation) {
                    warning(paste("Estimates trapped in the boundary space or failure of the package", package_name))
                    next
                  }
                  # normalize returned estimates
                  missing_estimated_theta_list <- missing_estimated_theta_list[c("p", "mu", "sigma")]
                  missing_estimated_theta_list$p <- missing_estimated_theta_list$p / sum(missing_estimated_theta_list$p)
                  missing_estimated_theta <- stats::setNames(unlist(missing_estimated_theta_list), bootstrap_colnames)

                  # store distribution of the estimates
                  distribution_parameters_temp <- distribution_parameters_temp %>%
                    dplyr::bind_rows(tibble::tibble(
                      package = package_name,
                      initialisation_method = init_algo,
                      tibble::as_tibble_row(missing_estimated_theta), OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                      entropy = entropy_value, nobservations = n, prop_outliers = prop_out, N.bootstrap = t))

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
                OVL = balanced_ovl, entropy = entropy_value,
                OVL_pairwise = pairwise_ovl,  nobservations = n, prop_outliers = prop_out
              )
            local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)
          } # number of observations
        }
      }
    } ##### theta configuration

  }
  return(list("distributions" = distribution_parameters, "local_scores" = local_scores))
}



#' @rdname benchmark_univariate_GMM_estimation
#' @export

benchmark_multivariate_GMM_estimation <- function(mixture_functions,
                                              mean_values, proportions, sigma_values,
                                              nobservations = c(2000), Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                              initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################

  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble() # store bias and mse for each parameter

  for (p in proportions) {
    for (mu in mean_values) {
      for (sigma in sigma_values) {
        #################################################################
        ##               simulation scenario description               ##
        #################################################################
        true_theta <- list(p = p, mu = mu, sigma = sigma)
        formatted_true_theta <- true_theta %>% format_theta_output()
        k <- length(p); bootstrap_colnames <- names(formatted_true_theta)
        message(paste("We aim at learning the following estimates: \n",
                      paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse=" // ")))
        # bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
        balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)# compute OVL
        pairwise_ovl <- MixSim::overlap(Pi = p, Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy


        for (n in nobservations) {
          filename <- paste0(n, "_observations_entropy", entropy_value, "_OVL_", balanced_ovl)

          distribution_parameters_temp <- tibble::tibble() # store temp summary scores
          for (t in 1:Nbootstrap) {
            simulated_distribution <- simulate_multivariate_GMM(theta = true_theta, n = n) # simulation
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##             estimation of the initial estimates             ##
              ##################################################################
              good_initialisation <- TRUE
              tryCatch({
              initial_estimates <- initialize_em_multivariate(
                x = simulated_distribution$x, k = k, nstart = nstart,
                short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
              )}, error = function(e) {
                  good_initialisation <<- FALSE # use this operator, <<-, as in the other case, its value is not updated
                  warning(paste("Error in the initialisation step occurs:", e,
                                                              "with init algorithm: ", init_algo))
                }
              )
              if(!good_initialisation) {
                dir.create("./errors/initialisation_failures", showWarnings = F, recursive = T)
                saveRDS(object = list(simulated_distribution=simulated_distribution, algo=init_algo),
                        file = glue::glue("./errors/initialisation_failures/init_algo_{init_algo}_step_{t}_{n}_observations.rds"))
                next # we skip the estimation of any package in that scenario
              }
              ##################################################################
              ##           estimation of GMMs with the EM algorithm           ##
              ##################################################################
              for (index in 1:length(mixture_functions)) {
                # retrieve function values
                mixture_function <- mixture_functions[[index]]$name_fonction; package_name <- names(mixture_functions)[index]
                success_estimation <- TRUE; good_parametrisation <- FALSE # checking that the algorithm is not trapped in boundary space or simply failed
                # message(glue::glue("Package is {package_name}, with following initialisation method {init_algo}."))
                tryCatch(
                  {
                    missing_estimated_theta_list <- do.call(mixture_function, c(
                      x = list(simulated_distribution$x), k = simulated_distribution$k,
                      epsilon = epsilon, itmax = itmax,
                      start = list(initial_estimates), mixture_functions[[index]]$list_params
                    ))
                    good_parametrisation <- check_parameters_validity_multivariate(missing_estimated_theta_list)
                  },
                  error = function(e) {
                    success_estimation <- FALSE; warning(paste("Error in the EM estimation step is:", e))
                  }
                )
                if (!good_parametrisation | !success_estimation) {
                  dir.create("./errors/estimation_failures", showWarnings = F, recursive = T)
                  saveRDS(object = list(simulated_distribution=simulated_distribution, initial_estimates=initial_estimates, package=package_name),
                          file = glue::glue("./errors/estimation_failures/package_{package_name}_init_algo_{init_algo}_step_{t}_{n}_observations.rds"))
                  warning(paste("Estimates trapped in the boundary space or failure of the package", package_name))
                  next # we skip the initialisation in that case
                }
                # normalize returned estimates
                missing_estimated_theta <- missing_estimated_theta_list %>% format_theta_output()

                # store distribution of the estimates
                distribution_parameters_temp <- distribution_parameters_temp %>%
                  dplyr::bind_rows(tibble::tibble(
                    package = package_name, initialisation_method = init_algo,
                    tibble::as_tibble_row(missing_estimated_theta), OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                    entropy = entropy_value, nobservations = n, N.bootstrap = t))

              } # estimation per package
            } # initialization algorithm
            # message(glue::glue("Bootstrap {t} has been performed.\n\n"))
          } # bootstrap repetitions

          #################################################################
          ##     summarize results obtained per scenario configuration   ##
          #################################################################

          distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_temp) %>%
            dplyr::mutate(formatted_true_parameters=list(as.list(formatted_true_theta)), true_parameters=list(as.list(true_theta)),
                          true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"),
                                                        " Nobservations: ", n))

          # summary scores per parameter
          local_scores_temp <- distribution_parameters_temp %>%
            dplyr::group_by(initialisation_method, package) %>%
            dplyr::summarize(get_local_scores(dplyr::cur_data() %>% dplyr::select(dplyr::all_of(bootstrap_colnames)),
                                              formatted_true_theta)) %>%
            tibble::add_column(
              OVL = balanced_ovl, entropy = entropy_value,
              formatted_true_parameters=list(as.list(formatted_true_theta)),
              true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"),
                                            " Nobservations: ", n),
              OVL_pairwise = pairwise_ovl,  nobservations = n)
          local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)
          # generate file name
          filename <- paste0("Nobservations_", n,
                              paste(names(formatted_true_theta), formatted_true_theta, collapse = ": "), ".rds")
          dir.create("./results", showWarnings = F, recursive = T)
          saveRDS(list(distribution=distribution_parameters_temp, local_scores=local_scores_temp),file.path("./results", filename))
          message("One configuration of parameter has been fully completed.\n\n")
        } # number of observations
      }
    }
  } ##### theta configuration

  return(list("distributions" = distribution_parameters, "local_scores" = local_scores))
}


#' @rdname benchmark_univariate_GMM_estimation
#' @export

benchmark_multivariate_GMM_estimation_parallel <- function(mixture_functions, mean_values, proportions, sigma_values,
                                                  cores = getOption("mc.cores", parallel::detectCores()),
                                                  nobservations = c(2000), Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                                  nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                                  initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################

  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble() # store bias and mse for each parameter

  for (p in proportions) {
    for (mu in mean_values) {
      for (sigma in sigma_values) {
        #################################################################
        ##               simulation scenario description               ##
        #################################################################
        true_theta <- list(p = p, mu = mu, sigma = sigma)
        formatted_true_theta <- true_theta %>% format_theta_output()
        k <- length(p); bootstrap_colnames <- names(formatted_true_theta)
        message(paste("We aim at learning the following estimates: \n",
                      paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse=" // ")))
        # bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
        balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)# compute OVL
        pairwise_ovl <- MixSim::overlap(Pi = p, Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy


        for (n in nobservations) {
          filename <- paste0(n, "_observations_entropy", entropy_value, "_OVL_", balanced_ovl)
          distribution_parameters_per_config <- parallel::mclapply(1:Nbootstrap, function(t) {
            simulated_distribution <- simulate_multivariate_GMM(theta = true_theta, n = n) # simulation
            distribution_parameters_per_run <- tibble::tibble()
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##             estimation of the initial estimates             ##
              ##################################################################
              good_initialisation <- TRUE
              tryCatch({
                initial_estimates <- initialize_em_multivariate(
                  x = simulated_distribution$x, k = k, nstart = nstart,
                  short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                )}, error = function(e) {
                  assign(good_initialisation, FALSE, pos=sys.nframe())
                  warning(paste("Error in the initialisation step occurs:", e, "with init algorithm: ", init_algo))
                }
              )
              if(!good_initialisation) {
                dir.create("./errors/initialisation_failures", showWarnings = F, recursive = T)
                saveRDS(object = list(simulated_distribution=simulated_distribution, algo=init_algo),
                        file = glue::glue("./errors/initialisation_failures/init_algo_{init_algo}_step_{t}_{n}_observations.rds"))
                next # we skip the estimation of any package in that scenario
              }
              ##################################################################
              ##           estimation of GMMs with the EM algorithm           ##
              ##################################################################
              for (index in 1:length(mixture_functions)) {
                # retrieve function values
                mixture_function <- mixture_functions[[index]]$name_fonction; package_name <- names(mixture_functions)[index]
                success_estimation <- TRUE; good_parametrisation <- FALSE # checking that the algorithm is not trapped in boundary space or simply failed
                # message(glue::glue("Package is {package_name}, with following initialisation method {init_algo}."))
                tryCatch(
                  {
                    missing_estimated_theta_list <- do.call(mixture_function, c(
                      x = list(simulated_distribution$x), k = simulated_distribution$k,
                      epsilon = epsilon, itmax = itmax,
                      start = list(initial_estimates), mixture_functions[[index]]$list_params
                    ))
                    good_parametrisation <- check_parameters_validity_multivariate(missing_estimated_theta_list)
                  },
                  error = function(e) {
                    success_estimation <- FALSE; warning(paste("Error in the EM estimation step is:", e))
                  }
                )
                if (!good_parametrisation | !success_estimation) {
                  dir.create("./errors/estimation_failures", showWarnings = F, recursive = T)
                  saveRDS(object = list(simulated_distribution=simulated_distribution, initial_estimates=initial_estimates, package=package_name),
                          file = glue::glue("./errors/estimation_failures/package_{package_name}_init_algo_{init_algo}_step_{t}_{n}_observations.rds"))
                  warning(paste("Estimates trapped in the boundary space or failure of the package", package_name))
                  next
                }
                # normalize returned estimates
                missing_estimated_theta <- missing_estimated_theta_list %>% format_theta_output()

                # store distribution of the estimates
                distribution_parameters_per_run <- distribution_parameters_per_run %>% dplyr::bind_rows(
                  tibble::tibble(package = package_name, initialisation_method = init_algo,
                    tibble::as_tibble_row(missing_estimated_theta), OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                    entropy = entropy_value, nobservations = n, N.bootstrap = t))

              } # estimation per package
            } # initialization algorithm
            message(glue::glue("Bootstrap {t} has been performed.\n\n"))
            return(distribution_parameters_per_run)
          }, mc.cores = cores) %>% dplyr::bind_rows()

          #################################################################
          ##     summarize results obtained per scenario configuration   ##
          #################################################################
          distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_per_config) %>%
            dplyr::mutate(formatted_true_parameters=list(as.list(formatted_true_theta)), true_parameters=list(as.list(true_theta)),
                          true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"),
                                                        " Nobservations: ", n))

          # summary scores per parameter
          local_scores_temp <- distribution_parameters_per_config %>%
            dplyr::group_by(initialisation_method, package) %>%
            dplyr::summarize(get_local_scores(dplyr::cur_data() %>% dplyr::select(dplyr::all_of(bootstrap_colnames)),
                                              formatted_true_theta)) %>%
            tibble::add_column(
              OVL = balanced_ovl, entropy = entropy_value,
              formatted_true_parameters=list(as.list(formatted_true_theta)),
              true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"),
                                            " Nobservations: ", n),
              OVL_pairwise = pairwise_ovl,  nobservations = n)
          local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)
          # generate file name
          filename <- paste0("Nobservations_", n,
                             paste(names(formatted_true_theta), formatted_true_theta, collapse = ": "), ".rds")
          dir.create("./results", showWarnings = F, recursive = T)
          saveRDS(list(distribution=distribution_parameters_per_config, local_scores=local_scores_temp),file.path("./results", filename))
          message("One configuration of parameter has been fully completed.\n\n")
        } # number of observations
      }
    }
  } ##### theta configuration

  return(list("distributions" = distribution_parameters, "local_scores" = local_scores))
}


#' @rdname benchmark_univariate_GMM_estimation
#' @export

benchmark_multivariate_GMM_estimation_with_socket <- function(mixture_functions, mean_values, proportions, sigma_values,
                                                              cores = getOption("mc.cores", parallel::detectCores()),
                                                              nobservations = c(2000), Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                                              initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################

  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble() # store bias and mse for each parameter

  for (p in proportions) {
    for (mu in mean_values) {
      for (sigma in sigma_values) {
        #################################################################
        ##               simulation scenario description               ##
        #################################################################
        true_theta <- list(p = p, mu = mu, sigma = sigma)
        formatted_true_theta <- true_theta %>% format_theta_output()
        k <- length(p); bootstrap_colnames <- names(formatted_true_theta)
        message(paste("We aim at learning the following estimates: \n",
                      paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse=" // ")))
        balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)# compute OVL
        pairwise_ovl <- MixSim::overlap(Pi = p, Mu = mu, S = sigma)$BarOmega %>%
          signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy


        for (n in nobservations) {
          filename <- paste0(n, "_observations_entropy", entropy_value, "_OVL_", balanced_ovl)
          cl <- parallel::makeCluster(cores, type="FORK"); doParallel::registerDoParallel(cl)
          parallel::clusterExport(cl, envir=sys.nframe(),
                                  varlist=c("initialisation_algorithms", "mixture_functions", "k",
                                            "epsilon", "itmax", "nstart", "short_iter", "short_eps",
                                            "balanced_ovl", "pairwise_ovl", "entropy_value", "n"))
          distribution_parameters_per_config <- parallel::parLapply(cl, 1:Nbootstrap, function(t) {
            simulated_distribution <- simulate_multivariate_GMM(theta = true_theta, n = n) # simulation
            distribution_parameters_per_run <- tibble::tibble()
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##             estimation of the initial estimates             ##
              ##################################################################
              good_initialisation <- TRUE
              tryCatch({
                initial_estimates <- initialize_em_multivariate(
                  x = simulated_distribution$x, k = k, nstart = nstart,
                  short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                )}, error = function(e) {
                  assign(good_initialisation, FALSE, pos=sys.nframe())
                  warning(paste("Error in the initialisation step occurs:", e, "with init algorithm: ", init_algo))
                }
              )
              if(!good_initialisation) next # we skip the initialisation, as it fails
              ##################################################################
              ##           estimation of GMMs with the EM algorithm           ##
              ##################################################################
              for (index in 1:length(mixture_functions)) {
                # retrieve function values
                mixture_function <- mixture_functions[[index]]$name_fonction; package_name <- names(mixture_functions)[index]
                success_estimation <- TRUE; good_parametrisation <- FALSE # checking that the algorithm is not trapped in boundary space or simply failed
                # message(glue::glue("Package is {package_name}, with following initialisation method {init_algo}."))
                tryCatch(
                  {
                    missing_estimated_theta_list <- do.call(mixture_function, c(
                      x = list(simulated_distribution$x), k = simulated_distribution$k,
                      epsilon = epsilon, itmax = itmax,
                      start = list(initial_estimates), mixture_functions[[index]]$list_params
                    ))
                    good_parametrisation <- check_parameters_validity_multivariate(missing_estimated_theta_list)
                  },
                  error = function(e) {
                    success_estimation <- FALSE; warning(paste("Error in the EM estimation step is:", e))
                  }
                )
                if (!good_parametrisation | !success_estimation) {
                  dir.create("./errors/estimation_failures", showWarnings = F, recursive = T)
                  saveRDS(object = list(simulated_distribution=simulated_distribution, initial_estimates=initial_estimates, package=package_name),
                          file = glue::glue("./errors/estimation_failures/package_{package_name}_init_algo_{init_algo}_step_{t}_{n}_observations.rds"))
                  warning(paste("Estimates trapped in the boundary space or failure of the package", package_name))
                  next
                }
                # normalize returned estimates
                missing_estimated_theta <- missing_estimated_theta_list %>% format_theta_output()

                # store distribution of the estimates
                distribution_parameters_per_run <- distribution_parameters_per_run %>% dplyr::bind_rows(
                  tibble::tibble(package = package_name, initialisation_method = init_algo,
                                 tibble::as_tibble_row(missing_estimated_theta), OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                                 entropy = entropy_value, nobservations = n, N.bootstrap = t))

              } # estimation per package
            } # initialization algorithm
            return(distribution_parameters_per_run)
          }) %>% dplyr::bind_rows()

          parallel::stopCluster(cl) # stop parallel computation
          #################################################################
          ##     summarize results obtained per scenario configuration   ##
          #################################################################
          distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_per_config) %>%
            dplyr::mutate(formatted_true_parameters=list(as.list(formatted_true_theta)), true_parameters=list(as.list(true_theta)),
                          true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"),
                                                        " Nobservations: ", n))

          # summary scores per parameter
          local_scores_temp <- distribution_parameters_per_config %>%
            dplyr::group_by(initialisation_method, package) %>%
            dplyr::summarize(get_local_scores(dplyr::cur_data() %>% dplyr::select(dplyr::all_of(bootstrap_colnames)),
                                              formatted_true_theta)) %>%
            tibble::add_column(
              OVL = balanced_ovl, entropy = entropy_value,
              formatted_true_parameters=list(as.list(formatted_true_theta)),
              true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";")),
              OVL_pairwise = pairwise_ovl,  nobservations = n)
          local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)
          # generate file name
          filename <- paste0("Nobservations_", n,
                             paste(names(formatted_true_theta), formatted_true_theta, collapse = ": "), ".rds")
          dir.create("./results", showWarnings = F, recursive = T)
          saveRDS(list(distribution=distribution_parameters_per_config, local_scores=local_scores_temp),file.path("./results", filename))
          message("One configuration of parameter has been fully completed.\n\n")
        } # number of observations
      }
    }
  } ##### theta configuration

  return(list("distributions" = distribution_parameters, "local_scores" = local_scores))
}










#' Launch the benchmark to compare computational performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param sigma_values,mean_values,proportions the true parameters to be retrieved
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


compute_microbenchmark_univariate <- function(mixture_functions,
                                   sigma_values, mean_values, proportions,
                                   prop_outliers = 0, nobservations = c(100, 1000, 10000),
                                   Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                   nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                   initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {

  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  time_data <- tibble::tibble()
  init_time_data <- tibble::tibble()
  for (prop_out in prop_outliers) {
    for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) {
            #################################################################
            ##               simulation scenario description               ##
            #################################################################
            true_theta <- true_theta;  k <- length(p) # number of components
            bootstrap_colnames <- names(unlist(true_theta[c("p", "mu", "sigma")])) # labels used for naming the parameters
            balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega %>%
              signif(digits = 2)# compute OVL
            pairwise_ovl <- MixSim::overlap(Pi = true_theta$p, Mu = as.matrix(true_theta$mu), S = as.matrix(true_theta$sigma))$BarOmega %>%
              signif(digits = 2)
            entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy

            for (n in nobservations) {
              filename <- paste0(
                n, "_observations_entropy", entropy_value,
                "_OVL_", balanced_ovl, "_prop_outliers_", prop_out
              )
              simulated_distribution <- simulate_univariate_GMM(n = n, theta = true_theta, prop_outliers = prop_out, interval = 2) # simulation of the experience
              mbm_temp <- tibble::tibble(); init_temp <- tibble::tibble()

              for (t in 1:Nbootstrap) {
                for (init_algo in initialisation_algorithms) {
                  ##################################################################
                  ##        estimation of time taken for the initialisation       ##
                  ##################################################################
                  initial_estimates <- initialize_em_univariate(
                    x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  )
                  # and determine time taken for the initialisation step
                  init_temp <- init_temp %>% dplyr::bind_rows(
                    microbenchmark::microbenchmark(initialisation_method = initialize_em_univariate(
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
                        epsilon = epsilon, itmax = itmax, initial_estimates = list(initial_estimates), mixture_functions[[index]]$list_params
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

              ##################################################################
              ##                    save time computations                    ##
              ##################################################################


              init_temp <- init_temp %>% dplyr::mutate(
                entropy = entropy_value,
                OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,nobservations = n, prop_outliers = prop_out)
              init_time_data <- init_time_data %>% dplyr::bind_rows(init_temp) # store time computation taken by the initialisation step

              mbm_temp <- mbm_temp %>% dplyr::mutate(
                entropy = entropy_value,
                OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                nobservations = n, prop_outliers = prop_out
              )
              time_data <- time_data %>% dplyr::bind_rows(mbm_temp) # store time computation taken by the EM algorithm
            } # number of observations
          } # theta configuration
        }
      }
  }
  return(list(time_data = time_data, init_time_data = init_time_data))
}


#' @rdname compute_microbenchmark_univariate
#' @export


compute_microbenchmark_multivariate <- function(mixture_functions,
                                              sigma_values, mean_values, proportions,
                                              cores = getOption("mc.cores", parallel::detectCores()),
                                              nobservations = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
                                              Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                              initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {

  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  time_data <- tibble::tibble(); init_time_data <- tibble::tibble()
  for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) {
        #################################################################
        ##               simulation scenario description               ##
        #################################################################
        true_theta <- list(p=p, mu=mu, sigma=sigma)
        formatted_true_theta <- true_theta %>% format_theta_output()
        k <- length(p); bootstrap_colnames <- names(formatted_true_theta)
        balanced_ovl <- MixSim::overlap(Pi = rep(1 / k, k), Mu = true_theta$mu, S = true_theta$sigma)$BarOmega %>%
          signif(digits = 2)# compute OVL
        pairwise_ovl <- MixSim::overlap(Pi = true_theta$p, Mu = true_theta$mu, S = true_theta$sigma)$BarOmega %>%
          signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy

        for (n in nobservations) {
          simulated_distribution <- simulate_multivariate_GMM(n = n, theta = true_theta) # simulation of the experience
          time_configurations <- parallel::mclapply(1:Nbootstrap, function(t) {
            init_temp <- tibble::tibble();  mbm_temp <- tibble::tibble() # store intemediate computations
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##        estimation of time taken for the initialisation       ##
              ##################################################################
              good_initialisation <- tryCatch({
                initial_estimates <- initialize_em_multivariate(
                  x = simulated_distribution$x, k = k, nstart = nstart,
                  short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                )}, error = function(e) {e})

              # if(inherits(good_initialisation, "error")) next # we skip the initialisation, as it fails
              if(!inherits(good_initialisation, "error")) {
                # time taken by the initialisation step
                init_temp <- init_temp %>% dplyr::bind_rows(
                  microbenchmark::microbenchmark(initialisation_method = initialize_em_multivariate(
                    x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  ), times = 1) %>%
                    tibble::as_tibble() %>% dplyr::rename(initialisation_method = expr) %>%
                    dplyr::mutate(time = microbenchmark:::convert_to_unit(time, "s"), initialisation_method = init_algo, N.bootstrap = t))

                # perform micro-benchmark
                for (index in 1:length(mixture_functions)) {
                  ##################################################################
                  ##        estimation of time taken by the EM algorithm          ##
                  ##################################################################

                  mixture_function <- mixture_functions[[index]]$name_fonction
                  package_name <- names(mixture_functions)[index]
                  mbm_temp <- mbm_temp %>% dplyr::bind_rows(microbenchmark::microbenchmark(package = do.call(mixture_function, c(
                      x = list(simulated_distribution$x), k = simulated_distribution$k,
                      epsilon = epsilon, itmax = itmax, initial_estimates = list(initial_estimates), mixture_functions[[index]]$list_params
                    )), times = 1) %>%
                      tibble::as_tibble() %>% dplyr::rename(package = expr) %>%
                      dplyr::mutate(time = microbenchmark:::convert_to_unit(time, "s"), initialisation_method = init_algo,
                        N.bootstrap = t, package = package_name))
                } # mixture package
              } # check error initialization algorithm


            } # initialization algorithm
            return(list(mbm_temp=mbm_temp, init_temp=init_temp))
          }, mc.cores=cores) # repeated Bootstraps

          ##################################################################
          ##                    save time computations                    ##
          ##################################################################

          # store time computation taken by the initialisation step
          init_time_data <- init_time_data %>% dplyr::bind_rows(
            time_configurations %>% purrr::map_dfr("init_temp") %>%
              dplyr::mutate(entropy = entropy_value, OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                            nobservations = n, formatted_true_parameters=list(as.list(formatted_true_theta)),
                            true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"))))
          # store time computation taken by the estimation part
          time_data <- time_data %>% dplyr::bind_rows(
            time_configurations %>% purrr::map_dfr ("mbm_temp") %>%
              dplyr::mutate(entropy = entropy_value, OVL = balanced_ovl, OVL_pairwise = pairwise_ovl,
                            nobservations = n, formatted_true_parameters=list(as.list(formatted_true_theta)),
                            true_parameters_factor=paste0(paste0(bootstrap_colnames, ": ", formatted_true_theta, collapse = ";"))))
        } # number of observations
      } # theta configuration
    }
  }
  return(list(time_data = time_data, init_time_data = init_time_data))
}
