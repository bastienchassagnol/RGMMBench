#' Launch the benchmark to compare statistical performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param id_scenario Possibility to set it to another number than one, to uniquely identify them
#' @param sigma_values,mean_values,proportions the true parameters to be retrieved
#' @param prop_outliers the proportion of outliers added in the simulation
#' @param nobservations the number of observations drawn to generate the random sample
#' @param Nbootstrap the number of bootstrap simulations and repetitions to perform
#' @param epsilon,itmax respectively criterion threshold and maximal number of iterations to reach it
#' @param nstart,short_iter,short_eps hyper-parameters to control the initialisation step
#' @param prior_prob add minimal uncertainty on the cluster assignment returned by hierarchical clustering method
#' @param initialisation_algorithms among 6 methods, which algorithms to be chosen for the initialisation phase
#' @param cores the number of cores to be used, by default all the available cores
#'
#' @return a list with the simulated distributions of the estimates, some summary scores per parameter and aggregated measures
#' as well as boxplot and Heatmap correlation representations of the estimates
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang ":="
#' @importFrom rlang .data
#' @import ggplot2
#'
#' @export

benchmark_univariate_GMM_estimation <- function(mixture_functions, sigma_values, mean_values, proportions,
                                                cores = getOption("mc.cores", parallel::detectCores()), id_scenario=NULL,
                                                prop_outliers = 0, nobservations = c(2000),
                                                Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                                nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                                initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  id_tibble <- tibble::tibble()
  if (is.null(id_scenario)) id_scenario <- 1 # to uniquely identify each run
  distribution_parameters <- tibble::tibble() # store empirical bootstrap distribution of the estimates
  local_scores <- tibble::tibble() # store bias and mse for each parameter


  for (prop_out in prop_outliers) {
    for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) {
          #################################################################
          ##               simulation scenario description               ##
          #################################################################
          true_theta <- list(p = p, mu = mu, sigma = sigma)
          formatted_true_theta <- true_theta %>% format_theta_output()
          k <- length(p)
          bootstrap_colnames <- names(formatted_true_theta)
          balanced_ovl <- compute_average_overlap(true_theta %>% magrittr::inset2("p", rep(1 / k, k))) %>% signif(digits = 2)
          pairwise_ovl <- compute_average_overlap(true_theta) %>% signif(digits = 2)
          entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy

          for (i in seq_along(nobservations)) {
            n <- nobservations[i]; letter_id <- letters[i]
            message(paste("We are at scenario ID:", paste0(id_scenario, letter_id), ".\n"))

            distribution_parameters_per_config <- parallel::mclapply(1:Nbootstrap, function(t) {
              distribution_parameters_per_run <- tibble::tibble()
              simulated_distribution <- simulate_univariate_GMM(n = n, theta = true_theta, prop_outliers = prop_out, interval = 2) # simulation
              for (init_algo in initialisation_algorithms) {
                ##################################################################
                ##             estimation of the initial estimates             ##
                ##################################################################
                good_initialisation <- tryCatch(
                  {
                    initial_estimates <- initialize_em_univariate(
                      x = simulated_distribution$x, k = k, nstart = nstart,
                      short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                    )
                  },
                  error = function(e) {
                    e
                  }
                )
                if (!inherits(good_initialisation, "error")) {
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

                    success_estimation <- tryCatch(
                      {
                        missing_estimated_theta_list <- do.call(mixture_function, c(
                          x = list(simulated_distribution$x), k = simulated_distribution$k,
                          epsilon = epsilon, itmax = itmax,
                          start = list(initial_estimates), mixture_functions[[index]]$list_params
                        ))
                        stopifnot(
                          "Parameters provided do not correspond to a fitted GMM estimation" =
                            check_parameters_validity_univariate(missing_estimated_theta_list) == TRUE
                        )
                      },
                      error = function(e) {
                        return(e)
                      }
                    )
                    if (!inherits(success_estimation, "error")) {
                      missing_estimated_theta_list <- missing_estimated_theta_list[c("p", "mu", "sigma")]
                      missing_estimated_theta <- stats::setNames(unlist(missing_estimated_theta_list), bootstrap_colnames)
                      distribution_parameters_per_run <- distribution_parameters_per_run %>%
                        dplyr::bind_rows(tibble::tibble(
                          package = package_name, initialisation_method = init_algo,
                          tibble::as_tibble_row(missing_estimated_theta), N.bootstrap = t
                        ))
                    }
                  } # estimation per package
                }
              } # initialization algorithm
              return(distribution_parameters_per_run)
            }, mc.cores = cores) %>% dplyr::bind_rows() # one config terminated

            #################################################################
            ##     summarize results obtained per scenario configuration   ##
            #################################################################
            iid_tibble_temp <- tibble::tibble(
              ID = paste0(id_scenario, letter_id), OVL = balanced_ovl, entropy = entropy_value,
              OVL_pairwise = pairwise_ovl, nobservations = n, prop_outliers = 0,
              formatted_true_parameters = list(as.list(formatted_true_theta)),
              true_parameters = list(as.list(true_theta)))
            id_tibble <- id_tibble %>% dplyr::bind_rows(id_tibble_temp)


            distribution_parameters_per_config <- distribution_parameters_per_config %>% tibble::add_column(ID = paste0(id_scenario, letter_id))
            distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_per_config)


            local_scores_temp <- distribution_parameters_per_config %>%
              dplyr::group_by(dplyr::across(c("initialisation_method", "package"))) %>%
              dplyr::summarize(get_local_scores(dplyr::cur_data() %>%
                                                  dplyr::select(dplyr::all_of(bootstrap_colnames)), formatted_true_theta)) %>%
              tibble::add_column(ID = paste0(id_scenario, letter_id))
            local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)

            # store temporary results
            dir.create("./results", showWarnings = F, recursive = T)
            saveRDS(list(
              distribution = distribution_parameters_per_config,
              local_scores = local_scores_temp,
              config = id_tibble_temp),
              file.path("./results", paste0("ID_scenario_", paste0(id_scenario, letter_id), ".rds")))
          } # number of observations
          message(paste("Scenario", id_scenario, "has been achieved, with the whole setting of observations.\n"))
          id_scenario <- id_scenario + 1
        }
      }
    }
  } ##### theta configuration (p, mu, sd, prop_outliers)
  return(list("distributions" = distribution_parameters, "local_scores" = local_scores, "config" = id_tibble))
}



#' @rdname benchmark_univariate_GMM_estimation
#' @export

benchmark_multivariate_GMM_estimation <- function(mixture_functions, mean_values, proportions, sigma_values,
                                                  id_scenario = NULL, cores = getOption("mc.cores", parallel::detectCores()),
                                                  nobservations = c(2000), Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                                  nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                                  initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {



  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  id_tibble <- tibble::tibble();
  if (is.null(id_scenario)) id_scenario <- 1 # to uniquely identify each run
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
        k <- length(p)
        bootstrap_colnames <- names(formatted_true_theta)
        balanced_ovl <- MixSim::overlap(Pi=rep(1/k, k), Mu=t(true_theta$mu), S=as.array(true_theta$sigma))$MaxOmega %>% signif(digits = 2)
        pairwise_ovl <- MixSim::overlap(Pi=true_theta$p, Mu=t(true_theta$mu), S=as.array(true_theta$sigma))$MaxOmega %>% signif(digits = 2)
        # balanced_ovl <- compute_average_overlap(true_theta %>% magrittr::inset2("p", rep(1 / k, k))) %>% signif(digits = 2)
        # pairwise_ovl <- compute_average_overlap(true_theta) %>% signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy


        for (i in seq_along(nobservations)) {
          n <- nobservations[i]; letter_id <- letters[i]
          message(paste("We are at scenario ID:", paste0(id_scenario, letter_id), ".\n"))
          # distribution_parameters_per_config <- tibble::tibble()
          # for (t in 1:Nbootstrap) {
          distribution_parameters_per_config <- parallel::mclapply(1:Nbootstrap, function(t) {
            simulated_distribution <- simulate_multivariate_GMM(theta = true_theta, n = n) # simulation
            distribution_parameters_per_run <- tibble::tibble()
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##             estimation of the initial estimates             ##
              ##################################################################
              good_initialisation <- tryCatch(
                {
                  initial_estimates <- initialize_em_multivariate(
                    x = simulated_distribution$x, k = k, nstart = nstart,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  )
                },
                error = function(e) {
                  return(e)
                }
              )

              if (!inherits(good_initialisation, "error")) {
                ##################################################################
                ##           estimation of GMMs with the EM algorithm           ##
                ##################################################################
                for (index in 1:length(mixture_functions)) {
                  # retrieve function values
                  mixture_function <- mixture_functions[[index]]$name_fonction
                  package_name <- names(mixture_functions)[index]
                  success_estimation <- tryCatch(
                    { # use of ternary operator structure
                      start <- if (any(isTRUE(all.equal(mixture_function, em_clustvarsel_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_pgmm_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_HDclassif_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_EMMIXmfa_multivariate)))) NULL else  list(initial_estimates)

                      missing_estimated_theta_list <- do.call(mixture_function, c(
                        x = list(simulated_distribution$x), k = simulated_distribution$k,
                        epsilon = epsilon, itmax = itmax, initialisation_algorithm=init_algo,
                        start = start, mixture_functions[[index]]$list_params
                      ))
                      stopifnot(
                        "Parameters provided do not correspond to a fitted GMM estimation" =
                          check_parameters_validity_multivariate(missing_estimated_theta_list) == TRUE
                      )
                    },
                    error = function(e) {
                      dir.create("./errors", showWarnings = FALSE, recursive = T)
                      saveRDS(list(x=simulated_distribution$x, k= simulated_distribution$k,
                                   epsilon = epsilon, itmax = itmax, start=start, error=e),
                              file = paste0("./errors/scenario_", id_scenario, "_init_algo_", init_algo,
                                            "_package_name_", package_name, "_bootstrap_", t,".rds"))
                      return(e)
                    }
                  )
                  if (!inherits(success_estimation, "error")) {
                    # normalize returned estimates
                    missing_estimated_theta <- missing_estimated_theta_list %>% format_theta_output()

                    # store distribution of the estimates
                    distribution_parameters_per_run <- distribution_parameters_per_run %>% dplyr::bind_rows(
                      tibble::tibble(
                        package = package_name, initialisation_method = init_algo,
                        tibble::as_tibble_row(missing_estimated_theta), N.bootstrap = t
                      )
                    )
                  }
                } # estimation per package
              } # if initialization is properly done
            } # initialization algorithm

          # distribution_parameters_per_config <- distribution_parameters_per_config %>%
          #   dplyr::bind_rows(distribution_parameters_per_run)
          # } # Bootstraps
          return(distribution_parameters_per_run)
          }, mc.cores = cores) %>% dplyr::bind_rows()

          #################################################################
          ##     summarize results obtained per scenario configuration   ##
          #################################################################
          id_tibble_temp <- tibble::tibble(
            ID = paste0(id_scenario, letter_id), OVL = balanced_ovl, entropy = entropy_value,
            OVL_pairwise = pairwise_ovl, nobservations = n, prop_outliers = 0,
            formatted_true_parameters = list(as.list(formatted_true_theta)),
            true_parameters = list(as.list(true_theta)))
          id_tibble <- id_tibble %>% dplyr::bind_rows(id_tibble_temp)


          distribution_parameters_per_config <- distribution_parameters_per_config %>% tibble::add_column(ID = paste0(id_scenario, letter_id))
          distribution_parameters <- distribution_parameters %>% dplyr::bind_rows(distribution_parameters_per_config)


          local_scores_temp <- distribution_parameters_per_config %>%
            dplyr::group_by(dplyr::across(c("initialisation_method", "package"))) %>%
            dplyr::summarize(get_local_scores(dplyr::cur_data() %>%
                                                dplyr::select(dplyr::all_of(bootstrap_colnames)), formatted_true_theta)) %>%
            tibble::add_column(ID = paste0(id_scenario, letter_id))
          local_scores <- local_scores %>% dplyr::bind_rows(local_scores_temp)

          # store temporary results
          dir.create("./results", showWarnings = F, recursive = T)
          saveRDS(list(
              distribution = distribution_parameters_per_config,
              local_scores = local_scores_temp,
              config = id_tibble_temp),
            file.path("./results", paste0("ID_scenario_", paste0(id_scenario, letter_id), ".rds")))
        } # number of observations
        message(paste("Scenario", id_scenario, "has been achieved, with the whole setting of observations.\n"))
        id_scenario <- id_scenario + 1
      }
    }
  } ##### theta configuration
  return(list("distributions" = distribution_parameters, "local_scores" = local_scores, "config" = id_tibble))
}





#' Launch the benchmark to compare computational performances between packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param mixture_functions List of the packages to be compared (Id:name of the package, value: its options)
#' @param id_scenario Possibility to set it to another number than one, to uniquely identify them
#' @param sigma_values,mean_values,proportions the true parameters to be retrieved
#' @param prop_outliers the proportion of outliers added in the simulation
#' @param nobservations the number of observations drawn to generate the random sample
#' @param Nbootstrap the number of bootstrap simulations and repetitions to perform
#' @param epsilon,itmax respectively criterion threshold and maximal number of iterations to reach it
#' @param nstart,short_iter,short_eps hyper-parameters to control the initialisation step
#' @param prior_prob add minimal uncertainty on the cluster assignment returned by hierarchical clustering method
#' @param initialisation_algorithms among 6 methods, which algorithms to be chosen for the initialisation phase
#' @param cores the number of cores to be used, by default all the available cores
#'
#' @return a list with the running time of the initialisation and the EM estimation itself, as well as corresponding time curve representations
#'
#' @export


compute_microbenchmark_univariate <- function(mixture_functions, id_scenario = NULL,
                                              sigma_values, mean_values, proportions,
                                              cores = getOption("mc.cores", parallel::detectCores()),
                                              prop_outliers = 0, nobservations = c(100, 1000, 10000),
                                              Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                              nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                              initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")) {

  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  id_tibble <- tibble::tibble()
  id_scenario <- 1
  time_data <- tibble::tibble()
  init_time_data <- tibble::tibble()
  for (prop_out in prop_outliers) {
    for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) { #################################################################
          ##               simulation scenario description               ##
          #################################################################
          true_theta <- list(p = p, mu = mu, sigma = sigma)
          formatted_true_theta <- true_theta %>% format_theta_output()
          k <- length(p)
          bootstrap_colnames <- names(formatted_true_theta)
          balanced_ovl <- MixSim::overlap(Pi=rep(1/k, k), Mu=t(true_theta$mu), S=as.array(true_theta$sigma))$MaxOmega %>%
            signif(digits = 2)
          pairwise_ovl <- MixSim::overlap(Pi=true_theta$p, Mu=t(true_theta$mu), S=as.array(true_theta$sigma))$MaxOmega %>%
            signif(digits = 2)
          entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy

          for (i in seq_along(nobservations)) {
            n <- nobservations[i]; letter_id <- letters[i]
            message(paste("We are at scenario ID:", paste0(id_scenario, letter_id), ".\n"))
            simulated_distribution <- simulate_univariate_GMM(n = n, theta = true_theta) # simulation of the experience
            time_configurations <- parallel::mclapply(1:Nbootstrap, function(t) {
              init_temp <- tibble::tibble()
              mbm_temp <- tibble::tibble() # store intermediate computations
              for (init_algo in initialisation_algorithms) {
                ##################################################################
                ##        estimation of time taken for the initialisation       ##
                ##################################################################
                good_initialisation <- tryCatch(
                  {
                    initial_estimates <- initialize_em_univariate(
                      x = simulated_distribution$x, k = k, nstart = nstart,
                      short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                    )
                  },
                  error = function(e) {
                    e
                  }
                )

                if (!inherits(good_initialisation, "error")) {
                  # time taken by the initialisation step
                  init_temp <- init_temp %>% dplyr::bind_rows(
                    microbenchmark::microbenchmark(initialisation_method = initialize_em_univariate(
                      x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                      short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                    ), times = 1) %>%
                      tibble::as_tibble() %>% dplyr::rename(initialisation_method = expr) %>%
                      dplyr::mutate(time = convert_to_unit(.data$time, "s"), initialisation_method = init_algo, N.bootstrap = t)
                  )

                  # perform micro-benchmark
                  for (index in 1:length(mixture_functions)) {
                    ##################################################################
                    ##        estimation of time taken by the EM algorithm          ##
                    ##################################################################

                    mixture_function <- mixture_functions[[index]]$name_fonction
                    package_name <- names(mixture_functions)[index]
                    good_estimation <- tryCatch(
                      {
                        mbm_temp_per_function <- microbenchmark::microbenchmark(package = do.call(mixture_function, c(
                          x = list(simulated_distribution$x), k = simulated_distribution$k,
                          epsilon = epsilon, itmax = itmax, initial_estimates = list(initial_estimates), mixture_functions[[index]]$list_params
                        )), times = 1) %>%
                          tibble::as_tibble() %>%
                          dplyr::rename(package = expr) %>%
                          dplyr::mutate(
                            time = convert_to_unit(.data$time, "s"), initialisation_method = init_algo,
                            N.bootstrap = t, package = package_name
                          )
                      },
                      error = function(e) {
                        e
                      }
                    )
                    if (!inherits(good_estimation, "error")) {
                      mbm_temp <- mbm_temp %>% dplyr::bind_rows(mbm_temp_per_function)
                    }
                  } # mixture package
                } # check error initialization algorithm
              } # initialization algorithm
              return(list(mbm_temp = mbm_temp, init_temp = init_temp))
            }, mc.cores = cores) # repeated Bootstraps

            ##################################################################
            ##                    save time computations                    ##
            ##################################################################

            id_tibble <- id_tibble %>% dplyr::bind_rows(tibble::tibble(
              ID = id_scenario, OVL = balanced_ovl, entropy = entropy_value,
              OVL_pairwise = pairwise_ovl, prop_outliers = 0, nobservations = n,
              formatted_true_parameters = list(as.list(formatted_true_theta)),
              true_parameters = list(as.list(true_theta))
            ))

            # store time computation taken by the initialisation step
            init_time_data_temp <- time_configurations %>%
              purrr::map_dfr("init_temp") %>%
              dplyr::mutate(ID = id_scenario, nobservations = n)
            init_time_data <- init_time_data %>% dplyr::bind_rows(init_time_data_temp)

            # store time computation taken by the estimation part
            time_data_temp <- time_configurations %>%
              purrr::map_dfr("mbm_temp") %>%
              dplyr::mutate(ID = id_scenario, nobservations = n)
            time_data <- time_data %>% dplyr::bind_rows(time_data_temp)
          } # number of observations

          # store temporary results (with all observations gathered this time)
          dir.create("./results", showWarnings = F, recursive = T)
          saveRDS(list(
            init_time_data = init_time_data %>% dplyr::filter(.data$ID == id_scenario),
            time_data = time_data %>% dplyr::filter(.data$ID == id_scenario),
            config = id_tibble %>% dplyr::filter(.data$ID == id_scenario)),
            file.path("./results", paste0("time_computation_ID_scenario_", id_scenario, ".rds")))

          message(paste("Scenario", id_scenario, "has been achieved, with the whole setting of observations.\n"))
          id_scenario <- id_scenario + 1
        } # theta configuration
      }
    }
  }
  return(list(time_data = time_data, init_time_data = init_time_data, "config" = id_tibble))
}


#' @rdname compute_microbenchmark_univariate
#' @export


compute_microbenchmark_multivariate <- function(mixture_functions, id_scenario = NULL,
                                                sigma_values, mean_values, proportions,
                                                cores = getOption("mc.cores", parallel::detectCores()),
                                                nobservations = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
                                                Nbootstrap = 100, epsilon = 10^-6, itmax = 1000,
                                                nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                                initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")) {

  #################################################################
  ##     name variables for storing parameters distribution      ##
  #################################################################
  id_tibble <- tibble::tibble();
  if (is.null(id_scenario)) id_scenario <- 1 # to uniquely identify each run
  time_data <- tibble::tibble();   init_time_data <- tibble::tibble()
  for (p in proportions) {
    for (mu in mean_values) {
      for (sigma in sigma_values) {
        #################################################################
        ##               simulation scenario description               ##
        #################################################################
        true_theta <- list(p = p, mu = mu, sigma = sigma)
        formatted_true_theta <- true_theta %>% format_theta_output()
        k <- length(p)
        bootstrap_colnames <- names(formatted_true_theta)
        balanced_ovl <- compute_average_overlap(true_theta %>% magrittr::inset2("p", rep(1 / k, k))) %>% signif(digits = 2)
        pairwise_ovl <- compute_average_overlap(true_theta) %>% signif(digits = 2)
        entropy_value <- compute_shannon_entropy(p) %>% signif(digits = 2) # compute entropy
        for (i in seq_along(nobservations)) {
          n <- nobservations[i]; letter_id <- letters[i]
          message(paste("We are at scenario ID:", paste0(id_scenario, letter_id), ".\n"))
          simulated_distribution <- simulate_multivariate_GMM(n = n, theta = true_theta) # simulation of the experience
          time_configurations <- parallel::mclapply(1:Nbootstrap, function(t) {
            init_temp <- tibble::tibble(); mbm_temp <- tibble::tibble() # store intermediate computations
            for (init_algo in initialisation_algorithms) {
              ##################################################################
              ##        estimation of time taken for the initialisation       ##
              ##################################################################
              good_initialisation <- tryCatch(
                {
                  initial_estimates <- initialize_em_multivariate(
                    x = simulated_distribution$x, k = k, nstart = nstart,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  )
                },
                error = function(e) {
                  return(e)
                }
              )

              if (!inherits(good_initialisation, "error")) {
                # time taken by the initialisation step
                init_temp <- init_temp %>% dplyr::bind_rows(
                  microbenchmark::microbenchmark(initialisation_method = initialize_em_multivariate(
                    x = simulated_distribution$x, k = k, nstart = nstart, prior_prob = prior_prob,
                    short_iter = short_iter, short_eps = short_eps, initialisation_algorithm = init_algo
                  ), times = 1) %>%
                    tibble::as_tibble() %>% dplyr::rename(initialisation_method = expr) %>%
                    dplyr::mutate(time = convert_to_unit(.data$time, "s"), initialisation_method = init_algo, N.bootstrap = t)
                )

                # perform micro-benchmark
                for (index in 1:length(mixture_functions)) {
                  ##################################################################
                  ##        estimation of time taken by the EM algorithm          ##
                  ##################################################################

                  mixture_function <- mixture_functions[[index]]$name_fonction
                  package_name <- names(mixture_functions)[index]
                  good_estimation <- tryCatch(
                    {
                      start <- if (any(isTRUE(all.equal(mixture_function, em_clustvarsel_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_pgmm_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_HDclassif_multivariate)),
                                       isTRUE(all.equal(mixture_function, em_EMMIXmfa_multivariate)))) NULL else  list(initial_estimates)

                      mbm_temp_per_function <- microbenchmark::microbenchmark(package = do.call(mixture_function, c(
                        x = list(simulated_distribution$x), k = simulated_distribution$k,
                        epsilon = epsilon, itmax = itmax, initialisation_algorithm=init_algo,
                        start = start, mixture_functions[[index]]$list_params
                      )), times = 1) %>%
                        tibble::as_tibble() %>%
                        dplyr::rename(package = expr) %>%
                        dplyr::mutate(
                          time = convert_to_unit(.data$time, "s"), initialisation_method = init_algo,
                          N.bootstrap = t, package = package_name
                        )
                    },
                    error = function(e) {
                      e
                    }
                  )
                  if (!inherits(good_estimation, "error")) {
                    mbm_temp <- mbm_temp %>% dplyr::bind_rows(mbm_temp_per_function)
                  }
                } # mixture package
              } # check error initialization algorithm
            } # initialization algorithm
            return(list(mbm_temp = mbm_temp, init_temp = init_temp))
          }, mc.cores = cores) # repeated Bootstraps

          ##################################################################
          ##                    save time computations                    ##
          ##################################################################

          # store time computation taken by the initialisation step
          init_time_data_temp <- time_configurations %>%
            purrr::map_dfr("init_temp") %>%
            dplyr::mutate(ID = id_scenario, nobservations = n)
          init_time_data <- init_time_data %>% dplyr::bind_rows(init_time_data_temp)

          # store time computation taken by the estimation part
          time_data_temp <- time_configurations %>%
            purrr::map_dfr("mbm_temp") %>%
            dplyr::mutate(ID = id_scenario, nobservations = n)
          time_data <- time_data %>% dplyr::bind_rows(time_data_temp)
          message(paste("\nScenario ID:", id_scenario,"with n=", n, "has been carried out.\n"))

        } # number of observations

        id_tibble <- id_tibble %>% dplyr::bind_rows(tibble::tibble(
          ID = id_scenario, OVL = balanced_ovl, entropy = entropy_value,
          OVL_pairwise = pairwise_ovl, prop_outliers = 0, nobservations = list(nobservations),
          formatted_true_parameters = list(as.list(formatted_true_theta)),
          true_parameters = list(as.list(true_theta))))

        # store temporary results (with all observations gathered this time)
        dir.create("./results", showWarnings = F, recursive = T)
        saveRDS(list(
          init_time_data = init_time_data %>% dplyr::filter(.data$ID == id_scenario),
          time_data = time_data %>% dplyr::filter(.data$ID == id_scenario),
          config = id_tibble %>% dplyr::filter(.data$ID == id_scenario)),
          file.path("./results", paste0("time_computation_ID_scenario_", id_scenario, ".rds")))

        message(paste("Scenario", id_scenario, "has been achieved, with the whole setting of observations.\n"))
        id_scenario <- id_scenario + 1
      } # theta configuration
    }
  }
  return(list(time_data = time_data, init_time_data = init_time_data, "config" = id_tibble))
}
