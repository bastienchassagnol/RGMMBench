# logLikelihood <- EMCluster::logL(as.matrix(simulated_distribution$x),
#                                  emobj = list(
#                                    pi = missing_estimated_theta_list$p, Mu = as.matrix(missing_estimated_theta_list$mu),
#                                    LTSigma = as.matrix(missing_estimated_theta_list$sigma^2)
#                                  ))


# # aggregate scores obtained per parameter
# global_bias <- local_scores_temp %>%
#   dplyr::filter(scores == "bias") %>%
#   dplyr::rowwise() %>%
#   dplyr::transmute(global_bias = sum(abs(dplyr::c_across(dplyr::all_of(bootstrap_colnames))))) # global bias
# global_mse <- local_scores_temp %>%
#   dplyr::filter(scores == "mse") %>%
#   dplyr::rowwise() %>%
#   dplyr::transmute(global_mse = sum(dplyr::c_across(dplyr::all_of(bootstrap_colnames)))) # global mse
# missed_cases <- distribution_parameters_temp %>%
#   dplyr::count(initialisation_method, package, name = "N.missed") %>%
#   dplyr::mutate(N.missed = Nbootstrap - N.missed) # count number of times the package failed in the estimation with the given number of clusters
#
# global_scores_temp <- global_scores_temp %>%
#   dplyr::group_by(initialisation_method, package) %>%
#   dplyr::summarise(
#     logLikelihood = mean(logLikelihood),
#     entropy = signif(entropy_value, digits = 2),
#     OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2), skew = skew[1], nobservations = n, prop_outliers = prop_out
#   ) %>%
#   dplyr::inner_join(missed_cases, by = c("initialisation_method", "package")) %>%
#   dplyr::inner_join(global_bias, by = c("initialisation_method", "package")) %>%
#   dplyr::inner_join(global_mse, by = c("initialisation_method", "package"))
# global_scores <- global_scores %>% dplyr::bind_rows(global_scores_temp)
