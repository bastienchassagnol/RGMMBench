
simulated_data <- rnmix_skewed_with_outliers (200, theta =list(p=c(0.5, 0.5),
                                                               mu=c(0, 4),sigma=c(0.3, 0.3),
                                                               skew=c(0, 0)))


relevant_mixture_functions <- list ("em R" = list(name_fonction=emnmix, list_params=list()),
                                    "Rmixmod" = list(name_fonction=em_Rmixmod, list_params=list()),
                                    "mixtools" = list(name_fonction=em_mixtools, list_params=list()),
                                    "bgmm"= list(name_fonction=em_bgmm, list_params=list()),
                                    "mclust" = list(name_fonction=em_mclust, list_params=list(prior = NULL)),
                                    "EMCluster" = list(name_fonction=em_EMCluster, list_params=list()),
                                    "GMKMcharlie"=list(name_fonction=em_GMKMcharlie, list_params=list()),
                                    "flexmix"= list(name_fonction=em_flexmix, list_params=list()),
                                    "DCEM"=list(name_fonction=em_DCEM, list_params=list()))

parameters_distribution <- benchmark_distribution_parameters(mixture_functions=relevant_mixture_functions,
                                                                       initialisation_algorithms = c("kmeans", "quantiles"),
                                                                       sigma_values=list("small OVL"= c(0.3, 0.3)),
                                                                       mean_values=list(c(0, 4)),
                                                                       proportions = list( "balanced"=c(0.5, 0.5)),
                                                                       skewness_values = list("null skewness"=rep(0, 2)),prop_outliers = c(0),
                                                                       Nbootstrap=5,  nobservations=c(150))

running_times <- compute_microbenchmark(mixture_functions=relevant_mixture_functions, sigma_values=list("null OVL"=rep(0.3, 4)),
                                        mean_values=list(c(0, 4)), proportions = list("balanced"=c(0.5, 0.5), "little unbalanced"=c(0.65, 0.35)),
                                        skewness_values = list("null skewness"=c(0, 0)),prop_outliers = c(0),
                                        Nbootstrap=5, nobservations=c(50, 100, 200), initialisation_algorithms = c("quantiles"))
