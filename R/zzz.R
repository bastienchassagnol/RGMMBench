
# em_mixsmsn <- function(x = x, k = 2, initialisation_algorithm = "hc", skew = rep(0, k),
#                        itmax = 5000, epsilon = 10^-12, start = NULL, ...) {
#   if (is.null(start)) {
#     start <- initialize_em_univariate(
#       x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
#       initialisation_algorithm = initialisation_algorithm, ...
#     )
#   }
#
#   # nu is the kurtosis, always equal to 3 for a classic Gaussian distribution(not skewed)
#   if (all(skew == 0)) {
#     # fit a classic Gaussian distribution
#     fit <- mixsmsn::smsn.mix(
#       y = x, nu = 3,
#       mu = start$mu, sigma2 = start$sigma^2, shape = skew, pii = start$p,
#       g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Normal",
#       error = 10^-12, iter.max = 5000, calc.im = FALSE
#     )
#   } else {
#     # fit a skewed distribution
#     fit <- mixsmsn::smsn.mix(
#       y = x, nu = 3,
#       mu = start$mu, sigma2 = start$sigma, shape = skew, pii = start$p,
#       g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Skew.normal",
#       error = 10^-12, iter.max = 5000, calc.im = FALSE
#     )
#   }
#
#
#   # return an ordered list by mean values
#   ordered_estimated_theta <- list(
#     p = fit$pii[order(fit$mu)],
#     mu = sort(fit$mu),
#     sigma = sqrt(fit$sigma2[order(fit$mu)])
#   )
#
#   ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
#   return(ordered_estimated_theta)
# }



# if (with_outliers)
#   time_plot <- time_plot +
#   facet_grid(~ factor(prop_outliers, levels=unique(prop_outliers), labels=unique(paste("Proportion of outliers:", unique(prop_outliers)))))
