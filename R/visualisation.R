#' Plot the boxplot representation of the estimated parameters
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param p,mu,sigma the true values of the parameters
#' @param with_outliers boolean: remove or not outlying estimates
#'
#' @return boxplot_parameters a ggplot object representing the boxplot distributions of the estimates per package and initialisation algorithm
#'
#' @export


plot_boxplots_parameters <- function(distribution_parameters, p, mu, sigma, with_outliers = FALSE) {

  # format true theta values, according to their use
  true_theta <- list(p = p, mu = mu, sigma = sigma) # true parameters of the distribution
  true_theta_df <- tibble::tibble(
    name_parameter = names(unlist(true_theta[c("p", "mu", "sigma")])),
    true_value = unlist(true_theta[c("p", "mu", "sigma")])
  ) %>%
    dplyr::mutate(name_parameter = factor(name_parameter, levels = unique(name_parameter)))
  k <- length(p) # number of components

  if (!with_outliers) {
    distribution_parameters <- distribution_parameters %>%
      dplyr::group_by(name_parameter) %>%
      dplyr::filter(value_parameter > (quantile(value_parameter, probs = c(0.25)) - 1.5 * IQR(value_parameter)) &
        value_parameter < (quantile(value_parameter, probs = c(0.75)) + 1.5 * IQR(value_parameter)))
  }

  # display parameters, removing outliers for visualisation purposes
  boxplot_parameters <- ggplot(distribution_parameters, aes(x = package, y = value_parameter, fill = initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "yellow",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~name_parameter, ncol = k, scales = "free_y") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, size = 18, vjust = 0.6),
      axis.ticks.length = unit(.1, "cm"),
      legend.text = element_text(size = 25),
      strip.text = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      title = element_blank(),
      panel.spacing = unit(.2, "pt")
    ) +
    scale_fill_viridis_d() +
    geom_hline(data = true_theta_df, aes(yintercept = true_value), col = "red", linetype = "dashed", size = 0.8)

  boxplot_parameters <- egg::tag_facet(boxplot_parameters,
    tag_pool = unique(distribution_parameters$name_parameter),
    open = "", close = "", hjust = -0.2, size = 5
  )
  return(boxplot_parameters)
}


#' Plot the running time of the initialisation methods
#'
#' @author Bastien CHASSAGNOL
#'
#' @param init_time_data the running time data
#'
#' @return initialisation_time_plots a ggplot object representing the running time curve distributions of the initialisation algorithms
#'
#' @export

plot_initialisation_time_computations <- function(init_time_data) {

  # get quantiles of the distribution
  init_time_data_summary <- init_time_data %>%
    dplyr::arrange(dplyr::desc(entropy), OVL) %>%
    dplyr::mutate(
      OVL = factor(OVL, labels = paste("Balanced OVL:", unique(OVL)), levels = unique(sort(OVL))),
      entropy = factor(entropy, labels = paste("Entropy:", unique(entropy)), levels = unique(sort(entropy, decreasing = TRUE)))
    ) %>%
    dplyr::group_by(
      OVL, entropy, nobservations, prop_outliers,
      skew, initialisation_method
    ) %>%
    dplyr::summarise(time_median = log10(median(time)), time_up = log10(quantile(time, probs = 0.95)), time_down = log10(quantile(time, probs = 0.05))) # use of log10 scale for computation purposes

  # plot quantiles of initialization time distribution
  initialisation_time_plots <- ggplot(init_time_data_summary, aes(
    x = nobservations, y = time_median, col = initialisation_method,
    linetype = initialisation_method, shape = initialisation_method
  )) +
    geom_point(size = 3) +
    geom_line(aes(y = time_median), size = 1.15) +
    geom_ribbon(alpha = 0.2, aes(ymin = time_down, ymax = time_up, fill = initialisation_method), colour = NA) +
    facet_grid(OVL ~ factor(entropy, levels = rev(unique(entropy)))) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.title = element_blank(),
      plot.title = element_blank(), legend.text = element_text(size = 14), axis.title = element_text(size = 12)
    ) +
    scale_shape_manual(values = c(1:length(unique(init_time_data_summary$initialisation_method)))) +
    labs(x = "Number of observations (log10)", y = "Time in seconds (log10)")
  return(initialisation_time_plots)
}

#' Plot the running time of the reviewed packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param time_data the running time data associated to the EM estimation performed by the reviewed packages
#'
#' @return time_plots a ggplot object representing the running time curve distributions of the reviewed packages
#'
#' @export

plot_time_computations <- function(time_data) {
  time_data_summary <- time_data %>%
    dplyr::group_by(
      package, OVL, entropy, nobservations, prop_outliers,
      skew, initialisation_method
    ) %>%
    dplyr::summarise(time_median = log10(median(time)), time_up = log10(quantile(time, probs = 0.95)), time_down = log10(quantile(time, probs = 0.05))) %>%
    # use of log for visualisation purposes
    dplyr::arrange(dplyr::desc(entropy), OVL) %>%
    dplyr::mutate(
      OVL = factor(OVL, labels = paste("Balanced OVL:", unique(OVL)), levels = unique(sort(OVL))),
      entropy = factor(entropy, labels = paste("Entropy:", unique(entropy)), levels = unique(sort(entropy, decreasing = TRUE)))
    )

  splitted_time_data <- split(time_data_summary, time_data_summary %>% dplyr::pull(initialisation_method))
  time_plots <- Map(function(data_per_algo, name_algo) {
    plot_per_algo <- ggplot(data_per_algo, aes(
      x = nobservations, y = time_median,
      col = package, linetype = package, shape = package
    )) +
      geom_point(size = 3) +
      geom_line(size = 1.15) +
      geom_ribbon(alpha = 0.1, aes(ymin = time_down, ymax = time_up, fill = package), colour = NA) +
      facet_grid(OVL ~ factor(entropy, levels = rev(unique(entropy)))) +
      # facet_grid(~ factor(prop_outliers, levels=unique(prop_outliers), labels=unique(paste("Proportion of outliers:", unique(prop_outliers))))) +
      labs(x = "Number of observations (log10)", y = "Time in seconds (log10)") +
      theme_bw() +
      theme(
        legend.position = "bottom", legend.title = element_blank(), axis.title = element_text(size = 18),
        plot.title = element_blank(), legend.text = element_text(size = 25)
      ) +
      scale_shape_manual(values = c(1:length(unique(data_per_algo[["package"]]))))
    return(plot_per_algo)
  }, splitted_time_data, names(splitted_time_data))

  return(time_plots)
}
