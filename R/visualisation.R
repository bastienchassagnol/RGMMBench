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


#' Plot the true density distributions
#'
#' @author Bastien CHASSAGNOL
#'
#' @param proportions,mean_values,sigma_values,skewness_values parameters of the several components represented,
#' represented by a list of vector parameters
#'
#' @return density_plots a ggplot object representing the distributions of the density functions
#'
#' @export

plot_density_distribution <- function(proportions, mean_values,
                                      sigma_values, skewness_values) {
  density_plots <- list()
  iteration <- 1
  #################################################################
  ##               generation of the distributions               ##
  #################################################################
  for (skew in skewness_values) {
    for (p in proportions) {
      for (mu in mean_values) {
        for (sigma in sigma_values) {
          true_theta <- list(p = p, mu = mu, sigma = sigma, skew = skew) # true parameters of the distribution
          k <- length(p) # number of components
          if (any(sapply(true_theta, length) != k)) {
            stop("One of the original componenent is not completely provided")
          }

          # compute 0.001 and 0.999 quantiles of min and max components
          min_value <- sapply(1:k, function(i) {
            sn::qsn(0.001, xi = true_theta$mu[i], omega = true_theta$sigma[i], alpha = true_theta$skew[i], tau = 0)
          }) %>% min()
          max_value <- sapply(1:k, function(i) {
            sn::qsn(0.999, xi = true_theta$mu[i], omega = true_theta$sigma[i], alpha = true_theta$skew[i], tau = 0)
          }) %>% max()

          distribution_data <- tibble::tibble(x = seq(min_value, max_value, length.out = 10000))
          # add a column for each component
          for (i in 1:k) {
            distribution_data <- distribution_data %>% tibble::add_column("component {i}" := true_theta$p[i] *
              sn::dsn(distribution_data$x, xi = true_theta$mu[i], omega = true_theta$sigma[i], alpha = true_theta$skew[i], tau = 0))
          }

          distribution_data <- distribution_data %>%
            dplyr::mutate(total = rowSums(dplyr::across(dplyr::starts_with("component ")))) %>%
            tidyr::pivot_longer(cols = dplyr::starts_with(c("component ", "total")), names_to = "components", values_to = "expression") %>%
            dplyr::mutate(skew = skew[1], components = as.factor(gsub("total", "mixture\ndistribution", gsub("component ", "", components))))


          ##################################################################
          ##               plot the simulated distributions               ##
          ##################################################################
          plot_title <- ""
          for (i in 1:k) {
            plot_title <- paste0(
              plot_title, "Component ", i, " has parameters ",
              "p", i, ": ", signif(true_theta$p[i], digits = 2), ", ",
              "mu", i, ": ", true_theta$mu[i], ", ",
              "sigma", i, ": ", true_theta$sigma[i], ", ",
              "skew", i, ": ", true_theta$skew[i], ".\n"
            )
          }

          density_plots[[paste("Plot ", iteration)]] <- ggplot(distribution_data, aes(
            x = x, color = components, y = expression,
            linetype = components, alpha = components, size = components
          )) +
            geom_line(key_glyph = draw_key_path) +
            theme_bw() +
            theme(
              plot.title = element_blank(), legend.position = "bottom", axis.title = element_blank(), plot.subtitle = element_text(hjust = 0.5),
              legend.title = element_blank(), legend.text = element_text(size = 25), legend.key.width = unit(1.5, "cm")
            ) +
            labs(subtitle = plot_title) +
            scale_linetype_manual(values = c(2:(k + 1), 1)) +
            scale_alpha_manual(values = c(rep(1, k), 0.25)) +
            scale_size_manual(values = c(rep(1, k), 2))
          iteration <- iteration + 1
        }
      }
    }
  }
  return(density_plots)
}
