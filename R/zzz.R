#' Plot the running time of the initialisation methods
#'
#' @author Bastien CHASSAGNOL
#'
#' @param init_time_data the running time data
#'
#' @return initialisation_time_plots a ggplot object representing the running time curve distributions of the initialisation packages
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


# if (with_outliers)
#   time_plot <- time_plot +
#   facet_grid(~ factor(prop_outliers, levels=unique(prop_outliers), labels=unique(paste("Proportion of outliers:", unique(prop_outliers)))))
