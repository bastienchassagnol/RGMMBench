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


plot_boxplots_parameters <- function(distribution_parameters, p, mu, sigma,
                                     remove_outliers = T, size_tag=8) {

  # format true theta values, according to their use
  true_theta <- list(p = p, mu = mu, sigma = sigma) # true parameters of the distribution
  true_theta_df <- tibble::tibble(
    name_parameter = names(unlist(true_theta[c("p", "mu", "sigma")])),
    true_value = unlist(true_theta[c("p", "mu", "sigma")])
  ) %>%
    dplyr::mutate(name_parameter = factor(name_parameter, levels = unique(name_parameter)))
  k <- length(p) # number of components

  if (remove_outliers) {
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


#' Plot the true density distributions, in the univariate context
#'
#' @author Bastien CHASSAGNOL
#'
#' @param proportions,mean_values,sigma_values parameters of the several components represented,
#' represented by a list of vector parameters
#'
#' @return density_plots a ggplot object representing the distributions of the density functions
#'
#' @export

plot_univariate_normal_density_distribution <- function(sigma_values, mean_values, proportions) {

  individual_plots <- list() # store for each plot its corresponding configuration
  for(p in proportions) {
    for (sigma in sigma_values) {
      for (mu in mean_values) {
        # format true theta values, according to their use
        true_theta <- list(p=p, mu=mu, sigma=sigma) # true parameters of the distribution
        k <- length(p) # number of components
        if (any(sapply(true_theta, length)!=k)) {
          stop("One of the original componenent is not completely provided")
        }
        # compute average OVL, average pairwise overlap including proportions and UR for each set of parameter
        balanced_ovl <- MixSim::overlap(Pi=rep(1/k, k), Mu=as.matrix(true_theta$mu), S=as.matrix(true_theta$sigma))$BarOmega
        entropy_value <- compute_shannon_entropy(p)
        # pairwise_ovl <- MixSim::overlap(Pi=true_theta$p, Mu=as.matrix(true_theta$mu), S=as.matrix(true_theta$sigma))$BarOmega

        filename <- paste0("entropy_", signif(entropy_value, digits = 2),"_OVL_", signif(balanced_ovl, digits = 2))

        # compute quantiles of min and max components
        min_comp <- which.min(true_theta$mu); max_comp <- which.max(true_theta$mu)
        min_value <- stats::qnorm(0.01, mean=true_theta$mu[min_comp], sd = true_theta$sigma[min_comp])
        max_value <- stats::qnorm(0.99, mean=true_theta$mu[max_comp], sd= true_theta$sigma[max_comp])
        # min_value <- -2; max_value <- 15

        temp_simu <- tibble::tibble(x=seq(min_value, max_value, length.out = 500))
        # add a column for each component
        for (i in 1:k) {
          temp_simu <- temp_simu %>% tibble::add_column("component {i}" := true_theta$p[i] *
                                                          stats::dnorm(temp_simu$x, mean=true_theta$mu[i], sd = true_theta$sigma[i]))
        }

        temp_simu <- temp_simu %>% dplyr::mutate(total=rowSums(dplyr::across(dplyr::starts_with("component ")))) %>%
          tidyr::pivot_longer(cols =dplyr::starts_with(c("component ", "total")), names_to = "components", values_to = "expression") %>%
          dplyr::mutate(entropy=signif(entropy_value, digits = 2), OVL=signif(balanced_ovl, digits = 2),
                        components=as.factor(gsub("total", "mixture\ndistribution",gsub("component ", "", components))))


        individual_plots[[filename]] <- ggplot(temp_simu, aes(x = x, color = components, y=expression,
                                                              linetype=components, alpha=components, size = components)) +
          geom_line(key_glyph = draw_key_path) +
          theme_bw() +
          theme(plot.title = element_blank(), legend.position = "bottom", axis.title=element_blank(),
                plot.subtitle = element_blank(), legend.title = element_blank(),
                legend.text = element_text(size = 25), legend.key.width = unit(1.5, 'cm'),
                plot.tag = element_text(size=18, face = "bold")) +
          # labs(subtitle = title_dp) +
          scale_linetype_manual(values = c(2:(k+1), 1)) +
          scale_alpha_manual(values = c(rep(1, k), 0.25)) +
          scale_size_manual(values = c(rep(1, k), 2))
      }
    }
  }
  return(individual_plots)
}

#' Plot Correlation Heatmap
#'
#' For each configuration of parameter, compare the similarity between the algorithms
#' for each initialisation algorithm, by computing the correlation between
#' their estimates
#'
#' @param distribution_parameters A tibble with the distribution of the parameters,
#' in the long format
#'
#' @export

plot_correlation_Heatmap <- function(distribution_parameters) {
  distribution_parameters <- distribution_parameters %>% dplyr::group_by(OVL, entropy, algorithm, name_parameter, initialisation_method) %>%
    mutate(index=dplyr::row_number()) %>% dplyr::ungroup() # add an index, ti uniquely identify each experiment


  # compute correlation matrix
  total_correlation_scores_global <- lapply(split(distribution_parameters, distribution_parameters$initialisation_method),
                                            function(x) {
                                              cor_data <- tidyr::pivot_wider(x, names_from = c("algorithm"), values_from="value_parameter") %>%
                                                dplyr::select(unique(x %>% dplyr::pull(algorithm))) %>% stats::cor(use="complete.obs")
                                              return(cor_data)
                                            })

  # generate associated Heatmap, for each initialisation algorithm
  total_correlation_scores_plots <- Map(function(cor_matrix, init_method) {
    complex_heatmap <- ComplexHeatmap::Heatmap(cor_matrix, name = "mat", heatmap_legend_param = list(title = ""),
                                               cluster_rows = TRUE, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix),
                                               row_title_gp = grid::gpar(fontsize = 4), column_names_rot = 45,
                                               cluster_columns = TRUE, column_names_gp = grid::gpar(fontsize = 8), column_labels = colnames(cor_matrix),
                                               width = unit(8, "cm"), height = unit(8, "cm"),
                                               column_title = paste("Initialization method:", init_method),
                                               column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"))
    return(grid::grid.grabExpr(ComplexHeatmap::draw(complex_heatmap)))
  }, total_correlation_scores_global, names(total_correlation_scores_global))

  return(total_correlation_scores_plots)
}
