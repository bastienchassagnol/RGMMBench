#' Plot the boxplot representation of the estimated parameters
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param true_theta the true values of the parameters (with p, mu and sigma)
#' @param remove_outliers boolean: remove or not outlying estimates, by default set to True
#' @param size_tag in each facet, controls the size of the parameter emphasised
#' @param num_col How to organise the facets? by default according to the number of components
#' @param match_symbol To discard the set of parameters you want to be removed
#'
#' @return boxplot_parameters a ggplot object representing the boxplot distributions of the estimates per package and initialisation package
#'
#' @export

plot_boxplots_parameters <- function(distribution_parameters, true_theta, remove_outliers = T,
                                     size_tag = 4, num_col = length(true_theta$p), match_symbol = "^p[[:digit:]]+|mu|sigma|sd") {

  # format true theta values, according to their use
  formatted_true_theta <- format_theta_output(true_theta)
  true_theta_df <- tibble::tibble(
    name_parameter = names(unlist(formatted_true_theta)) %>%
      factor(levels = unique(names(unlist(formatted_true_theta)))),
    true_value = unlist(formatted_true_theta)) %>%
    dplyr::filter(grepl(match_symbol, .data$name_parameter))


  # format distribution data
  distribution_parameters <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches(match_symbol),
      names_to = "name_parameter",
      values_to = "value_parameter"
    ) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter")) %>%
    dplyr::mutate(
      package = factor(.data$package, levels = unique(.data$package)),
      name_parameter = factor(.data$name_parameter, levels = unique(.data$name_parameter))
    )

  if (remove_outliers) {
    distribution_parameters <- distribution_parameters %>%
      dplyr::group_by(dplyr::across(c("name_parameter", "package"))) %>%
      dplyr::filter(.data$value_parameter > (quantile(.data$value_parameter, probs = c(0.25)) - 1.5 * stats::IQR(.data$value_parameter)) &
        .data$value_parameter < (quantile(.data$value_parameter, probs = c(0.75)) + 1.5 * stats::IQR(.data$value_parameter)))
  }
  # generate Boxplots
  boxplot_parameters <- ggplot(distribution_parameters, aes(x = .data$package, y = .data$value_parameter,
                                                            fill = .data$initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "yellow",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~ .data$name_parameter, ncol = num_col, scales = "free_y") +
    theme_bw() +
    theme(
      legend.position = "bottom",
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
    geom_hline(data = true_theta_df, aes(yintercept = .data$true_value), col = "red", linetype = "dashed", size = 0.8)

  if (rlang::is_installed("ggtext")) {
    boxplot_parameters <- boxplot_parameters + theme(axis.text.x = ggtext::element_markdown(angle = 90, size = 12, hjust = 0.5, vjust = 0.4))
  } else {
    message("For the best formatting of packages labels, please install the ggtext if possible.")
    boxplot_parameters <- boxplot_parameters + theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.4))
  }

  boxplot_parameters <- egg::tag_facet(boxplot_parameters,
    tag_pool = true_theta_df$name_parameter,
    open = "", close = "", hjust = -0.2, size = size_tag
  )

  return(boxplot_parameters)
}


#' Plot the confidence interval ellipses
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param true_theta the true values of the parameters (with p, mu and sigma)
#' @param alpha the confidence interval
#' @param npoints the number of points used to generate the ellipses
#'
#' @return a ggplot object representing for each configuration of parameter, in each panel the confidence
#' intervals associated to the median of the estimate of each package
#'
#' @export

plot_ellipses_bivariate <- function(distribution_parameters, true_theta, alpha = 0.05, npoints = 500) {

  # generate ellipses with the true parameters
  k <- length(true_theta$p)
  true_theta_mean <- purrr::map_dfr(1:k, ~ tibble::tibble(x = true_theta$mu[1, .x], y = true_theta$mu[2, .x], component = as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    true_theta_per_component <- list(p = true_theta$p[j], mu = true_theta$mu[, j], sigma = true_theta$sigma[, , j])
    ellipse_per_component <- generate_ellipse(true_theta_per_component, alpha = alpha, npoints = npoints) %>%
      tibble::add_column(component = as.character(j))
    return(ellipse_per_component)
  })


  ## generate data for ellipses with the estimated parameters
  # compute the average parameter
  distribution_parameters_mean <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("^p[[:digit:]]+|mu|sigma|sd"),
      names_to = "name_parameter",
      values_to = "value_parameter"
    ) %>%
    dplyr::mutate(dplyr::across(c("name_parameter", "package"), as.factor)) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter")) %>%
    dplyr::group_by(dplyr::across(c("package", "initialisation_method", "name_parameter"))) %>%
    dplyr::summarise(mean_parameter = mean(.data$value_parameter, na.rm = T)) %>%
    dplyr::ungroup()

  distribution_parameters_list <- purrr::map_dfr(1:k, function(j) {
    distribution_parameters_mean_per_component <- distribution_parameters_mean %>%
      dplyr::filter(stringr::str_detect(.data$name_parameter, paste0(j, "$"))) %>%
      dplyr::group_by(dplyr::across(c("package", "initialisation_method"))) %>%
      dplyr::summarise(
        mean_parameter = list(stats::setNames(.data$mean_parameter, nm = .data$name_parameter) %>%
          as.list() %>% unformat_theta_output()),
        component = j %>% as.character()
      )
    return(distribution_parameters_mean_per_component)
  })
  # generate data for the ellipses
  ellipse_data <- distribution_parameters_list %>%
    dplyr::group_by(dplyr::across(c("package", "initialisation_method", "component"))) %>%
    dplyr::summarize(generate_ellipse(.data$mean_parameter %>% unlist(recursive = F), alpha = alpha, npoints = npoints))


  # create ellipse plot
  ellipse_plot <- ggplot(ellipse_data, aes(x = .data$x, y = .data$y)) +
    geom_path(aes(linetype = .data$component, col = .data$package), size = 0.5) +
    theme_bw() +
    facet_wrap(~initialisation_method) +
    coord_fixed() +
    theme(
      legend.position = "bottom",
      title = element_blank(),
      legend.text = element_text(size = 18),
      axis.title = element_blank(),
      strip.text = element_text(size = 15, face = "bold")
    ) +
    geom_path(
      data = ellipse_standard_data, mapping = aes(x = .data$x, y = .data$y, linetype = .data$component),
      inherit.aes = F, col = "black", size = 2, alpha = 0.25
    ) +
    guides(linetype = guide_legend(override.aes = list(alpha = 0.25, shape = NA, size = 1))) +
    scale_color_discrete(name = NULL) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = true_theta_mean, mapping = aes(x = .data$x, y = .data$y, shape = .data$component, col = .data$component),
      size = 3, inherit.aes = F, show.legend = F
    ) +
    scale_color_manual(
      name = "centroid", values = c("red", "green"),
      guide = guide_legend(override.aes = list(
        shape = c(16, 17),
        color = c("red", "green")
      ))
    )



  return(ellipse_plot)
}

#' Generate the 14 configuration models
#'
#' @author Bastien CHASSAGNOL
#'
#' @param theta the true values of the parameters (with p, mu and sigma)
#'
#' @return a ggplot object representing for each configuration of parameter, in each panel the confidence
#' intervals associated to the median of the estimate of each package
#'
#' @export


plot_multi_parametrisation <- function(theta) {

  # retrieve the parameters of the multivariate distribution, and generate accordingly the ellipses
  k <- length(theta$p)
  theta_mean <- purrr::map_dfr(1:k, ~tibble::tibble(x=theta$mu[1,.x], y=theta$mu[2,.x], component=as.character(.x)))


  ellipse_eigen_decomposition <- purrr::map_dfr(1:k, function(j) {
    var_per_comp <- theta$sigma[,,j]
    eigen_values <- eigen(var_per_comp, symmetric = T)$values
    # from function DrawingEllipsesinR, to get the theta main angle of rotation
    # and https://cookierobotics.com/007/
    ellipse_equation <- tibble::tibble(a=sqrt(eigen_values[1]), b=sqrt(eigen_values[2]),
                                       angle=atan2(eigen_values[1] - var_per_comp[1, 1], var_per_comp[1, 2])) %>%
      dplyr::bind_cols(theta_mean %>% dplyr::filter(.data$component==as.character(j)))
  })

  ellipse_axis <- purrr::map_dfr(1:k, function(j) {
    ellipse_decomposition <- eigen(theta$sigma[,,j], symmetric = T)
    ellipse_coordinates_per_co <- matrix(t(ellipse_decomposition$vectors) * sqrt(ellipse_decomposition$values),
                                         nrow=2, dimnames = list(c("comp1", "comp2"), c("xend", "yend"))) %>%
      tibble::as_tibble() %>%
      dplyr::bind_cols(theta_mean %>% dplyr::filter(.data$component==as.character(j)))
    return(ellipse_coordinates_per_co)}) %>% # account for mean shift
    dplyr::mutate(xend=.data$x+.data$xend, yend=.data$y+.data$yend)


  isogradient <- ggplot(theta_mean) +
    theme_void() +
    theme(legend.position="none") +
    geom_segment(aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                 size = 2, arrow = arrow(length = unit(0.03, "inches")),
                 lineend = c('round'), linejoin = c('round'), linetype ="dotdash",
                 data=ellipse_axis, inherit.aes = F, alpha=0.8) +
    ggforce::geom_ellipse(aes(x0 = .data$x, y0 = .data$y, a = .data$a, b = .data$b, angle = .data$angle), size=2,
                          data=ellipse_eigen_decomposition, inherit.aes = F) +
    geom_point(aes(x=.data$x, y=.data$y), size=4, shape=4, stroke=2) +
    coord_fixed()


  return(isogradient)
}



#' Plot boxplot of the Hellinger distances per component
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param true_theta the true values of the parameters (with p, mu and sigma)
#' @param num_col How to organise the facets? by default according to the number of components
#'
#' @return a ggplot object representing for each configuration of parameter, in each panel the
#' Hellinger distance for each package and initialisation method, respectively to each component
#'
#' @export
plot_Hellinger <- function(distribution_parameters, true_theta, num_col = length(true_theta$p)) {

  # format the data
  k <- length(true_theta$p)
  distribution_parameters_list <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("^p[[:digit:]]+|mu|sigma|sd"),
      names_to = "name_parameter",
      values_to = "value_parameter"
    ) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(package = factor(.data$package, levels = unique(.data$package))) %>%
    dplyr::group_by(dplyr::across(c("package", "initialisation_method", "N.bootstrap"))) %>%
    dplyr::summarise(list_parameter = list(stats::setNames(.data$value_parameter, nm = .data$name_parameter) %>%
      as.list() %>% unformat_theta_output()))

  if (is.matrix(true_theta$mu)) { # multivariate case
    # compute Hellinger distance per component
    hellinger_data <- purrr::map_dfr(1:k, function(j) {
      hellinger_data_per_component <- distribution_parameters_list %>%
        dplyr::mutate(
          mu = lapply(.data$list_parameter, function(.x) .x$mu[, j]), component = as.character(j),
          sigma = lapply(.data$list_parameter, function(.x) .x$sigma[, , j]),
          hellinger_value = purrr::map2_dbl(.data$mu, .data$sigma, hellinger, true_theta$mu[, j], true_theta$sigma[, , j])
        ) %>%
        dplyr::select(-c("list_parameter", "N.bootstrap", "mu", "sigma"))
      return(hellinger_data_per_component)
    })
  } else { # univariate case
    hellinger_data <- purrr::map_dfr(1:k, function(j) {
      hellinger_data_per_component <- distribution_parameters_list %>%
        dplyr::mutate(
          mu = lapply(.data$list_parameter, function(.x) .x$mu[j]), component = as.character(j),
          sigma = lapply(.data$list_parameter, function(.x) .x$sigma[j]),
          hellinger_value = purrr::map2_dbl(.data$mu, .data$sigma, hellinger, true_theta$mu[j], true_theta$sigma[j])
        ) %>%
        dplyr::select(-c("list_parameter", "N.bootstrap", "mu", "sigma"))
      return(hellinger_data_per_component)
    })
  }

  # generate Hellinger boxplot
  boxplot_Hellinger <- ggplot(hellinger_data, aes(
    x = factor(.data$package, levels = unique(.data$package)),
    y = .data$hellinger_value, fill = .data$initialisation_method
  )) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "red",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~component, ncol = num_col) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 18),
      axis.ticks.length = unit(.1, "cm"),
      legend.text = element_text(size = 25),
      panel.spacing = unit(.2, "pt"),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_blank(), plot.title = element_blank(), plot.subtitle = element_blank()
    ) +
    scale_fill_viridis_d() +
    ylab("Hellinger distance")

  if (rlang::is_installed("ggtext")) {
    boxplot_Hellinger <- boxplot_Hellinger + theme(axis.text.x = ggtext::element_markdown(angle = 90, size = 12, hjust = 0.5, vjust = 0.4))
  } else {
    message("For the best formatting of packages labels, please install the ggtext if possible.")
    boxplot_Hellinger <- boxplot_Hellinger + theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 0.5, vjust = 0.4))
  }
  return(boxplot_Hellinger)
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
  time_summary <- time_data %>%
    dplyr::group_by(dplyr::across(c("package", "ID", "nobservations", "initialisation_method"))) %>%
    dplyr::summarise(time_median = median(.data$time), time_up = quantile(.data$time, probs = 0.95), time_down = quantile(.data$time, probs = 0.05)) %>%
    dplyr::mutate(time_median = log10(.data$time_median), time_up = log10(.data$time_up), time_down = log10(.data$time_down), nobservations = log10(.data$nobservations))

  time_plot <- ggplot(time_summary, aes(
    x = .data$nobservations, y = .data$time_median,
    col = .data$package, linetype = .data$package, shape = .data$package
  )) +
    geom_point(size = 3) +
    geom_line(size = 1.15) +
    geom_ribbon(alpha = 0.1, aes(ymin = .data$time_down, ymax = .data$time_up, fill = .data$package), colour = NA) +
    # facet_grid(OVL ~ factor(entropy, levels=rev(unique(entropy))))+
    labs(x = "Number of observations (log10)", y = "Time in seconds (log10)") +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.title = element_blank(), axis.title = element_text(size = 18),
      plot.title = element_blank(), legend.text = element_text(size = 22)
    ) +
    scale_shape_manual(values = c(1:length(unique(time_summary[["package"]]))))
  return(time_plot)
}




#' Plot the running time of the initialisation methods
#'
#' @author Bastien CHASSAGNOL
#'
#' @param init_time_data the running time data
#' @inheritDotParams egg::tag_facet x:family
#'
#' @return initialisation_time_plots a ggplot object representing the running time curve distributions of the initialisation packages
#'
#' @export

plot_initialisation_time_computations <- function(init_time_data, ...) {

  # get quantiles of the distribution
  init_time_data_summary <- init_time_data %>%
    dplyr::group_by(dplyr::across(c("ID", "nobservations", "initialisation_method"))) %>%
    dplyr::summarise(time_median = median(.data$time), time_up = quantile(.data$time, probs = 0.95), time_down = quantile(.data$time, probs = 0.05)) %>%
    dplyr::mutate(time_median = log10(.data$time_median), time_up = log10(.data$time_up), time_down = log10(.data$time_down), nobservations = log10(.data$nobservations))

  # plot quantiles of initialization time distribution
  initialisation_time_plots <- ggplot(init_time_data_summary, aes(
    x = .data$nobservations, y = .data$time_median, col = .data$initialisation_method,
    linetype = .data$initialisation_method, shape = .data$initialisation_method
  )) +
    geom_point(size = 3) +
    geom_line(aes(y = .data$time_median), size = 1.15) +
    geom_ribbon(alpha = 0.2, aes(ymin = .data$time_down, ymax = .data$time_up, fill = .data$initialisation_method), colour = NA) +
    facet_wrap(~ID) +
    theme_bw() +
    theme(
      legend.position = "bottom", legend.title = element_blank(), strip.text = element_blank(),
      plot.title = element_blank(), legend.text = element_text(size = 14), axis.title = element_text(size = 12)
    ) +
    # scale_shape_manual(values = c(1:length(unique(init_time_data_summary$initialisation_method)))) +
    labs(x = "Number of observations (log10)", y = "Time in seconds (log10)")


  initialisation_time_plots <- do.call(egg::tag_facet, list(p=initialisation_time_plots,
                                       tag_pool = init_time_data_summary$ID %>% unique(),
                                       open = "", close = "", ...))
  return(initialisation_time_plots)
}







#' Plot the true density distributions, in the univariate context
#'
#' @author Bastien CHASSAGNOL
#'
#' @param true_theta the true values of the parameters (with p, mu and sigma)
#' @param nobservations the number of bins used to regenerate the global mixture distribution
#' @param k the number of components, by default, the number of ratio elements stored in item p of true_theta
#'
#' @return density_plots a ggplot object representing the distributions of empirical univariate density functions
#'
#' @export

plot_univariate_normal_density_distribution <- function(true_theta, nobservations = 1000, k = length(true_theta$p)) {
  # compute 0.05 and 0.95 quantiles of the univariate GMM
  min_comp <- which.min(true_theta$mu)
  max_comp <- which.max(true_theta$mu)
  min_value <- stats::qnorm(0.001, mean = true_theta$mu[min_comp], sd = true_theta$sigma[min_comp])
  max_value <- stats::qnorm(0.999, mean = true_theta$mu[max_comp], sd = true_theta$sigma[max_comp])


  # generate density_data, for each component and then by merging them altogether to reconstitute the overall distribution
  density_data <- tibble::tibble(x = seq(min_value, max_value, length.out = nobservations))
  for (i in 1:k) {
    density_data <- density_data %>% tibble::add_column("component {i}" := true_theta$p[i] *
      stats::dnorm(density_data$x, mean = true_theta$mu[i], sd = true_theta$sigma[i]))
  }
  density_data <- density_data %>%
    dplyr::mutate(mixture = rowSums(dplyr::across(dplyr::starts_with("component ")))) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with(c("component ", "mixture")), names_to = "components", values_to = "expression") %>%
    dplyr::mutate(components = as.factor(gsub("component ", "", .data$components)))

  # generate the corresponding density plot
  density_plot <- ggplot(density_data, aes(
    x = .data$x, color = .data$components, y = expression,
    linetype = .data$components, alpha = .data$components, size = .data$components
  )) +
    geom_line(key_glyph = draw_key_path) +
    theme_bw() +
    theme(
      plot.title = element_blank(), legend.position = "bottom", axis.title = element_blank(), legend.title = element_blank(),
      plot.subtitle = element_blank(), legend.text = element_text(size = 25), legend.key.width = unit(1.5, "cm"),
      plot.tag = element_text(size = 18, face = "bold")
    ) +
    scale_linetype_manual(values = c(2:(k + 1), 1)) +
    scale_alpha_manual(values = c(rep(1, k), 0.25)) +
    scale_size_manual(values = c(rep(1, k), 2))

  return(density_plot)
}

#' @rdname plot_univariate_normal_density_distribution
#' @export

plot_bivariate_normal_density_distribution <- function(true_theta, nobservations = 1000, k = length(true_theta$p)) {

  # compute centroids and 95% confidence intervals
  theta_mean <- purrr::map_dfr(1:k, ~ tibble::tibble(x = true_theta$mu[1, .x], y = true_theta$mu[2, .x], component = as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    theta_per_component <- list(p = true_theta$p[j], mu = true_theta$mu[, j], sigma = true_theta$sigma[, , j])
    ellipse_per_component <- generate_ellipse(theta_per_component, alpha = 0.05, npoints = nobservations) %>%
      tibble::add_column(component = as.character(j))
    return(ellipse_per_component)
  })

  # generate the two-d density distribution data
  density_data <- simulate_multivariate_GMM(true_theta, n = nobservations)$x
  colnames(density_data) <- c("x", "y")
  density_data <- density_data %>% tibble::as_tibble()

  # plot the corresponding Ellipses
  isogradient <- ggplot(density_data, aes(x = .data$x, y = .data$y)) +
    stat_density_2d(
      geom = "raster", alpha = 0.5,
      aes(fill = after_stat(.data$density)),
      contour = FALSE
    ) +
    viridis::scale_fill_viridis() +
    theme_bw() +
    theme(
      axis.title.y = element_text(angle = 0, vjust = 0.5), legend.key = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    geom_point(
      mapping = aes(x = .data$x, y = .data$y, col = .data$component, shape = .data$component), data = theta_mean,
      inherit.aes = F, show.legend = T, size = 5
    ) +
    scale_color_manual(
      name = "Cluster", values = c("red", "green"),
      guide = guide_legend(override.aes = list(
        fill = NA,
        shape = c(16, 17),
        color = c("red", "green")
      ))
    ) +
    geom_path(
      data = ellipse_standard_data, mapping = aes(x = .data$x, y = .data$y, colour = .data$component),
      inherit.aes = F, size = 1.5, show.legend = F
    ) +
    coord_fixed() +
    labs(fill = "Density", colour = "Cluster", shape = "Cluster")

  return(isogradient)
}

#' @importFrom adegraphics s.class
#' @rdname plot_univariate_normal_density_distribution
#' @export
plot_HD_density_distribution <- function(true_theta, nobservations = 10^3,
                                         k = length(true_theta$p), ade_plot=FALSE) {

  tibble_dataset <- purrr::map_dfr(1:k, function(j) {
    x_per_component <- MASS::mvrnorm(n = nobservations, mu = true_theta$mu[,j],
                                     Sigma = true_theta$sigma[, , j], empirical = FALSE)
    tibble_per_component <- x_per_component %>% as_tibble() %>%
      tibble::add_column(component=paste("Comp", j) %>% as.factor())
    return(tibble_per_component)
  })


  pca1 <- ade4::dudi.pca(tibble_dataset %>% dplyr::select(-component),
                         scannf = FALSE, nf = 2, center=FALSE, scale=FALSE)
  mypalette <- rainbow(k); D <- nrow(true_theta$mu)
  if (requireNamespace("adegraphics", quietly = TRUE) & ade_plot) {
    eigen_plot <- adegraphics::s1d.barchart(pca1$eig, p1d.horizontal = F,
                                            ppolygons.col = "blue", plot = F)

    ind_plot <- s.class(pca1$li, fac = tibble_dataset$component,
                        plot=FALSE, col=mypalette,
                        pellipses.lwd = 2, pellipses.border = mypalette,
                        pellipses.col = mypalette,
                        starSize = 0, ppoints.cex = 0.2)
    var_plot <- adegraphics::s.corcircle(pca1$co,
                            lab = names(tibble_dataset %>% dplyr::select(-component)),
                            fullcircle = FALSE, plot = FALSE)


    ###  Attempt, failed, to convert the ade4 plot into a good looking ggplot

    # ind_plot_t <- ind_plot %>% adegraphics::gettrellis() %>% ggplotify::as.ggplot()
    # ind_plot_t + annotation_custom(grob=eigen_plot %>% adegraphics::gettrellis() %>% ggplotify::as.grob(),
    #                                xmin = 1, xmax = 3, ymin = -0.3, ymax = 0.6)
    # test <- gridExtra::arrangeGrob(grobs=list(ggplotify::as.ggplot(ind_plot %>% adegraphics::gettrellis())),
    #                                bottom=ggplotify::as.ggplot(ind_plot %>% adegraphics::gettrellis()),
    #                                left=ggplotify::as.ggplot(eigen_plot %>% adegraphics::gettrellis()))
    # ggsave("test_HD_adegraphics.pdf", ind_plot %>% adegraphics::gettrellis() %>% ggplotify::as.grob())


    HD_dp <- adegraphics::insert(eigen_plot, ind_plot,  ratio=0.25, plot=F,
                                 posi = "topleft",  inset = c(0, -0.05))
    HD_dp <- adegraphics::insert(var_plot, HD_dp,  ratio=0.3,
                                    posi = "bottomright", plot=F)
  }
  else { # with ggplot object
    eigen_plot <- factoextra::fviz_eig(pca1, main=NULL, ggtheme = theme_bw(), addlabels = F) +
      ggtitle(element_blank()) +
      ylab("") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "null"),
            panel.spacing = unit(c(0, 0, 0, 0), "null"))
    indiv_plot <- factoextra::fviz_pca_biplot(pca1, # Individuals
                    geom.ind = "point",
                    fill.ind = tibble_dataset$component,
                    pointshape = 21, pointsize = 1, ind.var=1,
                    repel = TRUE, addEllipses = TRUE, ellipse.level=0.95, # Variables
                    col.var = "contrib", gradient.cols = c("blue", "white", "red"), geom.var = "arrow",
                    ggtheme = theme_bw())+
      labs(fill = "Clusters", color = "Contribution", title = element_blank())+ # Change legend title
      theme(axis.title = element_text(size = 14, face = 'bold'),
            axis.text = element_text(size = 14, face = 'bold')) +
      # hrbrthemes::theme_ipsum()+
      coord_fixed()

    HD_dp <- cowplot::plot_grid(eigen_plot, indiv_plot,nrow=2, align = "hv", axis = "tblr")

    # modified_eigen_plot <- eigen_plot + # set the background as transparent
    #   xlab("") + ylab("") + hrbrthemes::theme_ipsum() +
    #   theme(panel.background = element_rect(fill=alpha("white", 0.2)),
    #         plot.background =  element_rect(fill = "transparent"),
    #         plot.margin = margin(t = 0, r = 0,   b = 0,  l = 0),
    #         axis.title=element_blank(),
    #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #         axis.text.x=element_blank(), axis.text.y=element_blank(),
    #         axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
    #
    # limits_ggplot <- get_plot_limits(indiv_plot)
    # range_y <- limits_ggplot$ymax - limits_ggplot$ymin; range_x <- limits_ggplot$xmax - limits_ggplot$xmin
    # HD_dp <- indiv_plot +
    #   annotation_custom(grob = ggplotify::as.grob(modified_eigen_plot),
    #                     xmin = limits_ggplot$xmin + 0.6 * range_x, xmax = limits_ggplot$xmax,
    #                     ymin = limits_ggplot$ymin + 0.6 * range_y, ymax = limits_ggplot$ymax)
  }

  # generate parallel distirbution plots
  sampled_tibble_dataset <- tibble_dataset %>% dplyr::slice_sample(n=50)
  parallel_plot_unscaled <- GGally::ggparcoord(sampled_tibble_dataset, showPoints = T,
             columns = 1:D, groupColumn = D + 1, scale="globalminmax", alphaLines = 0.3) +
    scale_color_discrete(type=mypalette) +
    # hrbrthemes::theme_ipsum()+
    theme_bw() +
    labs(title = "Parallel Coordinate Plot",
         subtitle = "No scaling") +
    xlab("") + ylab("")


  parallel_plot_scaled <- GGally::ggparcoord(sampled_tibble_dataset, showPoints = T,
                                             columns = 1:D, groupColumn = D + 1,
                                             scale="std", alphaLines = 0.3) +
    scale_color_discrete(type=mypalette) +
    theme(title=element_blank()) +
    labs(subtitle = "Univariately normalised") +
    theme_bw() +
    # hrbrthemes::theme_ipsum()
    xlab("") + ylab("")
  parallel_plot <- cowplot::plot_grid(parallel_plot_unscaled, parallel_plot_scaled,
                                      align="hv", axis="tblr",nrow=2)


  # final concatenation of plot
  final_dp <- cowplot::plot_grid(HD_dp + theme(plot.tag = element_text(size=24, vjust = 1.5, hjust=-0.5, face = "bold")) +
                                   labs(tag = "A"),
                                 parallel_plot+ theme(plot.tag = element_text(size=24, face = "bold",
                                                                              vjust = 1.5, hjust=-0.5)) + labs(tag = "B"),
                                 ncol=2, align = "hv", axis="tblr")

  # final_dp <- cowplot::plot_grid(HD_dp + theme(plot.tag = element_text(size=24, face = "bold")) + labs(tag = "A") +
  #                                  labs(title = "Individual and variable projection"),
  #                                parallel_plot_unscaled + theme(plot.tag = element_text(size=24, face = "bold")) + labs(tag = "B"),
  #                                ncol=2, align = "v", axis="tblr", rel_widths = c(1, 2))

  return(final_dp)
}







#' Plot Correlation Heatmap
#'
#' For each configuration of parameter, compare the similarity between the packages
#' for each initialisation package, by computing the correlation between
#' their estimates
#'
#' @param distribution_parameters A tibble with the distribution of the parameters,
#' in the long format
#'
#' @export

plot_correlation_Heatmap <- function(distribution_parameters) {
  # format the data
  distribution_parameters_long <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("^p[[:digit:]]+|mu|sigma|sd"),
      names_to = "name_parameter",
      values_to = "value_parameter"
    ) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(dplyr::across(c("name_parameter", "package"), as.factor))

  # compute correlation matrix
  total_correlation_scores_global <- lapply(
    split(distribution_parameters_long, distribution_parameters_long$initialisation_method),
    function(x) {
      cor_data <- tidyr::pivot_wider(x, names_from = c("package"), values_from = "value_parameter") %>%
        dplyr::select(unique(x %>% dplyr::pull("package"))) %>%
        stats::cor(use = "complete.obs")
      return(cor_data)
    }
  )

  # generate associated Heatmap, for each initialisation package
  total_correlation_scores_plots <- Map(function(cor_matrix, init_method) {
    complex_heatmap <- ComplexHeatmap::Heatmap(cor_matrix,
      name = "mat", heatmap_legend_param = list(title = ""),
      cluster_rows = TRUE, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix),
      row_title_gp = grid::gpar(fontsize = 4), column_names_rot = 45,
      cluster_columns = TRUE, column_names_gp = grid::gpar(fontsize = 8), column_labels = colnames(cor_matrix),
      width = unit(8, "cm"), height = unit(8, "cm"), column_title = init_method
    )
    return(complex_heatmap)
    # return(grid::grid.grabExpr(ComplexHeatmap::draw(complex_heatmap)))
  }, total_correlation_scores_global, names(total_correlation_scores_global))

  return(total_correlation_scores_plots)
}
