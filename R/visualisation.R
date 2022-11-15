#' Plot the boxplot representation of the estimated parameters
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param true_theta the true values of the parameters (with p, mu and sigma)
#' @param remove_outliers boolean: remove or not outlying estimates, by default set to True
#' @param size_tag in each facet, controls the size of the parameter emphasised
#' @param num_col How to organise the facets? by default according to the number of components
#'
#' @return boxplot_parameters a ggplot object representing the boxplot distributions of the estimates per package and initialisation package
#'
#' @export

plot_boxplots_parameters <- function(distribution_parameters, true_theta, remove_outliers = T, size_tag=4, num_col=length(true_theta$p)) {

  # format true theta values, according to their use
  formatted_true_theta <- format_theta_output(true_theta)
  true_theta_df <- tibble::tibble(name_parameter = names(unlist(formatted_true_theta)) %>%
                                    factor(levels=unique(names(unlist(formatted_true_theta)))),
                                  true_value = unlist(formatted_true_theta))


  # format distribution data
  distribution_parameters <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sigma|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    select(c("package", "initialisation_method", "name_parameter", "value_parameter")) %>%
    dplyr::mutate(package=factor(package, levels = unique(package)),
                  name_parameter=factor(name_parameter, levels=unique(name_parameter)))

  if (remove_outliers) {
    distribution_parameters <- distribution_parameters %>%
      dplyr::group_by(name_parameter) %>%
      dplyr::filter(value_parameter > (quantile(value_parameter, probs = c(0.25)) - 1.5 * IQR(value_parameter)) &
                      value_parameter < (quantile(value_parameter, probs = c(0.75)) + 1.5 * IQR(value_parameter)))
  }
  # generate Boxplots
  boxplot_parameters <- ggplot(distribution_parameters, aes(x = package, y = value_parameter, fill = initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "yellow",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~name_parameter, ncol = num_col, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom",
      axis.ticks.length = unit(.1, "cm"),
      legend.text = element_text(size = 25),
      strip.text = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      title = element_blank(),
      panel.spacing = unit(.2, "pt")) +
    scale_fill_viridis_d() +
    geom_hline(data = true_theta_df, aes(yintercept = true_value), col = "red", linetype = "dashed", size = 0.8)

  if (rlang::is_installed("ggtext")) {
    boxplot_parameters <- boxplot_parameters + theme(axis.text.x = ggtext::element_markdown(angle = 90, size = 12, hjust=0.5, vjust = 0.4))}
  else {
    message("For the best formatting of packages labels, please install the ggtext if possible.")
    boxplot_parameters <- boxplot_parameters + theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.5, vjust = 0.4))
  }

  boxplot_parameters <- egg::tag_facet(boxplot_parameters,
                                       tag_pool = true_theta_df$name_parameter,
                                       open = "", close = "", hjust = -0.2, size = size_tag)

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

plot_ellipses_bivariate <- function(distribution_parameters, true_theta, alpha=0.05, npoints=500) {

  # generate ellipses with the true parameters
  k <- length(true_theta$p)
  true_theta_mean <- purrr::map_dfr(1:k, ~tibble::tibble(x=true_theta$mu[1,.x], y=true_theta$mu[2,.x], component=as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    true_theta_per_component <- list(p=true_theta$p[j], mu=true_theta$mu[,j], sigma=true_theta$sigma[,,j])
    ellipse_per_component <- generate_ellipse(true_theta_per_component, alpha = alpha, npoints = npoints) %>%
      tibble::add_column(component=as.character(j))
    return(ellipse_per_component)
  })


  ## generate data for ellipses with the estimated parameters
  # compute the average parameter
  distribution_parameters_mean <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sigma|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::mutate(across(c("name_parameter", "package"), as.factor)) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter")) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "name_parameter")))) %>%
    dplyr::summarise(mean_parameter=mean(value_parameter, na.rm = T)) %>% dplyr::ungroup()

  distribution_parameters_list <- purrr::map_dfr(1:k, function(j) {
    distribution_parameters_mean_per_component <- distribution_parameters_mean %>%
      dplyr::filter(stringr::str_detect(name_parameter, paste0(j, "$"))) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method")))) %>%
      dplyr::summarise(mean_parameter=list(stats::setNames(mean_parameter, nm = name_parameter) %>%
                                             as.list() %>% unformat_theta_output()),
                       component=j %>% as.character())
    return(distribution_parameters_mean_per_component)
  })
  # generate data for the ellipses
  ellipse_data <- distribution_parameters_list %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "component")))) %>%
    dplyr::summarize(generate_ellipse(mean_parameter %>% unlist(recursive = F), alpha = alpha, npoints = npoints))


  # create ellipse plot
  ellipse_plot <- ggplot(ellipse_data, aes(x=x, y=y)) +
    geom_path(aes(linetype=component, col=package), size=0.5) +
    theme_bw() +
    facet_wrap(~initialisation_method) +
    coord_fixed()+
    theme(legend.position = "bottom",
          title = element_blank(),
          legend.text = element_text(size = 18),
          axis.title = element_blank(),
          strip.text = element_text(size = 15, face = "bold")) +
    geom_path(data = ellipse_standard_data,mapping=aes(x=x, y=y, linetype=component),
              inherit.aes = F, col="black", size=2, alpha=0.25) +
    guides(linetype = guide_legend(override.aes = list(alpha = 0.25, shape=NA, size=1)))+
    scale_color_discrete(name=NULL) +
    ggnewscale::new_scale_color()+
    geom_point(data = true_theta_mean, mapping = aes(x=x, y=y, shape=component, col=component),
               size=3, inherit.aes = F, show.legend = F) +
    scale_color_manual(name="centroid", values=c("red", "green"),
                       guide = guide_legend(override.aes = list(shape = c(16, 17),
                                                                color = c("red", "green"))))



  return(ellipse_plot)
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
plot_Hellinger <- function(distribution_parameters, true_theta, num_col=length(true_theta$p)) {

  # format the data
  k <- length(true_theta$p)
  distribution_parameters_list <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sigma|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(package=factor(package, levels = unique(package))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "N.bootstrap"))),) %>%
    dplyr::summarise(list_parameter=list(stats::setNames(value_parameter, nm = name_parameter)%>%
                                           as.list() %>% unformat_theta_output()))

  if (is.matrix(true_theta$mu)) { # multivariate case
  # compute Hellinger distance per component
  hellinger_data <- purrr::map_dfr(1:k, function(j) {
    hellinger_data_per_component <- distribution_parameters_list %>%
      dplyr::mutate(mu=lapply(list_parameter, function(.x) .x$mu[,j]), component=as.character(j),
                    sigma=lapply(list_parameter, function(.x) .x$sigma[,,j]),
                    hellinger_value=purrr::map2_dbl(mu, sigma, hellinger, true_theta$mu[,j], true_theta$sigma[,,j])) %>%
      select(-c("list_parameter", "N.bootstrap", "mu", "sigma"))
    return(hellinger_data_per_component)
  })}
  else { # univariate case
    hellinger_data <- purrr::map_dfr(1:k, function(j) {
      hellinger_data_per_component <- distribution_parameters_list %>%
        dplyr::mutate(mu=lapply(list_parameter, function(.x) .x$mu[j]), component=as.character(j),
                      sigma=lapply(list_parameter, function(.x) .x$sigma[j]),
                      hellinger_value=purrr::map2_dbl(mu, sigma, hellinger, true_theta$mu[j], true_theta$sigma[j])) %>%
        select(-c("list_parameter", "N.bootstrap", "mu", "sigma"))
      return(hellinger_data_per_component)
    })
  }

  # generate Hellinger boxplot
  boxplot_Hellinger <- ggplot(hellinger_data, aes(x = factor(package, levels = unique(package)),
                                                  y = hellinger_value, fill = initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 1, colour = "red",
                 position = position_dodge(width = 0.9), show.legend = FALSE) +
    facet_wrap(~component, ncol = num_col) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 18),
          axis.ticks.length = unit(.1, "cm"),
          legend.text = element_text(size = 25),
          panel.spacing = unit(.2, "pt"),
          legend.title = element_blank(),
          axis.title.y = element_text(size=15),
          axis.title.x = element_blank(), plot.title = element_blank(), plot.subtitle = element_blank()) +
    scale_fill_viridis_d() +
    ylab("Hellinger distance")

  if (rlang::is_installed("ggtext"))
    boxplot_Hellinger <- boxplot_Hellinger + theme(axis.text.x = ggtext::element_markdown(angle = 90, size = 12, hjust=0.5, vjust = 0.4))
  else {
    message("For the best formatting of packages labels, please install the ggtext if possible.")
    boxplot_Hellinger <- boxplot_Hellinger + theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.5, vjust = 0.4))
  }
  return(boxplot_Hellinger)
}





#' Plot the running time of the reviewed packages
#'
#' @author Bastien CHASSAGNOL
#'
#' @param time_data the running time data associated to the EM estimation performed by the reviewed packages
#' @param grouping_colnames the names of the column, to establish the facets
#'
#' @return time_plots a ggplot object representing the running time curve distributions of the reviewed packages
#'
#' @export

plot_time_computations <- function(time_data){

  time_summary <- time_data %>% dplyr::group_by(across(all_of(c("package", "ID", "nobservations", "initialisation_method")))) %>%
    dplyr::summarise(time_median= median(time), time_up=quantile(time, probs = 0.95), time_down=quantile(time, probs =0.05)) %>%
    dplyr::mutate(time_median=log10(time_median), time_up=log10(time_up), time_down=log10(time_down), nobservations=log10(nobservations))

  time_plot <- ggplot(time_summary, aes(x = nobservations, y = time_median,
                                        col=package, linetype=package, shape=package)) +
    geom_point(size=3) +
    geom_line(size=1.15) +
    geom_ribbon(alpha=0.1, aes(ymin=time_down, ymax=time_up, fill=package), colour=NA)+
    #facet_grid(OVL ~ factor(entropy, levels=rev(unique(entropy))))+
    labs(x="Number of observations (log10)", y="Time in seconds (log10)") +
    theme_bw() +
    theme(legend.position="bottom", legend.title = element_blank(), axis.title = element_text(size=18),
          plot.title = element_blank(), legend.text = element_text(size=22)) +
    scale_shape_manual(values=c(1:length(unique(time_summary[["package"]]))))
  return(time_plot)
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

plot_univariate_normal_density_distribution <- function(true_theta, nobservations=1000, k=length(true_theta$p)) {
  # compute 0.05 and 0.95 quantiles of the univariate GMM
  min_comp <- which.min(true_theta$mu); max_comp <- which.max(true_theta$mu)
  min_value <- stats::qnorm(0.001, mean=true_theta$mu[min_comp], sd = true_theta$sigma[min_comp])
  max_value <- stats::qnorm(0.999, mean=true_theta$mu[max_comp], sd= true_theta$sigma[max_comp])


  # generate density_data, for each component and then by merging them altogether to reconstitute the overall distribution
  density_data <- tibble::tibble(x=seq(min_value, max_value, length.out = nobservations))
  for (i in 1:k) {
    density_data <- density_data %>% tibble::add_column("component {i}" := true_theta$p[i] *
                                                    stats::dnorm(density_data$x, mean=true_theta$mu[i], sd = true_theta$sigma[i]))
  }
  density_data <- density_data %>% dplyr::mutate(mixture=rowSums(dplyr::across(dplyr::starts_with("component ")))) %>%
    tidyr::pivot_longer(cols =dplyr::starts_with(c("component ", "mixture")), names_to = "components", values_to = "expression") %>%
    dplyr::mutate(components=as.factor(gsub("component ", "", components)))

  # generate the corresponding density plot
  density_plot <- ggplot(density_data, aes(x = x, color = components, y=expression,
                                        linetype=components, alpha=components, size = components)) +
    geom_line(key_glyph = draw_key_path) +
    theme_bw() +
    theme(plot.title = element_blank(), legend.position = "bottom", axis.title=element_blank(), legend.title = element_blank(),
          plot.subtitle = element_blank(), legend.text = element_text(size = 25), legend.key.width = unit(1.5, 'cm'),
          plot.tag = element_text(size=18, face = "bold")) +
    scale_linetype_manual(values = c(2:(k+1), 1)) +
    scale_alpha_manual(values = c(rep(1, k), 0.25)) +
    scale_size_manual(values = c(rep(1, k), 2))

  return (density_plot)
}

#' @rdname plot_univariate_normal_density_distribution
#' @export

plot_bivariate_normal_density_distribution <- function(true_theta, nobservations=1000, k=length(true_theta$p)) {

  # compute centroids and 95% confidence intervals
  theta_mean <- purrr::map_dfr(1:k, ~tibble::tibble(x=true_theta$mu[1,.x], y=true_theta$mu[2,.x], component=as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    theta_per_component <- list(p=true_theta$p[j], mu=true_theta$mu[,j], sigma=true_theta$sigma[,,j])
    ellipse_per_component <- generate_ellipse(theta_per_component, alpha = 0.05, npoints = nobservations) %>%
      tibble::add_column(component=as.character(j))
    return(ellipse_per_component)})

  # generate the two-d density distribution data
  density_data <- simulate_multivariate_GMM (true_theta, n=nobservations)$x
  colnames(density_data) <- c("x", "y"); density_data <- density_data %>% tibble::as_tibble()

  # plot the corresponding Ellipses
  isogradient <- ggplot(density_data, aes(x = x, y = y)) +
    stat_density_2d(geom = "tile", alpha=0.5,
                    aes(fill = ..density..),
                    contour = FALSE) +
    viridis::scale_fill_viridis() +
    theme_bw() +
    theme(axis.title.y = element_text(angle = 0, vjust=0.5), legend.key = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_point(mapping = aes(x=x, y=y, col=component, shape=component), data = theta_mean,
               inherit.aes = F, show.legend = T, size=5) +
    scale_color_manual(name="Cluster", values=c("red", "green"),
                       guide = guide_legend(override.aes = list(fill = NA,
                                                                shape = c(16, 17),
                                                                color = c("red", "green")))) +
    geom_path(data = ellipse_standard_data,mapping=aes(x=x, y=y, colour=component),
              inherit.aes = F, size=1.5, show.legend = F) +
    coord_fixed() +
    labs(fill="Density", colour="Cluster",shape="Cluster")

  return(isogradient)
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
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sigma|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(across(c("name_parameter", "package"), as.factor))

  # compute correlation matrix
  total_correlation_scores_global <- lapply(split(distribution_parameters_long, distribution_parameters_long$initialisation_method),
                                            function(x) {
                                              cor_data <- tidyr::pivot_wider(x, names_from = c("package"), values_from="value_parameter") %>%
                                                dplyr::select(unique(x %>% dplyr::pull(package))) %>% stats::cor(use="complete.obs")
                                              return(cor_data)
                                            })

  # generate associated Heatmap, for each initialisation package
  total_correlation_scores_plots <- Map(function(cor_matrix, init_method) {
    complex_heatmap <- ComplexHeatmap::Heatmap(cor_matrix, name = "mat", heatmap_legend_param = list(title = ""),
                                               cluster_rows = TRUE, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix),
                                               row_title_gp = grid::gpar(fontsize = 4), column_names_rot = 45,
                                               cluster_columns = TRUE, column_names_gp = grid::gpar(fontsize = 8), column_labels = colnames(cor_matrix),
                                               width = unit(8, "cm"), height = unit(8, "cm"), column_title = init_method)
    return(complex_heatmap)
    # return(grid::grid.grabExpr(ComplexHeatmap::draw(complex_heatmap)))
  }, total_correlation_scores_global, names(total_correlation_scores_global))

  return(total_correlation_scores_plots)
}



