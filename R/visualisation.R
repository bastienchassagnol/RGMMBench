#' Plot the boxplot representation of the estimated parameters
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param p,mu,sigma the true values of the parameters
#' @param with_outliers boolean: remove or not outlying estimates
#'
#' @return boxplot_parameters a ggplot object representing the boxplot distributions of the estimates per package and initialisation package
#'
#' @export


plot_boxplots_parameters_univariate <- function(distribution_parameters, p, mu, sigma,
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



#' @rdname plot_boxplots_parameters_univariate
#' @export

plot_boxplots_parameters_bivariate <- function(distribution_parameters, remove_outliers = T,
                                               size_tag=8, num_col=NULL) {

  # format true theta values, according to their use
  true_theta <- distribution_parameters$formatted_true_parameters[[1]]
  if(missing(num_col)) num_col <- length(distribution_parameters$true_parameters[[1]]$p)
  true_theta_df <- tibble::tibble(name_parameter = names(unlist(true_theta)), true_value = unlist(true_theta))


  # format the data
  distribution_parameters <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::mutate(across(c("name_parameter", "package"), as.factor)) %>%
    select(c("package", "initialisation_method", "name_parameter", "value_parameter"))

  if (remove_outliers) {
    distribution_parameters <- distribution_parameters %>%
      dplyr::group_by(name_parameter) %>%
      dplyr::filter(value_parameter > (quantile(value_parameter, probs = c(0.25)) - 1.5 * IQR(value_parameter)) &
                      value_parameter < (quantile(value_parameter, probs = c(0.75)) + 1.5 * IQR(value_parameter)))
  }

  # my_packages <- unique(distribution_parameters %>% pull(package))
  # my_labels <- stats::setNames(paste0("style='color:",
  # dplyr::if_else(my_packages %in% c("em R", "Rmixmod", "mixtools"), 'red', 'green'), "'> ", my_packages), nm=my_packages)
  # my_labels <- purrr::imap_chr(my_labels, ~ ifelse(.y %in% c("mclust", "flexmix", "mixtools"), paste0("<b ", .x, "</b>"), .x))

  # display parameters, removing outliers for visualisation purposes
  boxplot_parameters <- ggplot(distribution_parameters, aes(x = factor(package, levels = unique(package)),
                                                            y = value_parameter, fill = initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "yellow",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~name_parameter, ncol = num_col, scales = "free_y") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, size = 15, vjust = 0.6),
      # axis.text.x = ggtext::element_markdown(angle = 35, size = 15, vjust =0.6),
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

  boxplot_parameters <- egg::tag_facet(boxplot_parameters, tag_pool = distribution_parameters$name_parameter %>% unique() %>% sort(),
                                       open = "", close = "", hjust = -0.2, size = 4)

  return(boxplot_parameters)
}


#' Plot the confidence interval ellipses
#'
#' @author Bastien CHASSAGNOL
#'
#' @param distribution_parameters the estimated bootstrap distributions
#' @param alpha the confidence interval
#' @param npoints the number of points used to generate the ellipses
#'
#' @return a ggplot object representing for each configuration of parameter, in each panel the confidence
#' intervals associated to the median of the estimate of each package
#'
#' @export

plot_ellipses_bivariate <- function(distribution_parameters, alpha=0.05, npoints=500) {

  # format true theta values, according to their use
  true_theta <- distribution_parameters$true_parameters[[1]]; k <- length(true_theta$p)
  true_theta_mean <- purrr::map_dfr(1:k, ~tibble::tibble(x=true_theta$mu[1,.x], y=true_theta$mu[2,.x], component=as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    true_theta_per_component <- list(p=true_theta$p[j], mu=true_theta$mu[,j], sigma=true_theta$sigma[,,j])
    ellipse_per_component <- generate_ellipse(true_theta_per_component, alpha = alpha, npoints = npoints) %>%
      tibble::add_column(component=as.character(j))
    return(ellipse_per_component)
  })


  # format the data
  distribution_parameters_mean <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::mutate(across(c("name_parameter", "package"), as.factor)) %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter")) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "name_parameter")))) %>%
    dplyr::summarise(mean_parameter=mean(value_parameter, na.rm = T)) %>%
    dplyr::ungroup()

  # format parameters as an extended list
  distribution_parameters_list <- purrr::map_dfr(1:k, function(j) {
    distribution_parameters_mean_per_component <- distribution_parameters_mean %>%
      dplyr::filter(stringr::str_detect(name_parameter, paste0(j, "$"))) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method")))) %>%
      dplyr::summarise(mean_parameter=list(stats::setNames(mean_parameter, nm = name_parameter) %>%
                                                as.list() %>% unformat_theta_output()),
                       component=j %>% as.character())
    return(distribution_parameters_mean_per_component)
  })

  # generate ellipse data
  ellipse_data <- distribution_parameters_list %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "component")))) %>%
    dplyr::summarize(generate_ellipse(mean_parameter %>% unlist(recursive = F), alpha = alpha, npoints = npoints))



  # create ellipse plot
  ellipse_plot <- ggplot(ellipse_data, aes(x=x, y=y)) +
    geom_path(aes(linetype=component, col=package), size=0.5) +
    theme_bw() +
    facet_wrap(~initialisation_method) +
    coord_fixed()+
    theme(legend.position = "bottom", title = element_blank(),
          legend.text = element_text(size = 18),
          axis.title.y = element_text(angle = 0, vjust=0.5),
          strip.text.x = element_text(size = 12, face = "bold")) +
    geom_path(data = ellipse_standard_data,mapping=aes(x=x, y=y, linetype=component),
              inherit.aes = F, col="black", size=2, alpha=0.25) +
    # scale_linetype_manual(name="component", values=1:2,
    #                       guide = guide_legend(override.aes = list(alpha = 0.25, shape=NA, size=1)))+
    guides(linetype = guide_legend(override.aes = list(alpha = 0.25, shape=NA, size=1)))+
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
#'
#' @return a ggplot object representing for each configuration of parameter, in each panel the
#' Hellinger distance for each package and initialisation method, respectively to each component
#'
#' @export
plot_Hellinger <- function(distribution_parameters, num_col=2) {
  # format the data
  true_theta <- distribution_parameters$true_parameters[[1]]; k <- length(true_theta$p)

  # format the data
  distribution_parameters_list <- distribution_parameters %>%
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(package=factor(package, levels = unique(package))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("package", "initialisation_method", "N.bootstrap"))),) %>%
    dplyr::summarise(list_parameter=list(stats::setNames(value_parameter, nm = name_parameter)%>%
                                           as.list() %>% unformat_theta_output()))

  # compute Hellinger distance per component
  hellinger_data <- purrr::map_dfr(1:k, function(j) {
    hellinger_data_per_component <- distribution_parameters_list %>%
      dplyr::mutate(mu=lapply(list_parameter, function(.x) .x$mu[,j]), component=as.character(j),
                    sigma=lapply(list_parameter, function(.x) .x$sigma[,,j]),
                    hellinger_value=purrr::map2_dbl(mu, sigma, hellinger, true_theta$mu[,j], true_theta$sigma[,,j])) %>%
      select(-c("list_parameter", "N.bootstrap", "mu", "sigma"))
    return(hellinger_data_per_component)
  })

  # generate Hellinger boxplot
  boxplot_Hellinger <- ggplot(hellinger_data, aes(x = factor(package, levels = unique(package)),
                                                  y = hellinger_value, fill = initialisation_method)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 1, colour = "red",
      position = position_dodge(width = 0.9), show.legend = FALSE) +
    facet_wrap(~component, ncol = num_col) +
    theme_bw() +
    theme(legend.position = "bottom",
      axis.text.x = element_text(angle = 35, size = 15, vjust = 0.6),
      axis.ticks.length = unit(.1, "cm"),
      legend.text = element_text(size = 25),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      title = element_blank(),
      panel.spacing = unit(.2, "pt")) +
    scale_fill_viridis_d()

  return(boxplot_Hellinger)
}



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

plot_time_computations <- function(time_data, grouping_colnames = c("package", "OVL", "entropy", "nobservations"),
                                   with_outliers=F){

  time_summary <- time_data %>% dplyr::group_by(across(all_of(grouping_colnames))) %>%
    dplyr::summarise(time_median= median(time), time_up=quantile(time, probs = 0.95), time_down=quantile(time, probs =0.05)) %>%
    # dplyr::mutate(OVL=factor(OVL, labels= paste("Balanced OVL:", unique(OVL)), levels= unique(sort(OVL))),
    #               entropy=factor(entropy, labels=paste("Entropy:", unique(entropy)), levels=unique(sort(entropy, decreasing = TRUE)))) %>%
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

  # if (with_outliers)
  #   time_plot <- time_plot +
  #   facet_grid(~ factor(prop_outliers, levels=unique(prop_outliers), labels=unique(paste("Proportion of outliers:", unique(prop_outliers)))))

  return(time_plot)
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

#' @rdname plot_univariate_normal_density_distribution
#' @export

plot_bivariate_normal_density_distribution <- function(theta, n=5000) {

  # retrieve the parameters of the multivariate distribution, and generate accordingly the ellipses
  k <- length(theta$p)
  theta_mean <- purrr::map_dfr(1:k, ~tibble::tibble(x=theta$mu[1,.x], y=theta$mu[2,.x], component=as.character(.x)))
  ellipse_standard_data <- purrr::map_dfr(1:k, function(j) {
    theta_per_component <- list(p=theta$p[j], mu=theta$mu[,j], sigma=theta$sigma[,,j])
    ellipse_per_component <- generate_ellipse(theta_per_component, alpha = 0.05, npoints = n) %>%
      tibble::add_column(component=as.character(j))
    return(ellipse_per_component)})
  ellipse_eigen_decomposition <- purrr::map_dfr(1:k, function(j) {
    var_per_comp <- theta$sigma[,,j]
    eigen_values <- eigen(var_per_comp, symmetric = T)$values
    # from function DrawingEllipsesinR, to get the theta main angle of rotation
    # and https://cookierobotics.com/007/
    ellipse_equation <- tibble::tibble(a=sqrt(eigen_values[1]), b=sqrt(eigen_values[2]),
                                       angle=atan2(eigen_values[1] - var_per_comp[1, 1], var_per_comp[1, 2])) %>%
      dplyr::bind_cols(theta_mean %>% filter(component==as.character(j)))
  })


  ellipse_axis <- purrr::map_dfr(1:k, function(j) {
    ellipse_decomposition <- eigen(theta$sigma[,,j], symmetric = T)
    ellipse_coordinates_per_co <- matrix(t(ellipse_decomposition$vectors) * sqrt(ellipse_decomposition$values),
                                         nrow=2, dimnames = list(c("comp1", "comp2"), c("xend", "yend"))) %>%
      tibble::as_tibble() %>%
      dplyr::bind_cols(theta_mean %>% filter(component==as.character(j)))
    return(ellipse_coordinates_per_co)}) %>% # account for mean shift
    dplyr::mutate(xend=x+xend, yend=y+yend)

  ellipse_labels <- purrr::map_dfr(1:k, function(j) {
    var_per_comp <- theta$sigma[,,j]
    ellipse_decomposition <- eigen(var_per_comp, symmetric = T); eigen_values <- ellipse_decomposition$values
    main_angle <- atan2(eigen_values[1] - var_per_comp[1, 1], var_per_comp[1, 2]) *360 /(2*pi)
    main_angle <- ifelse(abs(main_angle) > 90, main_angle -180,  main_angle)
    minor_angle <- ifelse(abs(main_angle + 90) > 90, main_angle - 90 , main_angle + 90)

    ellipse_coordinates_per_co <- matrix(t(ellipse_decomposition$vectors) * sqrt(eigen_values),
                                         nrow=2, dimnames = list(c("comp1", "comp2"), c("xend", "yend"))) %>%
      tibble::as_tibble() %>%
      dplyr::bind_cols(theta_mean %>% filter(component==as.character(j))) %>%
      dplyr::mutate(angle=c(main_angle, minor_angle), x=xend+x, y=yend + y,
                     eigen=eigen_values,
                     label=lapply(paste0("$\\sqrt{\\lambda ", 1:2, "}= ", round(sqrt(eigen_values), digits = 3)),
                                  latex2exp::TeX, output="character"))
    return(ellipse_coordinates_per_co)})


  # simulate the two-d density distribution
  # add text to characterize the ellipses?
  density_data <- simulate_multivariate_GMM (theta, n=n)$x
  colnames(density_data) <- c("x", "y")
  density_data <- density_data %>% tibble::as_tibble()


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
    scale_color_manual(name="component", values=c("red", "green"),
                       guide = guide_legend(override.aes = list(fill = NA,
                                                                shape = c(16, 17),
                                                                color = c("red", "green")))) +
    geom_path(data = ellipse_standard_data,mapping=aes(x=x, y=y, colour=component),
              inherit.aes = F, size=2, show.legend = F, linetype=2) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = component),
                 size = 1, arrow = arrow(length = unit(0.03, "inches")),
                 lineend = c('round'), linejoin = c('round'),
                 data=ellipse_axis, inherit.aes = F, show.legend = F) +
    ggforce::geom_ellipse(aes(x0 = x, y0 = y, a = a, b = b, angle = angle,
                              col=component), data=ellipse_eigen_decomposition,
                          show.legend = F, inherit.aes = F) +
    geom_text(mapping = aes(label=label, x=x, y=y, angle=angle), parse=T, col="black",
              data=ellipse_labels, show.legend = F, size=4, nudge_x = 0, nudge_y = 0.2, inherit.aes = F) +
    coord_fixed()


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
    tidyr::pivot_longer(dplyr::matches("p[[:digit:]]+|mu|sd"), names_to = "name_parameter",
                        values_to = "value_parameter") %>%
    dplyr::select(c("package", "initialisation_method", "name_parameter", "value_parameter", "N.bootstrap")) %>%
    dplyr::mutate(across(c("name_parameter", "package"), as.factor))
    # dplyr::group_by(OVL, entropy, package, name_parameter, initialisation_method) %>%
    # mutate(index=dplyr::row_number()) %>% dplyr::ungroup() # add an index, that uniquely identify each experiment


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
                                               width = unit(8, "cm"), height = unit(8, "cm"),
                                               column_title = paste("Initialization method:", init_method),
                                               column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"))
    return(grid::grid.grabExpr(ComplexHeatmap::draw(complex_heatmap)))
  }, total_correlation_scores_global, names(total_correlation_scores_global))

  return(total_correlation_scores_plots)
}
