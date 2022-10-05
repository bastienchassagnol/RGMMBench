library(ggplot2); library(dplyr)

distribution_parameters_files <- list.files("../mixture_models//multivariate_parallel_2022-09-16_20-17-30/job_1/results/",
                                            pattern = "\\.rds$",full.names = T)

distribution_parameters <- purrr::map_dfr(distribution_parameters_files, ~readRDS(.x)[["distribution"]]) %>%
  select(-c(nobservations, N.bootstrap))
# distribution_parameters %>% write.csv("./results/distribution_multivariate.csv",  row.names=F)
# distribution_parameters <- readr::read_csv("./results/distribution_multivariate.csv")


distribution_parameters_long <- distribution_parameters %>%
  tidyr::pivot_longer(dplyr::matches("^(p(1|2)+|mu|sd)"), names_to = "name_parameter",
                      values_to = "value_parameter") %>%
  dplyr::mutate(across(c("name_parameter", "package"), as.factor))

correlation_scores_plots <- distribution_parameters_long  %>% plot_correlation_Heatmap()
correlation_scores_plots <- gridExtra::arrangeGrob(grobs=correlation_scores_plots, top="", ncol = 2)
ggsave("./images/global_correlation_scores.pdf",
       correlation_scores_plots, width = 10.5, height = 13.5,dpi = 300)




##################################################################
##                      boxplot parameters                      ##
##################################################################
distribution_parameters <- readRDS("../mixture_models/results/reduced_parameters_distribution_parallel_corrected.rds")$distributions %>%
  filter(true_parameters_factor==true_parameters_factor[1])
theta <- distribution_parameters$true_parameters[[1]]

bivariate_boxplot <- plot_boxplots_parameters_bivariate(distribution_parameters, num_col = 4)
bivariate_ellipse <- plot_ellipses_bivariate(distribution_parameters)


density_plots_four <- plot_density_distribution(sigma_values=list("null OVL"=rep(0.3, 4), "aevrage OVL"=rep(1, 4), "high OVL"= rep(2, 4)),
                                                mean_values=list(c(0, 4, 8, 12)), skewness_values = list("null skewness"=rep(0, 4)),
                                                proportions=list("balanced"=rep(1/4, 4), "small unbalance"= c(0.2, 0.4, 0.2, 0.2), "highly unbalanced"=c(0.1, 0.7, 0.1, 0.1)))
## four components unbalanced overlapping
upper_part <- cowplot::plot_grid(density_plots_four$individual_plots$UR_0.9_skewness_0_OVL_0.08 + labs(tag="A") +
                                   theme(plot.tag = element_text(size=18, face = "bold")),
                                 plot_time_computations(unbalanced_overlapping_microbenchmark_computations$init_focus_time_plots$kmeans$data)$kmeans +
                                   theme(strip.text = element_blank(), plot.title = element_blank(), plot.tag = element_text(size=18, face = "bold"),
                                         legend.text = element_text(size=20), axis.title=element_text(size=20)) +
                                   labs(tag = "B", x="Number of observations (log10)", y="Time in seconds (log10)"),
                                 ncol = 2, align = 'h', axis="tblr")
four_components_unbalanced_overlapping <- gridExtra::arrangeGrob(upper_part,
                                                                 plot_boxplots_parameters(benchmark_scores_2000_observations_UR_0.9_skewness_0_OVL_0.08_prop_outliers_0$plots$data, p=c(0.1, 0.7, 0.1, 0.1), mu= c(0, 4, 8, 12), sigma=rep(2, 4)) +
                                                                   theme(plot.tag = element_text(size=18, face = "bold")) + labs(tag = "C"),
                                                                 top="", nrow = 2, heights = c(1.2, 2), padding = unit(1, "line"))
ggsave("paper_figures/four_components_unbalanced_overlapping.pdf", four_components_unbalanced_overlapping, width = 20, height = 18,dpi = 600)





distributions_plots <- purrr::imap(split(distribution_parameters, distribution_parameters$true_parameters_factor),
                                   function(dist_per_fact, title) {
                                return(list("boxplot_parameter"=plot_boxplots_parameters_bivariate(dist_per_fact, num_col = 4) + ggtitle(title),
                                "ellipse"=plot_ellipses_bivariate(dist_per_fact) + ggtitle(title),
                                "hellinger"=plot_Hellinger(dist_per_fact) + ggtitle(title))) })

ellipse_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("ellipse"), ncol = 1, nrow = 1, top=NULL)
ggsave("./images/ellipse_plots.pdf", ellipse_plots, width = 15, height = 8,dpi = 300)

boxplot_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("boxplot_parameter"), ncol = 1, nrow = 1, top=NULL)
ggsave("./images/boxplot_parameters.pdf", boxplot_plots, width = 15, height = 15,dpi = 300)

hellinger_plots <- gridExtra::marrangeGrob(distributions_plots %>% purrr::map("hellinger"), ncol = 1, nrow = 1, top=NULL)
ggsave("./images/hellinger_plot.pdf", hellinger_plots, width = 15, height = 15,dpi = 300)

hellinger_plot <- plot_Hellinger(splitted_distribution[[11]])
ggsave("./images/hellinger_plot_specific.pdf", hellinger_plot, width = 15, height = 15,dpi = 300)


distribution_parameters_long <- splitted_distribution[[11]] %>%
  tidyr::pivot_longer(dplyr::matches("^(p(1|2)+|mu|sd)"), names_to = "name_parameter",
                      values_to = "value_parameter") %>%
  dplyr::mutate(across(c("name_parameter", "package"), as.factor))

correlation_scores_plots <- distribution_parameters_long  %>% plot_correlation_Heatmap()
correlation_scores_plots <- gridExtra::arrangeGrob(grobs=correlation_scores_plots, top="", ncol = 2)
ggsave("./images/highly_correlation_scores.pdf",
       correlation_scores_plots, width = 10.5, height = 13.5,dpi = 300)









