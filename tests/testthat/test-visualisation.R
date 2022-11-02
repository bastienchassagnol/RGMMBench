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
library(dplyr); library(ggplot2)
distribution_parameters <- readRDS("../mixture_models/results/reduced_parameters_distribution_parallel.rds")$distributions
time_computations <- readRDS("../mixture_models/results/reduced_time_computation_verbose.rds")$time_data


#################################################################
##          negative correlated, strong OVL, unbalanced          ##
#################################################################
unique_scenarios <- unique(distribution_parameters$true_parameters_factor)
specific_distribution_parameters <- distribution_parameters %>%
  filter(true_parameters_factor==unique_scenarios[11])
specific_time_computations <- time_computations %>%
  filter(true_parameters_factor=="p1_0.9_p2_0.1_mu_var1_comp1_0_mu_var1_comp2_2_mu_var2_comp1_2_mu_var2_comp2_0_sd_var1_var1_comp1_1_sd_var2_var1_comp1_-0.8_sd_var2_var2_comp1_1_sd_var1_var1_comp2_1_sd_var2_var1_comp2_-0.8_sd_var2_var2_comp2_1")
theta <- specific_distribution_parameters$true_parameters[[1]]



bivariate_boxplot <- plot_boxplots_parameters_bivariate(distribution_parameters, num_col = 4)
bivariate_ellipse <- plot_ellipses_bivariate(distribution_parameters)
time_computations <- plot_time_computations(specific_time_computations)






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



distribution_parameters <- readRDS("../mixture_models/results/reduced_parameters_distribution_parallel.rds")$distributions
unique_scenarios <- unique(distribution_parameters$true_parameters_factor)
specific_distribution_parameters <- distribution_parameters %>%
  filter(true_parameters_factor==unique_scenarios[11])

correlation_scores_plots <- specific_distribution_parameters  %>% plot_correlation_Heatmap()
# correlation_scores_plots <- gridExtra::arrangeGrob(grobs=correlation_scores_plots, top="", ncol = 2)
ggsave("./images/highly_correlation_scores.pdf", correlation_scores_plots, width = 10.5, height = 13.5,dpi = 300)









