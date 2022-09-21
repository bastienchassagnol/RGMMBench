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

quantiles_0_95 <- mvtnorm::qmvnorm(p=0.95, tail = c("both.tails"), mean = c(0, 2),
                                   sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2))
