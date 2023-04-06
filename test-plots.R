
library(dplyr)
installed_packages <- pkgapi::map_package(path = ".")
installed_packages_plot_data <- installed_packages$calls %>%
  mutate(package=stringr::str_extract(to, "(.*)(?=::)"),
         file=tools::file_path_sans_ext(file) %>% basename(),
         shape=forcats::fct_recode(file, triangle = "benchmark", ellipse = "visualisation",
                                   square = "mixture",  star = "utils")) %>%
  select(from, str, package, shape) %>% rename(to=str)  %>% filter(package=="RGMMBench")



unique_functions <- tibble::as_tibble_col(union(installed_packages_plot_data$from, installed_packages_plot_data$to),
                                          column_name = "id")

nodes_visNet <- tibble::tibble( "id"=unique_functions$id, label = unique_functions$id,
                                shape=unique_functions %>% dplyr::left_join(installed_packages_plot_data %>% select(-c(to, package)),
                                                                            by = c("id"="from")) %>%  dplyr::distinct() %>%
                                  pull(shape),
                                title = paste0("<p><b>", unique_functions$id,"</b><br></p>"))

visNetwork::visNetwork(nodes= nodes_visNet, height = "500px",
                       edges=installed_packages_plot_data %>% select(-package), width = "100%") %>%
  visNetwork::visEdges(arrows =list(to = list(enabled = TRUE, scaleFactor = 2)),
           color = list(color = "lightblue", highlight = "red"))

# graph_packages <- igraph::graph_from_data_frame(installed_packages_plot_data %>% select(-package),
#                                                 directed = TRUE)




