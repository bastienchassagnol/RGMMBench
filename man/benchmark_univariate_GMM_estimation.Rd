% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{benchmark_univariate_GMM_estimation}
\alias{benchmark_univariate_GMM_estimation}
\alias{benchmark_multivariate_GMM_estimation}
\title{Launch the benchmark to compare statistical performances between packages}
\usage{
benchmark_univariate_GMM_estimation(
  mixture_functions,
  sigma_values,
  mean_values,
  proportions,
  cores = getOption("mc.cores", parallel::detectCores()),
  id_scenario = NULL,
  prop_outliers = 0,
  nobservations = c(2000),
  Nbootstrap = 100,
  epsilon = 10^-4,
  itmax = 500,
  nstart = 10L,
  short_iter = 200,
  short_eps = 10^-2,
  prior_prob = 0.05,
  initialisation_algorithms = c("kmeans", "quantiles", "random", "hc", "rebmix")
)

benchmark_multivariate_GMM_estimation(
  mixture_functions,
  mean_values,
  proportions,
  sigma_values,
  id_scenario = NULL,
  cores = getOption("mc.cores", parallel::detectCores()),
  nobservations = c(2000),
  Nbootstrap = 100,
  epsilon = 10^-4,
  itmax = 500,
  nstart = 10L,
  short_iter = 200,
  short_eps = 10^-2,
  prior_prob = 0.05,
  initialisation_algorithms = c("kmeans", "random", "hc", "rebmix")
)
}
\arguments{
\item{mixture_functions}{List of the packages to be compared (Id:name of the package, value: its options)}

\item{sigma_values, mean_values, proportions}{the true parameters to be retrieved}

\item{cores}{the number of cores to be used, by default all the available cores}

\item{id_scenario}{Possibility to set it to another number than one, to uniquely identify them}

\item{prop_outliers}{the proportion of outliers added in the simulation}

\item{nobservations}{the number of observations drawn to generate the random sample}

\item{Nbootstrap}{the number of bootstrap simulations and repetitions to perform}

\item{epsilon, itmax}{respectively criterion threshold and maximal number of iterations to reach it}

\item{nstart, short_iter, short_eps}{hyper-parameters to control the initialisation step}

\item{prior_prob}{add minimal uncertainty on the cluster assignment returned by hierarchical clustering method}

\item{initialisation_algorithms}{among 6 methods, which algorithms to be chosen for the initialisation phase}
}
\value{
a list with the simulated distributions of the estimates, some summary scores per parameter and aggregated measures
as well as boxplot and Heatmap correlation representations of the estimates
}
\description{
Launch the benchmark to compare statistical performances between packages
}
\author{
Bastien CHASSAGNOL
}
