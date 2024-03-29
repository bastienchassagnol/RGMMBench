% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixture.R
\name{emnmix_univariate}
\alias{emnmix_univariate}
\alias{em_Rmixmod_univariate}
\alias{em_EMCluster_univariate}
\alias{em_bgmm_univariate}
\alias{em_flexmix_univariate}
\alias{em_mixtools_univariate}
\alias{em_mclust_univariate}
\alias{em_GMKMcharlie_univariate}
\alias{em_otrimle}
\title{Custom R implementation of the EM algorithm in the univariate context}
\usage{
emnmix_univariate(
  x,
  k,
  epsilon = 10^-4,
  itmax = 500,
  nstart = 10L,
  start = NULL,
  initialisation_algorithm = "kmeans",
  ...
)

em_Rmixmod_univariate(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  epsilon = 10^-4,
  itmax = 500,
  start = NULL,
  ...
)

em_EMCluster_univariate(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  epsilon = 10^-4,
  itmax = 500,
  start = NULL,
  ...
)

em_bgmm_univariate(
  x = x,
  k = 2,
  epsilon = 10^-4,
  itmax = 500,
  initialisation_algorithm = "kmeans",
  start = NULL,
  ...
)

em_flexmix_univariate(
  x = x,
  k = 2,
  epsilon = 10^-4,
  itmax = 500,
  minprior = 0.05,
  initialisation_algorithm = "kmeans",
  start = NULL,
  ...
)

em_mixtools_univariate(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  epsilon = 10^-4,
  itmax = 500,
  start = NULL,
  ...
)

em_mclust_univariate(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  start = NULL,
  epsilon = 10^-4,
  itmax = 500,
  ...
)

em_GMKMcharlie_univariate(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  embedNoise = 1e-06,
  epsilon = 10^-4,
  itmax = 500,
  start = NULL,
  parallel = FALSE,
  ...
)

em_otrimle(
  x = x,
  k = 2,
  initialisation_algorithm = "kmeans",
  epsilon = 10^-4,
  itmax = 500,
  ...
)
}
\arguments{
\item{x}{the vector of the observations}

\item{k}{the number of components}

\item{epsilon}{the criterion threshold considered as the tolerance between two consecutive log-likelihoods}

\item{itmax}{the maximal number of iterations to reach the threshold}

\item{start}{a list of initial estimates provided by the user, with 3 entries:
\itemize{
\item The proportions \code{p}: \eqn{p} of each component (must be included between 0 and 1, and sum to one overall)
\item The mean matrix \code{mu}: \eqn{\mathrm{\mu}=(\mu_{i,j}) \in \mathbb{R}^{n \times k}}, with each column
giving the mean values of the variables within a given component
\item The 3-dimensional covariance matrix array \code{Sigma}: \eqn{\mathrm{\Sigma}=(\Sigma_{i,j,l}) \in \mathbb{R}^{n \times n \times k}}, with each matrix
\eqn{\Sigma_{..l}, l \in \{ 1, \ldots, k\}} storing the covariance matrix of a given component,
whose diagonal terms correspond to the variance of each variable, and off-terms diagonal elements return the covariance matrix
}}

\item{initialisation_algorithm, nstart}{hyper-parameters, when the user rather uses
one of our implemented initialization algorithms}

\item{...}{additional parameters for the reviewed packages}

\item{minprior}{Minimum prior probability of clusters, components falling below this threshold are removed during the iteration.}

\item{embedNoise}{
A small constant added to the diagonal entries of all covariance matrices. This may prevent covariance matrices collapsing prematurely. A suggested value is 1e-6. Covariance degeneration is detected during Cholesky decomposition, and will lead the trainer to remove the corresponding mixture component. For high-dimensional problem, setting \code{embedNoise} to nonzero may pose the illusion of massive log-likelihood, all because one or more mixture components are so close to singular, which makes the densities around them extremely high.
}

\item{parallel}{only relevant for GMKMCharlie package which has a native parallel implementation (by default, takes half of the available clusters)}
}
\value{
a list of the estimated parameters, ordered by increasing mean for identifiability issues
}
\description{
Custom R implementation of the EM algorithm in the univariate context
}
\seealso{
\code{\link[=emnmix_multivariate]{emnmix_multivariate()}}
}
\author{
Bastien CHASSAGNOL
}
