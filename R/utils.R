#' Compute the intersection points between two Gaussian mixtures (where one component goes above the other)
#'
#' @author Bastien CHASSAGNOL
#'
#' @param p1,p2,mu1,mu2,sigma1,sigma2 parameters of the two Gaussian components
#' @return the overlap roots
#' @export

compute_overlap_roots <- function(p1, p2, mu1, mu2, sigma1, sigma2) {
  # homoscedastic case, only one intersection point
  if (sigma1 == sigma2) {
    return((sigma1^2 * log(p2 / p1)) / (mu1 - mu2) + (mu1 + mu2) / 2)
  }
  # heteroscedastic case, existence of two intersections points
  # conditioned on the proportions
  else if ((sigma2 > sigma1 & p1 > (sigma1 / (sigma1 + sigma2))) |
           (sigma2 < sigma1 & p1 < (sigma1 / (sigma1 + sigma2)))) {
    variant_part <- sigma1 * sigma2 * sqrt((mu1 - mu2)^2 + 2 * (sigma2^2 - sigma1^2) * (log(p1 / p2) + log(sigma2 / sigma1)))
    return(c(
      (sigma1^2 * mu2 - sigma2^2 * mu1 - variant_part) / (sigma1^2 - sigma2^2),
      (sigma1^2 * mu2 - sigma2^2 * mu1 + variant_part) / (sigma1^2 - sigma2^2)
    ))
  } else {
    stop("There's no intersection point in that configuration")
  }
}

#' Compute a set of points, in order to draw 95% confidence intervals
#'
#' @author Bastien CHASSAGNOL
#'
#' @param theta the list of the estimates returned by the EM algorithm
#' @param alpha the confidence interval
#' @param npoints the number of points used to generate the ellipes
#' @return a tibble with `n` points associated to x and y locations associated to the `1-alpha` confidence interval


generate_ellipse <- function (theta, alpha = 0.05, npoints = 500) {

# control input
mu <- theta$mu %>% as.vector(); sigma <- theta$sigma
if(!is.matrix(sigma))
  sigma <- sigma[,,1]


es <- eigen(sigma)
e1 <- es$vec %*% diag(sqrt(es$val))
r1 <- sqrt(qchisq(1 - alpha, 2))
theta <- seq(0, 2 * pi, len = npoints)
v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
pts <- t(mu - (e1 %*% t(v1))) %>% tibble::as_tibble()
colnames(pts) <- c("x", "y")
return(pts)
}



#' Compute bias and mse for each parameter of a GMM distribution
#'
#' @author Bastien CHASSAGNOL
#'
#' @param estimated_theta the list of the estimates returned by the EM algorithm
#' @param true_theta the true parameters
#' @return a tibble with mean, standard deviation, bias and mse for each parameter
#' @export


get_local_scores <- function(estimated_theta, true_theta) {
  # names_param <- names(unlist(true_theta[c("p", "mu", "sigma")]))
  mean_parameter <- apply(estimated_theta, 2, mean) # mean of the distribution
  empiric_sds <- apply(estimated_theta, 2, sd) # sd of the distribution
  bias <- apply(estimated_theta, 2, mean) - unlist(true_theta) # bias of the distribution
  mse <- apply(estimated_theta, 2, var) + bias^2 # mse of the distribution
  return(tibble::tibble(
    scores = c("mean", "sd", "bias", "mse"),
    dplyr::bind_rows(mean_parameter, empiric_sds, bias, mse)
  ))
}

#' The Maha distance
maha <- function(mu, sigma) {
  inv.sigma <- solve(sigma)
  d <- t(mu) %*% inv.sigma %*% mu
  return(as.vector(d))
}

#' Compute the shannon entropy of a discrete distribution, normalised from 0 to 1 (equibalanced classes)
#'
#' @author Bastien CHASSAGNOL
#'
#' @param ratios vector of the proportions of the mixture
#' @return the entropy score
#' @export


compute_shannon_entropy <- function(ratios) {
  if (min(ratios) < 0 | max(ratios) > 1) {
    stop("Probabilities must be stricly included between 0 and 1")
  }

  # normalization process + remove NULL components, as information fisher is not modified by empty classes
  ratios <- ratios[ratios != 0]
  ratios <- ratios / sum(ratios)

  # entropy included between 0 (one component storing all information) and 1 (uniform distribution, balanced classes)
  return(-sum(ratios * logb(ratios, base = length(ratios))))
}



#' Check the positive definiteness of a symmetric real matrix
#'
#'To do so, we compute the eigen values of the corresponding matrix with `eigen` function,
#'and if all of them are positive above a given threshold, then we return true.
#' @author Bastien CHASSAGNOL
#'
#' @param sigma a symmetric matrix with real values
#' @param tol the numerical maximal tolerance threshold,
#' to which an eigen value below it is considered negative
#'
#' @return a boolean, whether or not the matrix can be considered
#' positive definite or not
#' @export
#'

is_positive_definite <- function(sigma, tol=1e-6) {
  eigen_values <- eigen(sigma, symmetric = TRUE)$values
  return(all(eigen_values >= -tol))
}

#' Convert vector to symmetric matrix
#'
#' Reverse operation of \code{upper.tri} or \code{lower.tri}. Given a vector of
#' upper or lower triangular elements of a matrix
#' (optionally including the diagonal elements), returns the
#' corresponding symmetric matrix. The elements of the vector can be arranged by
#' column (the default) or by row.
#'
#' @details Many thanks for the author of the `patr1ckm` package for providing us this method,
#' [Handy R functions](https://rdrr.io/github/patr1ckm/patr1ckm/man/vec2sym.html),
#'  enabling to revert an `upper.tri` or `lower.tri` operation.
#'
#' @param x vector containing upper or lower triangular elements of a symmetric matrix
#' @param If NULL, x contains diagonal elements. Otherwise a value or vector of the appropriate length to be placed on the diagonal.
#' @param lower x is from the lower triangle (default = \code{TRUE})
#' @param byrow the elements in x are ordered row-wise (default = \code{FALSE})
#' @return Symmetric matrix
#' @details Note that if x is a vector of the lower triagular elements given by column,
#' this is equivalent to the upper triangular elements given by row. Similarly, if
#' x is a vector of the lower triangular elements given by row, this is equivalent
#' to the upper triangular elements given by column.
#' @examples
#' x <- c(1,1,1,2,2,3)
#' check <- matrix(c(1,1,1,1,1,1,2,2,1,2,1,3,1,2,3,1),4,4)
#' identical(vec2sym(x,diag=1,lower=T,byrow=F), check)
#' x <- c(1,1,1,2,2,3)
#' check <- matrix(c(1,1,1,1,2,2,1,2,3),3,3)
#' identical(vec2sym(x,lower=T,byrow=F),check)
#' @author patr1ckm
#' @export
vec2sym <- function(x,diagonal=NULL,lower=TRUE,byrow=FALSE){
  if(is.null(diagonal)){
    ## Assume that x contains the diagonal elements as well
    p <- (sqrt(1 + 8 * length(x)) - 1)/2
    S <- diag(p)
    if((!lower & byrow) | (lower & !byrow)){
      S[lower.tri(S,diag=T)] <- x
      S[which(lower.tri(t(S),diag=T),arr.ind=T)[,c(2,1)]] <- x
    }
    if((lower & byrow) | (!lower & !byrow)) {
      S[upper.tri(S,diag=T)] <- x
      S[which(upper.tri(t(S),diag=T),arr.ind=T)[,c(2,1)]] <- x
    }
  } else{
    ## diagonal elements are given by 'diagonal'
    p <- (sqrt(1 + 8 * length(x)) + 1)/2
    S <- diag(p)
    if((!lower & byrow) | (lower & !byrow)){
      S[lower.tri(S)] <- x
      S[which(lower.tri(t(S)),arr.ind=T)[,c(2,1)]] <- x
    }
    if((lower & byrow) | (!lower & !byrow)) {
      S[upper.tri(S)] <- x
      S[which(upper.tri(t(S)),arr.ind=T)[,c(2,1)]] <- x
    }
    diag(S) <- diagonal
  }

  return(S)
}


#' Convert condensed covariance matrix to large array format
#'
#' Shift from the lower triangular matrices \eqn{k \times \frac{p(p+1)}{2}} to large array, of dimension
#' \eqn{p \times p \times k}
#'
#' @param x a covariance matrix, formatted as \eqn{\frac{p(p+1)}{2} \times k},
#' each column \eqn{j \in \{1, \ldots, k \}} describing the covariance matrix associated to component \eqn{j}
#' @param transposed if transposed, the input is of dimension \eqn{k \times \frac{p(p+1)}{2}},
#' each row \eqn{j} describing the covariance matrix associated to component \eqn{j}
#' @param ... additional parameters to pass to function [vec2sym()], enabling to deal with
#' distinct way to provide covariance information
#'
#' @return a large array, of dimension \eqn{p \times p \times k}: one array of covariance matrix,
#' on the third dimension, for each component
#' @export

trig_mat_to_array <- function(x, transposed=TRUE, ...) {

  if (!transposed)
    x <- t(x)

  k <- nrow(x); dim_gaussian <- (sqrt(8 * ncol(x) + 1) - 1)/2
  stopifnot("Dimension of the input must be a Pascal number (to have a triangular valid matrix)" =
              is_integer(dim_gaussian))
  sigma <- array(0, dim=c(dim_gaussian, dim_gaussian, k))

  for (j in 1:k) {
      sigma[,,j] <- vec2sym(x[j,], ...)
  }

  return(sigma)
}

#' Convert an array of covariance matrices to short lower triangular format
#'
#' @seealso [trig_mat_to_array()], which is the reciprocal operation
#' @export

array_to_trig_mat <- function(x, transposed=TRUE) {
  k <- dim(x)[3]; dim_gaussian <- dim(x)[2]
  dim_triangular <- ((dim_gaussian + 1) * dim_gaussian)/2

  trig_matrix <- matrix(0, nrow = k, ncol = dim_triangular)
  for (j in 1:k) {
    # recover elements from the lower traingular part of the covariance matrix, left to right, up to bottom
    trig_matrix[j, ] <- x[,,j][lower.tri(x[,,j], diag = TRUE)] %>% as.vector()
  }
  if (!transposed)
    trig_matrix <- t(trig_matrix)
  return(trig_matrix)
}

#' Test whether a number is an integer
#'
#' @param x (a numeric input of size 1 (real number))
#'
#' @return F if the Euclidean division by 1 does not return 0

is_integer <- function(x) {
  return(x%%1 == 0)
}

#' Control parameters output
#'
#' This step ensures that the estimates returned are uniquely ordered by
#' partial ordering on the means, and that the sum-o-one constraint, that may
#' be violated by numerical artefacts, is enforced
#'
#'
#' @param theta estimation of the parameters returned either by an initialisation
#' algorithm or by an EM algorithm on an univariate or multivariate GMM MLE estimation
#' * The proportions `p`: \eqn{p} of each component (must be included between 0 and 1, and sum to one overall)
#' * The mean matrix `mu`: \eqn{\mathrm{\mu}=(\mu_{i,j}) \in \mathbb{R}^{n \times k}}, with each column
#' giving the mean values of the variables within a given component
#' * The 3-dimensional covariance matrix array `Sigma`: \eqn{\mathrm{\Sigma}=(\Sigma_{i,j,l}) \in \mathbb{R}^{n \times n \times k}}, with each matrix
#' \eqn{\Sigma_{..l}, l \in \{ 1, \ldots, k\}} storing the covariance matrix of a given component,
#' whose diagonal terms correspond to the variance of each variable, and off-terms diagonal elements return the covariance matrix
#' @param n the number of observations to be drawn
#'
#' @return a list of the estimates, uniquely identified, by ranking each component
#' based on the ordering of their means

enforce_identifiability <- function(theta) {
  k <- length(theta$p)

  if (is.array(theta$sigma)) {
    # in that case, we are in a multivariate context
    ordered_components <- do.call(order, theta$mu %>% as.data.frame())
    ordered_theta <- list(
      p = theta$p[ordered_components],
      mu = theta$mu[, ordered_components],
      sigma = theta$sigma[,,ordered_components]
    )
  }
  else {
    # univariate context
    ordered_theta <- list(
      p = theta$p[order(theta$mu)],
      mu = sort(theta$mu),
      sigma = theta$sigma[order(theta$mu)]
    )
  }

  # enforce sum-to-one constraint
  ordered_theta <- ordered_theta %>% purrr::map(function(x) unname(x))
  ordered_theta$p <- ordered_theta$p / sum(ordered_theta$p)
  ordered_theta$p[k] <- 1 - sum(ordered_theta$p[-k])
  return(ordered_theta)
}

test_function <- function(x, w) {
  return (apply(x, 2, function(column) {
    prin(column)
    stats::weighted.mean(column, w)}
    ))
}


#' Format the estimated parameters
#'
#' Especially, we remove redundant pairwise correlations, and name adequatly each parameter
#'
#' @author Bastien CHASSAGNOL
#'
#' @param theta a list with 3 entries:
#' * The proportions `p`: \eqn{p} of each component (must be included between 0 and 1, and sum to one overall)
#' * The mean matrix `mu`: \eqn{\mathrm{\mu}=(\mu_{i,j}) \in \mathbb{R}^{n \times k}}, with each column
#' giving the mean values of the variables within a given component
#' * The 3-dimensional covariance matrix array `Sigma`: \eqn{\mathrm{\Sigma}=(\Sigma_{i,j,l}) \in \mathbb{R}^{n \times n \times k}},
#'  with each matrix \eqn{\Sigma_{..l}, l \in \{ 1, \ldots, k\}} storing the covariance matrix of a given component,
#' whose diagonal terms correspond to the variance of each variable, and off-terms diagonal elements return the covariance matrix
#' @export


format_theta_output <- function(theta) {
  mu <- theta$mu; sigma <- theta$sigma; p <- theta$p
  if(!is.matrix(mu)) { # univariate setting
    return(unlist(theta))
  }
  else { # multivariate setting
  dim_gaussian <- nrow(mu); k <- length(p)

  formatted_p <- stats::setNames(p, nm=paste0("p", 1:k))
  formatted_mu <- stats::setNames(mu %>% t() %>% as.vector(),
                                  nm = tidyr::crossing(variable=paste0("mu_var", 1:dim_gaussian), component=paste0("comp", 1:k)) %>%
                                    tidyr::unite(col="names_mu", variable, component) %>% dplyr::pull(names_mu))

  names_sigma <- c()
  for (j in 1:k) {
    for(i in 1:dim_gaussian) {
      for (l in i:dim_gaussian) {
        names_sigma <- c(names_sigma, glue::glue("sd_var{l}_var{i}_comp{j}"))
      }
    }
  }
  formatted_sigma <- stats::setNames(sigma %>% array_to_trig_mat(transposed = F) %>% as.vector(), names_sigma)
  return(c(formatted_p, formatted_mu, formatted_sigma))
  }
}

#' @rdname format_theta_output
unformat_theta_output <- function(formatted_theta) {
names_theta <- names(formatted_theta); theta <- list()
# deal with proportions
theta$p <- formatted_theta[stringr::str_detect(names_theta, "^p[[:digit:]]+$")] %>% unlist() %>% unname()
k <- stringr::str_detect(names_theta, "^p[[:digit:]]+$") %>% sum()

# deal specifically per component parameter's distribution
if (any(grepl("sd", names_theta))) { # multivariate setting
  dim_gaussian <- stringr::str_detect(names_theta, "mu") %>% sum() / k
  # deal with mean vector
  theta$mu <- matrix(formatted_theta[stringr::str_detect(names_theta, "mu")] %>% unlist(), nrow=dim_gaussian, ncol=k)
  # deal with sigma
  theta$sigma <- matrix(formatted_theta[stringr::str_detect(names_theta, "sd")] %>% unlist(), nrow=k, byrow = T) %>% trig_mat_to_array()
}
else { # univariate setting
  # deal with mean vector
  theta$mu <- formatted_theta[stringr::str_detect(names_theta, "mu")] %>% unlist() %>% unname()
  # deal with sigma
  theta$sigma <- formatted_theta[stringr::str_detect(names_theta, "sigma")] %>% unlist() %>% unname()
}
return(theta)
}



#' Compute the Hellinger distance
#'
#' @return the Hellinger distance between two multivariate Gaussian distributions
hellinger <- function (mu1, Sigma1, mu2, Sigma2) {
  p <- length(mu1);   d <- mu1 - mu2
  vars <- (Sigma1 + Sigma2)/2
  # in univariate dimension
  if (p == 1) {
    d <- sqrt(Sigma1 * Sigma2/ vars) *
      exp((-1/4) * d^2/(2*vars))
    return(sqrt(1 - d))
  }
  # in multivariate dimension
  else {
    hell_dist <- det(Sigma1)^(1/4) * det(Sigma2)^(1/4) / det(vars)^(1/2) *
      exp((-1/8) * maha(d, vars))
    return(sqrt(1 - hell_dist) %>% as.numeric())
  }
}





# # aggregate scores obtained per parameter
# global_bias <- local_scores_temp %>%
#   dplyr::filter(scores == "bias") %>%
#   dplyr::rowwise() %>%
#   dplyr::transmute(global_bias = sum(abs(dplyr::c_across(dplyr::all_of(bootstrap_colnames))))) # global bias
# global_mse <- local_scores_temp %>%
#   dplyr::filter(scores == "mse") %>%
#   dplyr::rowwise() %>%
#   dplyr::transmute(global_mse = sum(dplyr::c_across(dplyr::all_of(bootstrap_colnames)))) # global mse
# missed_cases <- distribution_parameters_temp %>%
#   dplyr::count(initialisation_method, package, name = "N.missed") %>%
#   dplyr::mutate(N.missed = Nbootstrap - N.missed) # count number of times the package failed in the estimation with the given number of clusters
#
# global_scores_temp <- global_scores_temp %>%
#   dplyr::group_by(initialisation_method, package) %>%
#   dplyr::summarise(
#     logLikelihood = mean(logLikelihood),
#     entropy = signif(entropy_value, digits = 2),
#     OVL = signif(balanced_ovl, digits = 2), OVL_pairwise = signif(pairwise_ovl, digits = 2), skew = skew[1], nobservations = n, prop_outliers = prop_out
#   ) %>%
#   dplyr::inner_join(missed_cases, by = c("initialisation_method", "package")) %>%
#   dplyr::inner_join(global_bias, by = c("initialisation_method", "package")) %>%
#   dplyr::inner_join(global_mse, by = c("initialisation_method", "package"))
# global_scores <- global_scores %>% dplyr::bind_rows(global_scores_temp)




#' #' @describeIn emnmix_univariate EM implementation with mixsmsn package (designed to deal with skewed GMMs especially)
#' #' @export
#' em_mixsmsn <- function(x = x, k = 2, initialisation_algorithm = "hc", skew = rep(0, k),
#'                        itmax = 5000, epsilon = 10^-12, start = NULL, ...) {
#'   if (is.null(start)) {
#'     start <- initialize_em_univariate(
#'       x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
#'       initialisation_algorithm = initialisation_algorithm, ...
#'     )
#'   }
#'
#'   # nu is the kurtosis, always equal to 3 for a classic Gaussian distribution(not skewed)
#'   if (all(skew == 0)) {
#'     # fit a classic Gaussian distribution
#'     fit <- mixsmsn::smsn.mix(
#'       y = x, nu = 3,
#'       mu = start$mu, sigma2 = start$sigma^2, shape = skew, pii = start$p,
#'       g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Normal",
#'       error = 10^-12, iter.max = 5000, calc.im = FALSE
#'     )
#'   } else {
#'     # fit a skewed distribution
#'     fit <- mixsmsn::smsn.mix(
#'       y = x, nu = 3,
#'       mu = start$mu, sigma2 = start$sigma, shape = skew, pii = start$p,
#'       g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Skew.normal",
#'       error = 10^-12, iter.max = 5000, calc.im = FALSE
#'     )
#'   }
#'
#'
#'   # return an ordered list by mean values
#'   ordered_estimated_theta <- list(
#'     p = fit$pii[order(fit$mu)],
#'     mu = sort(fit$mu),
#'     sigma = sqrt(fit$sigma2[order(fit$mu)])
#'   )
#'
#'   ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
#'   return(ordered_estimated_theta)
#' }





