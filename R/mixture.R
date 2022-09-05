#' Generation of a mixture drawn from an univariate GMM, with possibility to add outliers and skewness
#'
#' @details We do not implement the possibility of adding outliers in the multivariate context
#'
#' @author Bastien CHASSAGNOL
#'
#' @param theta a list with 3 entries:
#' * The proportions `p`: \eqn{p} of each component (must be included between 0 and 1, and sum to one overall)
#' * The mean matrix `mu`: \eqn{\mathrm{\mu}=(\mu_{i,j}) \in \mathbb{R}^{n \times k}}, with each column
#' giving the mean values of the variables within a given component
#' * The 3-dimensional covariance matrix array `Sigma`: \eqn{\mathrm{\Sigma}=(\Sigma_{i,j,l}) \in \mathbb{R}^{n \times n \times k}}, with each matrix
#' \eqn{\Sigma_{..l}, l \in \{ 1, \ldots, k\}} storing the covariance matrix of a given component,
#' whose diagonal terms correspond to the variance of each variable, and off-terms diagonal elements return the covariance matrix
#' @param n the number of observations to be drawn
#'
#' @return Depending on the univariate or multivariate context:
#' * a list with the number of components k, the true parameters p, mu, sigma,
#' the observed variables x, the hidden observations s and an indicator of the outliers s_outliers
#' * a list with the number of components k, the true parameters p, mu, sigma,
#' the observations `X` living in space \eqn{\mathbb{R}^k}, the hidden indicator variables `s` (vector of size \eqn{n})
#'
#' @export


simulate_univariate_GMM <- function(theta = list(p = c(0.40, 0.60), mu = c(175, 165), sigma = c(10, 12)), n=100,
                                       prop_outliers = 0, interval = 2) {
  # get values from theta parameter
  p <- theta$p
  mu <- theta$mu
  sigma <- theta$sigma
  k <- length(p)

  # generate hidden variables s set
  s <- sample(1:k, size = n, replace = TRUE, prob = p)

  # generate observed variable set
  x <- stats::rnorm(n, mean = mu[s], sd = sigma[s])

  # select randomly prop_outliers points to be drawn for uniform distribution
  # choice of points such that they are outside quantiles 0.025 of each distribution
  outliers_indexes <- sample(1:n, size = round(prop_outliers * n), replace = FALSE)
  s_outliers <- s
  s_outliers[outliers_indexes] <- 0

  # choice of enough ranged random distribution
  min_component <- which.min(sapply(1:k, function(j) stats::qnorm(0.05, mean = mu[j], sd = sigma[j])))
  max_component <- which.max(sapply(1:k, function(j) stats::qnorm(0.95, mean = mu[j], sd = sigma[j])))
  length_interval <- interval * (stats::qnorm(0.95, mean = mu[max_component], sd = sigma[max_component]) -
    stats::qnorm(0.05, mean = mu[min_component], sd = sigma[min_component]))
  x[outliers_indexes] <- runif(
    n = round(prop_outliers * n),
    min = stats::qnorm(0.01, mean = mu[min_component], sd = sigma[min_component]) - length_interval,
    max = stats::qnorm(0.99, mean = mu[max_component], sd = sigma[max_component]) + length_interval
  )


  return(list(k = k, p = p, mu = mu, sigma = sigma, x = x, s = s, s_outliers = s_outliers))
}


#' @rdname simulate_univariate_GMM
#' @export


simulate_multivariate_GMM <- function(theta, n=500) {
  ##################################################################
  ##                        check validity                        ##
  ##################################################################

  p <- theta$p; mu <- theta$mu; sigma <- theta$sigma; k <- length(p) # get values from theta parameter
  dimension_gaussian <- dim(mu)[1]


  stopifnot("Parameters provided do not correspond to a fitted GMM estimation"=
              check_parameters_validity_multivariate(theta)==TRUE)


  #################################################################
  ##                    generate observations                    ##
  #################################################################
  # generate hidden variables s set
  s <- sample(1:k, size = n, replace = TRUE, prob = p)

  # store multi-dimensional matrix, as mvrnorm does not support indexing
  x <- matrix(0, nrow=n, ncol=dimension_gaussian)

  # generate observed variable set, with each observation drawn from a
  # multivariate GMM indexed by vector s
  for (i in 1:n)
    x[i,] <- MASS::mvrnorm(n = 1, mu=mu[,s[i]], Sigma=sigma[,,s[i]],
                     tol = 1e-12, empirical = FALSE)

  return(list(k = k, p = p, mu = mu, sigma = sigma, x = x, s = s))
}


#'  Return initial estimates to the EM algorithm for GMM estimation
#'
#'  One of the main drawback of the EM algorithm is that it requires initial guess as starting point.
#'  And so, careful initialisation, depending on the properties of the mixture, is required:
#'  * `initialize_em_univariate` returns the initial estimates in the univariate dimension.
#'  * `initialize_em_multivariate`returns the initial estimates in a multivariate context. It's worth noting
#'  that *quantiles* initialisation method is not available in the multivariate context, as
#'  no unique set of parametrisation could be returned.
#'
#' @author Bastien CHASSAGNOL
#'
#' @param x the vector of the observations
#' @param k the number of components
#' @param nstart the number of random restarts with kmeans, random and small EM method
#' @param prior_prob minimal uncertainity added to the minor components of each observation assigned by hierarchical clustering
#' @param short_iter,short_eps hyperparameters of the small EM method
#' @param initialisation_algorithm the choice of the initialisation method, between kmeans, quantiles, random, hc, small em and rebmix method
#' @param ... additional hyperparameters supplied with some of the initialisation methods
#'
#' @return a list of the estimated parameters, ordered by increasing mean for identifiability issues
#'
#' @export

initialize_em_univariate <- function(x = NULL, k = 2, nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                          initialisation_algorithm = c("kmeans", "quantiles", "random", "hc", "small em", "rebmix"), ...) {
  if (is.null(x)) {
    stop("Please provide data for x.")
  }

  initialisation_algorithm <- match.arg(
    initialisation_algorithm,
    c("kmeans", "quantiles", "random", "hc", "small em", "rebmix")
  )
  n <- length(x)

  ### initialize parametrization
  if (initialisation_algorithm == "kmeans") {
    fit <- stats::kmeans(x = x, centers = k, nstart = nstart, iter.max = short_iter)
    estimated_theta <- list(
      p = fit$size / n,
      mu = as.vector(fit$centers),
      sigma = sqrt(fit$withinss / fit$size)
    )
  } else if (initialisation_algorithm == "quantiles") {
    fit <- bgmm::init.model.params(X = x, k = k, method = "all")
    estimated_theta <- list(
      p = fit$pi,
      mu = as.vector(fit$mu),
      sigma = sqrt(as.vector(fit$cvar))
    )
  } else if (initialisation_algorithm == "random") {
    all_logs <- lapply(1:nstart, function(y) {
      fit <- EMCluster::simple.init(as.matrix(x), nclass = k)
      logLikelihood_per_random <- EMCluster::logL(as.matrix(x), fit)
      return(list(parameters = fit, logLikelihood = logLikelihood_per_random))
    })

    best_model <- which.max(all_logs %>% purrr::map_dbl("logLikelihood"))
    fit <- purrr::map(all_logs, "parameters")[[best_model]]
    estimated_theta <- list(
      p = fit$pi,
      mu = as.vector(fit$Mu),
      sigma = sqrt(as.vector(fit$LTSigma))
    )
  } else if (initialisation_algorithm == "hc") {
    clustering_hc <- mclust::hcV(data = x, minclus = 1, ...)
    partition_hc <- as.vector(mclust::hclass(clustering_hc, G = k))
    # generate posterior probability, adding artificially min probability to be assigned to each cluster
    z <- mclust::unmap(partition_hc)
    z[z == 1] <- 1 - prior_prob * (k - 1)
    z[z == 0] <- prior_prob
    # generate first iteration, starting with m step
    fit <- mclust::meV(
      data = x, z = z, prior = NULL,
      control = mclust::emControl(itmax = 1, equalPro = FALSE)
    )
    estimated_theta <- list(p = fit$parameters$pro, mu = fit$parameters$mean, sigma = sqrt(fit$parameters$variance$sigmasq))
  } else if (initialisation_algorithm == "small em") {
    all_logs <- lapply(1:nstart, function(y) {
      # take some random points using EMCluster simple method
      start <- EMCluster::simple.init(as.matrix(x), nclass = k)
      # small runs of EM from mixtools package, as implementing the true EM algorithm
      fit <- em_Rmixmod_univariate(
        x = x,
        start = list(p = start$pi, mu = as.vector(start$Mu), sigma = sqrt(as.vector(start$LTSigma))),
        short_eps = short_eps, short_iter = short_iter, k = k
      )

      logLikelihood_per_em <- EMCluster::logL(as.matrix(x),
        emobj = list(pi = fit$p, Mu = as.matrix(fit$mu), LTSigma = as.matrix(fit$sigma^2))
      )
      return(list(parameters = fit, logLikelihood = logLikelihood_per_em))
    })

    best_model <- which.max(all_logs %>% purrr::map_dbl("logLikelihood"))
    estimated_theta <- purrr::map(all_logs, "parameters")[[best_model]]
  } else if (initialisation_algorithm == "rebmix") {
    EM_control <- methods::new("EM.Control",
      strategy = "exhaustive", variant = "EM", acceleration = "fixed",
      acceleration.multiplier = 1, tolerance = short_eps, maximum.iterations = 1
    )

    fit <- suppressMessages(rebmix::REBMIX(
      object = "REBMIX", Dataset = list(data.frame(x)), Preprocessing = "kernel density estimation",
      Restraints = "loose", cmax = k, cmin = k, Criterion = "BIC", pdf = "normal", model = "REBMIX", EMcontrol = EM_control
    ))

    estimated_theta <- list(p = unlist(fit@w), mu = unlist(fit@Theta[[1]])[seq(2, 3 * k, by = 3)], sigma = unlist(fit@Theta[[1]])[seq(3, 3 * k, by = 3)])
  }

  ####### order list by mean values, to deal with identifiability issues
  ordered_estimated_theta <- list(
    p = estimated_theta$p[order(estimated_theta$mu)],
    mu = sort(estimated_theta$mu),
    sigma = estimated_theta$sigma[order(estimated_theta$mu)]
  )
  # ensure that probabilities sum up to one
  ordered_estimated_theta$p <- ordered_estimated_theta$p / sum(ordered_estimated_theta$p)
  ordered_estimated_theta$p[k] <- 1 - sum(ordered_estimated_theta$p[-k])
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(function(x) as.numeric(unname(x)))
  return(ordered_estimated_theta)
}


#' @rdname initialize_em_univariate
initialize_em_multivariate <- function(x, k = 2, nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
                                     initialisation_algorithm = c("kmeans", "random", "hc", "small em", "rebmix"), ...) {

  initialisation_algorithm <- match.arg(
    initialisation_algorithm,
    c("kmeans", "random", "hc", "small em", "rebmix")
  )
  n <- nrow(x) # number of observations
  dim_gaussian <- ncol(x) # dimension of the observation space

  ### initialize parametrization
  if (initialisation_algorithm == "kmeans") {
    fit <- stats::kmeans(x = x, centers = k, nstart = nstart, iter.max = short_iter)
    estimated_theta <- estimate_supervised_multivariate_GMM(x=x,s = fit$cluster, k = k)

  }
  else if (initialisation_algorithm == "random") {
    all_logs <- lapply(1:nstart, function(y) {
      fit <- EMCluster::simple.init(x, nclass = k)
      logLikelihood_per_random <- EMCluster::logL(x, fit)
      return(list(parameters = fit, logLikelihood = logLikelihood_per_random))
    })

    best_model <- which.max(all_logs %>% purrr::map_dbl("logLikelihood"))
    fit <- purrr::map(all_logs, "parameters")[[best_model]]
    estimated_theta <- list(p = fit$pi, mu = t(fit$Mu),
                            sigma = trig_mat_to_array(fit$LTSigma))
  }
  else if (initialisation_algorithm == "hc") {
    clustering_hc <- mclust::hcVVV(data = x, minclus = 1, ...)
    estimated_theta <- estimate_supervised_multivariate_GMM(x=x,
                                                            s = c(mclust::hclass(clustering_hc, G = k)), k = k)
  }
  else if (initialisation_algorithm == "small em") {
    all_logs <- lapply(1:nstart, function(y) {
      # take some random points using EMCluster simple method
      start <- EMCluster::simple.init(as.matrix(x), nclass = k)
      # small runs of EM from mixtools package, as implementing the true EM algorithm
      fit <- em_Rmixmod_multivariate(
        x = x,
        start = list(p = start$pi, mu = t(start$Mu), sigma = trig_mat_to_array(start$LTSigma)),
        short_eps = short_eps, short_iter = short_iter, k = k
      )

      logLikelihood_per_em <- EMCluster::logL(x, pi = fit$p, Mu = t(fit$mu), LTSigma = fit$sigma)
      return(list(parameters = fit, logLikelihood = logLikelihood_per_em))
    })

    best_model <- which.max(all_logs %>% purrr::map_dbl("logLikelihood"))
    estimated_theta <- purrr::map(all_logs, "parameters")[[best_model]]
  }
  else if (initialisation_algorithm == "rebmix") {
    EM_control <- methods::new("EM.Control",
                               strategy = "exhaustive", variant = "EM", acceleration = "fixed",
                               acceleration.multiplier = 1, tolerance = short_eps, maximum.iterations = 1
    )

    # one step of EM is required, just for the computation of the estimated parameters
    fit <- suppressMessages(rebmix::REBMIX(
      model = "REBMVNORM", Dataset = list(data.frame(x)), Preprocessing = "kernel density estimation",
      Restraints = "loose", cmax = k , cmin = k - 1, Criterion = "BIC", EMcontrol = EM_control
    ))

    mu <- matrix(fit@Theta[[1]][grep("^theta1\\.",names(fit@Theta[[1]]))] %>% unlist() %>% unname(),
            nrow=dim_gaussian, ncol=k) # estimates are returned as an unique long vector

    sigma <- array(fit@Theta[[1]][seq(3, 3 * k, by = 3)] %>% unlist() %>% unname(),
                   dim=c(dim_gaussian, dim_gaussian, k))
    estimated_theta <- list(p = unlist(fit@w), mu = mu, sigma = sigma)
  }

  ####### order list by mean values, to deal with identifiability issues
  ordered_components <- do.call(order, estimated_theta$mu %>% as.data.frame())
  ordered_estimated_theta <- list(
    p = estimated_theta$p[ordered_components],
    mu = estimated_theta$mu[, ordered_components],
    sigma = estimated_theta$sigma[,,ordered_components]
  )

  # ensure that probabilities sum up to one
  ordered_estimated_theta$p <- ordered_estimated_theta$p / sum(ordered_estimated_theta$p)
  ordered_estimated_theta$p[k] <- 1 - sum(ordered_estimated_theta$p[-k])
  ordered_estimated_theta$p <- ordered_estimated_theta$p %>% as.numeric() %>% unname()
  return(ordered_estimated_theta)
}





#'  Custom R implementation of the EM algorithm
#'
#' @author Bastien CHASSAGNOL
#'
#' @param x the vector of the observations
#' @param k the number of components
#' @param itmax the maximal number of iterations to reach the threshold
#' @param epsilon the criterion threshold considered as the tolerance between two consecutive log-likelihoods
#' @param start list of initial estimates provided by the user
#' @param initialisation_algorithm,nstart hyper-parameters, when the user rather uses one of our implemented initialization algorithms
#' @param parallel only relevant for GMKMCharlie package which has a native parallel implementation (by default, takes half of the available clusters)
#' @param prior the mclust object used to store the supposed prior distributions of the parameters' components (only relevant if a Bayesian implementation is required)
#' @param ... additional parameters for the reviewed packages
#'
#' @return a list of the estimated parameters, ordered by increasing mean for identifiability issues
#'
#' @importFrom stats IQR dnorm median quantile runif sd var
#'
#' @export




emnmix <- function(x, k, itmax = 5000, epsilon = 10^-12, nstart = 10L, start = NULL,
                   initialisation_algorithm = "kmeans", ...) {

  ### retrieve the initial configuration
  x <- as.vector(as.matrix(x))
  n <- length(x) # get number of observations

  # get initial estimates if not provided by the user
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }
  p <- start$p
  mu <- start$mu
  sigma <- start$sigma
  old_loglik <- +Inf

  for (iter in 1:itmax) {
    ###  E-step  ###
    # compute the posterior probabilities
    expection_results <- predict_posterior_probability(x, list(p = p, mu = mu, sigma = sigma))
    eta <- expection_results$eta
    new_loglik <- expection_results$loglik

    ### M-step ###

    # intermediate calculations
    S0 <- apply(eta, 2, sum)
    S1 <- apply(eta * x, 2, sum)
    S2 <- apply(eta * x^2, 2, sum)

    # update parameter vector theta
    p <- S0 / n
    mu <- S1 / S0
    sigma <- sqrt((S2 - 2 * mu * S1 + mu^2 * S0) / S0)

    # deal with underflow or removal of components
    if (!check_parameters_validity_univariate(list(p = p, mu = mu, sigma = sigma), k = k)) {
      break
    } ## early stop in case the algorithm is trapped in boundary space

    else if (abs(new_loglik - old_loglik) < epsilon) { # early stop if criterion threshold is met
      break
    }


    old_loglik <- new_loglik # update the newly computed log-likelihood
  }

  # return parameters by mean values for identifiability issues
  ordered_estimated_theta <- list(
    p = p[order(mu)],
    mu = sort(mu),
    sigma = sigma[order(mu)]
  )
  return(ordered_estimated_theta)
}



#' @describeIn emnmix EM implementation with Rmixmod package
#' @importClassesFrom Rmixmod Strategy GaussianParameter
#' @export
em_Rmixmod_univariate <- function(x = x, k = 2, initialisation_algorithm = "kmeans",
                       itmax = 5000, epsilon = 10^-12, start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  # define strategy
  strategy_Rmixmod <- methods::new("Strategy",
    algo = "EM", nbTryInInit = 10L,
    initMethod = "parameter", nbIterationInAlgo = itmax,
    epsilonInAlgo = epsilon,
    parameter = methods::new("GaussianParameter",
      proportions = start$p,
      mean = as.matrix(start$mu),
      variance = lapply(start$sigma, function(x) as.matrix(x^2))
    )
  )

  # fit the model
  fit <- Rmixmod::mixmodCluster(
    data = as.vector(x), nbCluster = k, dataType = "quantitative",
    strategy = strategy_Rmixmod
  )@bestResult@parameters

  fit <- list(
    p = fit@proportions,
    mu = as.vector(fit@mean),
    sigma = sqrt(as.vector(unlist(fit@variance)))
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$p[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = fit$sigma[order(fit$mu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}


#' @rdname em_Rmixmod_univariate  EM implementation with Rmixmod package in the mulvariate context
#' @importClassesFrom Rmixmod Strategy GaussianParameter
#' @export


em_Rmixmod_multivariate <- function(x = x, k = 2, initialisation_algorithm = "kmeans",
                                  itmax = 5000, epsilon = 10^-12, start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  # define strategy
  strategy_Rmixmod <- methods::new("Strategy",
                                   algo = "EM", nbTryInInit = 10L,
                                   initMethod = "parameter", nbIterationInAlgo = itmax,
                                   epsilonInAlgo = epsilon,
                                   parameter = methods::new("GaussianParameter",
                                                            proportions = start$p,
                                                            mean = start$mu,
                                                            variance = lapply(seq(dim(start$sigma)[3]), function(x) start$sigma[ , , x])
                                   )
  )

  # fit the model
  fit <- Rmixmod::mixmodCluster(
    data = x %>% as.data.frame(), nbCluster = k, dataType = "quantitative",
    strategy = strategy_Rmixmod
  )@bestResult@parameters

  fit <- list(
    p = fit@proportions,
    mu = fit@mean,
    sigma = fit@variance %>% simplify2array()
  )

  # return an ordered list by mean values
  ordered_components <- do.call(order, fit$mu %>% as.data.frame())
  ordered_estimated_theta <- list(
    p = fit$p[ordered_components],
    mu = fit$mu[, ordered_components],
    sigma = fit$sigma[,,ordered_components]
  )

  return(ordered_estimated_theta)
}





#' @describeIn emnmix EM implementation with EMCluster package
#' @export
em_EMCluster <- function(x = x, k = 2, initialisation_algorithm = "kmeans",
                         itmax = 5000, epsilon = 10^-12, start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }


  # fit the model
  fit <- EMCluster::emcluster(
    x = as.matrix(x),
    pi = start$p, Mu = as.matrix(start$mu), LTSigma = as.matrix(start$sigma),
    EMC = EMCluster::.EMControl(em.iter = itmax, em.eps = epsilon)
  )

  fit <- list(
    p = fit$pi,
    mu = as.vector(fit$Mu),
    sigma = sqrt(as.vector(fit$LTSigma))
  )


  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$p[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = fit$sigma[order(fit$mu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}



#' @describeIn emnmix EM implementation with bgmm package
#' @export
em_bgmm <- function(x = x, k = 2, itmax = 5000, epsilon = 10^-12,
                    initialisation_algorithm = "kmeans", start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  # fit the model
  fit <- bgmm::unsupervised(
    X = x, k = k,
    init.params = list(
      pi = start$p,
      mu = as.matrix(start$mu),
      cvar = array(start$sigma^2, dim = c(k, 1, 1))
    ),
    stop.likelihood.change = epsilon, stop.max.nsteps = itmax, ...
  )

  fit <- list(
    p = fit$pi,
    mu = as.vector(fit$mu),
    sigma = sqrt(as.vector(fit$cvar))
  )


  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$p[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = fit$sigma[order(fit$mu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}


#' @describeIn emnmix EM implementation with flexmix package
#' @export
em_flexmix <- function(x = x, k = 2, itmax = 5000, epsilon = 10^-12,
                       initialisation_algorithm = "kmeans", start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }
  # predicted class from initial guess
  predicted_classes_original <- predict_posterior_probability(x, start)$eta

  # fit the model
  fit <- flexmix::flexmix(x ~ 1,
    k = k, cluster = predicted_classes_original,
    model = flexmix::FLXMCnorm1(),
    control = list(minprior = 0, tolerance = epsilon, iter.max = itmax, verbose = 0)
  )

  fit <- list(
    p = fit@size / length(x),
    mu = flexmix::parameters(fit)["mean", ],
    sigma = flexmix::parameters(fit)["sd", ]
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$p[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = fit$sigma[order(fit$mu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}


#' @describeIn emnmix EM implementation with mixtools package
#' @export
em_mixtools <- function(x = x, k = 2, initialisation_algorithm = "hc",
                        itmax = 5000, epsilon = 10^-12, start = NULL, ...) {
  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  fit <- mixtools::normalmixEM(
    x = x, lambda = start$p, mu = start$mu, sigma =
      start$sigma, k = k, epsilon = epsilon, maxit = itmax
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$lambda[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = fit$sigma[order(fit$mu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}

#' @describeIn emnmix EM implementation with mclust package
#' @export
em_mclust <- function(x = x, k = 2, initialisation_algorithm = "hc", start = NULL,
                      itmax = 5000, epsilon = 10^-12, prior = NULL, ...) {

  # set relevant parameters to be done
  control <- mclust::emControl(tol = epsilon, itmax = itmax)
  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  start_em <- list(pro = start$p, mean = start$mu, variance = list(sigmasq = start$sigma^2))
  estimated_theta <- mclust::emV(
    data = x, parameters = start_em,
    modelNames = "V", warn = FALSE, control = control, prior = prior
  )$parameters
  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = estimated_theta$pro[order(estimated_theta$mean)],
    mu = sort(estimated_theta$mean),
    sigma = sqrt(estimated_theta$variance$sigmasq[order(estimated_theta$mean)])
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}


#' @describeIn emnmix EM implementation with DCEM package
#' @export
em_DCEM <- function(x = x, k = 2, initialisation_algorithm = "hc",
                    itmax = 5000, epsilon = 10^-12, start = NULL, ...) {
  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  fit <- DCEM:::dcem_cluster_uv(as.matrix(x),
    meu = start$mu, sigma = start$sigma, prior = start$p,
    num_clusters = k, iteration_count = itmax, threshold = epsilon, num_data = length(x), numcols = 1
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$prior[order(fit$meu)],
    mu = sort(fit$meu),
    sigma = fit$sigma[order(fit$meu)]
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}

# em algorithm for mixture models using GMKMcharlie package
#' @describeIn emnmix EM implementation with GMKMcharlie package
#' @export
em_GMKMcharlie_univariate <- function(x = x, k = 2, initialisation_algorithm = "hc",
                           itmax = 5000, epsilon = 10^-12, start = NULL, parallel = FALSE, ...) {
  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  avalaible_processors <- ifelse(parallel, parallel::detectCores() %/% 2, 1)
  fit <- GMKMcharlie::GM(t(as.matrix(x)),
    alpha = start$p, mu = t(as.matrix(start$mu)), sigma = t(as.matrix(start$sigma)), G = k,
    convergenceEPS = epsilon, maxIter = itmax, maxCore = avalaible_processors, verbose = FALSE
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$alpha[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = sqrt(unlist(fit$sigma[order(fit$mu)]))
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}


# em algorithm for multivariate mixture models using GMKMcharlie package
#' @describeIn emnmix EM implementation with GMKMcharlie package
#' @export
em_GMKMcharlie_multivariate <- function(x = x, k = 2, initialisation_algorithm = "hc",
                                      itmax = 5000, epsilon = 10^-12, start = NULL, parallel = FALSE, ...) {
  # initialization section
  if (is.null(start)) {
    start <- initialize_em_univariate(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  avalaible_processors <- ifelse(parallel, parallel::detectCores() %/% 2, 1)
  fit <- GMKMcharlie::GM(t(as.matrix(x)),
                         alpha = start$p, mu = t(as.matrix(start$mu)), sigma = t(as.matrix(start$sigma)), G = k,
                         convergenceEPS = epsilon, maxIter = itmax, maxCore = avalaible_processors, verbose = FALSE
  )

  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$alpha[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = sqrt(unlist(fit$sigma[order(fit$mu)]))
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}





#' @describeIn emnmix EM implementation with otrimle package (designed to deal with outliers especially)
#' @export
em_otrimle <- function(x = x, k = 2, initialisation_algorithm = "hc",
                       itmax = 5000, epsilon = 10^-12, ...) {
  if (initialisation_algorithm == "hc") {
    start <- otrimle::InitClust(as.matrix(x), G = 2, modelName = "V")
  } else {
    stop("This initialisation method is not enabled with otrimle package.")
  }
  avalaible_processors <- parallel::detectCores() %/% 2
  fit <- otrimle::otrimle(
    data = as.matrix(x), G = k, initial = start,
    iter.max = itmax, tol = epsilon, ncores = avalaible_processors, monitor = FALSE
  )

  # return an ordered list by mean values + proportion of outliers
  ordered_estimated_theta <- list(
    p = fit$pi[-1][order(fit$mean)],
    mu = sort(fit$mean),
    sigma = sqrt(fit$cov[order(fit$mean)]),
    prop_extra = fit$exproportion
  )
  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}






#' Estimate the parameters in supervised case (all labels associated to the observations are available)
#'
#' @author Bastien CHASSAGNOL
#'
#' @param x the univariate or multivariate distribution of observed variables
#' @param s the index of the clusters, ranging from 1 to \eqn{k}, the number of clusters
#' @param s_outliers binary index, with 0 corresponding to outlying points
#' @param k the number of classes (by default, the number of unique values within the indicator vector)
#'
#' @return a list with three arguments, respectively the proportions, the means, and the variance / covariance
#'
#' @seealso  \code{\link{simulate_univariate_GMM}}
#'
#' @export

estimate_supervised_univariate_GMM <- function(x, s, s_outliers=s) {

  ##### retrieve relevant values
  x <- simulated_distribution$x
  s_outliers <- simulated_distribution$s_outliers
  k <- simulated_distribution$k

  # remove outliers
  x <- x[s_outliers != 0]
  s_outliers <- s_outliers[s_outliers != 0]

  parameters_colnames <- c(
    names(unlist(simulated_distribution[c("p", "mu", "sigma")])),
    paste0("mean", 1:k), paste0("sd", 1:k)
  )
  ##### retrieve the parameters
  parameters_per_component <- lapply(split(x, f = as.factor(s_outliers)), function(x_subset) {
    # fit model for each component
    model_fitted <- sn::selm(x_subset ~ 1, family = "SN", fixed.param = skewed)
    coefs_dp <- sn::coef(model_fitted, param.type = "dp")
    coefs_cp <- sn::coef(model_fitted, param.type = "cp")
    return(list(
      mu = coefs_dp["mean"], sigma = coefs_dp["sd"],
      mean = coefs_cp["mean"], sd = coefs_cp["s.d."]
    ))
  })

  observed_estimated_theta <- stats::setNames(c(
    as.vector(table(s_outliers)) / length(s_outliers),
    parameters_per_component %>% purrr::map_dbl("mean"),
    parameters_per_component %>% purrr::map_dbl("sd")
  ), nm = parameters_colnames[1:(3 * k)])
  return(observed_estimated_theta)
}


#' @rdname estimate_supervised_univariate_GMM
#' @export

estimate_supervised_multivariate_GMM <- function(x, s, k=2) {
  k <- length(unique(s)); dim_gaussian <- ncol(x); n <- nrow(x)
  mu <- matrix(0, nrow = dim_gaussian, ncol=k); sigma <- array(0, dim=c(dim_gaussian, dim_gaussian, k))

  # estimate the ratios
  p <- table(s)/n

  for (j in 1:k) {
    mu[,j] <- apply(x[s==j,], 2, mean) # estimate the means
    sigma[,,j] <- cov_MLE(x[s==j, ]) # estimate covariances
  }

  return(list(p=p, mu=mu, sigma=sigma))
}




#' @title Helper functions for the parameter estimation of GMMs
#'
#' @description These functions are small chunks of code designed to decompose the computation of the EM algorithm into simpler steps.
#' `logsumexp` returns the computation of equation \eqn{\log(\exp(sum(x)))}, avoiding numerical overflows
#'
#' @param l a vector of numeric terms
#' @return a numeric scalar value, result of the previously described equation

logsumexp <- function(l) {
  i <- which.max(l)
  res <- l[i] + log1p(sum(exp(l[-i] - l[i])))
  if (is.nan(res)) res <- -Inf
  return(res)
}


#' @describeIn logsumexp `predict_posterior_probability` returns the expected probability for each observation
#' to belong to any of the \eqn{k} clusters set a priori, given the estimated parameters
#'
#' @param x the vector of observed values, of size \eqn{n}
#' @param estimated_theta the estimated parameters
#' @return a list with two elements:
#' * the posterior probability matrix, `eta`: \eqn{\eta=(\eta_{i,j}) \in [0, 1]^{n \times k}}, with \eqn{\eta_{i,j}}
#' giving the posterior probability of observation \eqn{i} to belong to cluster \eqn{j}
#' * `loglik` returns the expected log-likelihood of our experiment

predict_posterior_probability <- function(x, estimated_theta) {
  # get relevant parameters
  x <- as.vector(x)
  k <- length(estimated_theta$p)
  n <- length(x)
  p <- estimated_theta$p
  mu <- estimated_theta$mu
  sigma <- estimated_theta$sigma
  eta <- matrix(NA, nrow = n, ncol = k) # store posterior distribution for each observation (P(S=i|X))

  # compute posterior probability, using log tip
  for (j in 1:k) {
    eta[, j] <- log(p[j]) + dnorm(x, mu[j], sigma[j], log = TRUE)
  }
  aux <- apply(eta, 1, logsumexp)
  eta <- exp(eta - aux)
  return(list(eta = eta, loglik = sum(aux)))
}


#' Check whether the estimation has been trapped in the boundary space
#'
#' * Function `check_parameters_validity_univariate` asserts at each step of the EM algorithm that
#' it doesn't fall in a degenerate case (either the package performing the EM computation has failed, and returns an
#' error message, or the algorithm is trapped in the boundary space, leading to inconsistent division by zero).
#' * Function `check_parameters_validity_multivariate` has the same functionality, but
#' is adjusted to multivariate parametrisation, and includes additionally a checking whether
#' the covariance matrix is positive definite or not
#'
#' @seealso [logsumexp()]
#' @export

check_parameters_validity_univariate <- function(theta, k = 2) {
  machine_limit <- .Machine$double.eps
  machine_max <- .Machine$double.xmax
  theta <- theta[c("p", "mu", "sigma")]
  if (any(is.na(unlist(theta))) | length(unlist(theta)) == 0) {
    warning(paste0(
      "Missing data with estimated theta, with ratios: ", paste(theta$p, collapse = " / "), ", ",
      "mu: ", paste(theta$mu, collapse = " / "), " and sigma: ", paste(theta$sigma, collapse = " / ")
    ))
    return(FALSE)
  } else if (any(theta$p < machine_limit | theta$p > 1 - machine_limit) |
    any(theta$sigma < machine_limit) | any(theta$sigma > machine_max) | any(sapply(theta, length) != k)) {
    warning(paste0(
      "Numerical overflow with estimated theta, with ratios: ", paste(theta$p, collapse = " / "), ", ",
      "mu: ", paste(theta$mu, collapse = " / "), " and sigma: ", paste(theta$sigma, collapse = " / ")
    ))
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' @rdname check_parameters_validity_univariate
check_parameters_validity_multivariate <- function(theta, k = length(theta$p)) {
  machine_limit <- .Machine$double.eps; machine_max <- .Machine$double.xmax
  is_valid_parametrisation <- TRUE # store whether the parametrisation is correctly performed

  p <- theta$p; mu <- theta$mu; sigma <- theta$sigma; k <- length(p) # get values from theta parameter
  dimension_gaussian <- dim(mu)[1]

  if (any(is.na(unlist(theta))) | length(unlist(theta)) == 0) {
    warning("NA values or no output from the estimation"); is_valid_parametrisation <- FALSE
  }

  # check the parametrisation of proportions (sum-to-one constraint)
  if (sum(p)!=1 | any(p<machine_limit) | any(p>1 - machine_limit)) {
    warning("One at least of your proportions does not enforce the sum-to-one constraint"); is_valid_parametrisation <- FALSE
  }

  # check parametrisation of means
  if (!all.equal(dim(mu),c(dimension_gaussian, k)) | any(c(mu) > machine_max)) {
    warning("The mean parameter must be of dimension ndim * k, with k the number of clusters"); is_valid_parametrisation <- FALSE
  }

  # check the parametrisation of covariance
  if (!all.equal(dim(sigma), c(dimension_gaussian, dimension_gaussian, k)) | any(c(sigma) > machine_max)) {
    warning("the covariance array stores for each component a variance matrix of fixed size k"); is_valid_parametrisation <- FALSE
  }
  else {
    for (j in 1:k) {
      if (!is_positive_definite(sigma[,,j])) {
        warning(glue::glue("Covariance matrix for component {j} is not positive definite")); is_valid_parametrisation <- FALSE
        break
      }
    }
  }

  return(is_valid_parametrisation)
}










