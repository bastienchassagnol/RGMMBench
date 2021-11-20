#' Generation of a mixture drawn from an univariate GMM, with possibility to add outliers and skewness
#'
#' @author Bastien CHASSAGNOL
#'
#' @param n the number of observations to be drawn
#' @param theta a list with 4 entries, corresponding to the true values of the parameters (proportion p, mean mu, deviation sigma, skewness skew)
#' @param prop_outliers,interval the proportion of outliers in the mixture, and its range respective to the 0.05 and 0.95 quantiles of the global distribution
#'
#' @return a list with the number of components k, the true parameters p, mu, sigma, skew, the observed variables x, the hidden observations s and an indicator of the outliers s_outliers
#'
#' @export


rnmix_skewed_with_outliers <- function(n,
                                       theta = list(p = c(0.40, 0.60), mu = c(175, 165), sigma = c(10, 12), skew = c(0, 0)),
                                       prop_outliers = 0, interval = 2) {
  # get values from theta parameter
  p <- theta$p
  mu <- theta$mu
  sigma <- theta$sigma
  skew <- theta$skew
  k <- length(p)

  # generate hidden variables s set
  s <- sample(1:k, size = n, replace = TRUE, prob = p)

  # generate observed variable set
  x <- sn::rsn(n, xi = mu[s], omega = sigma[s], alpha = skew[s])

  # select randomly prop_outliers points to be drawn for uniform distribution
  # choice of points such that they are outside quantiles 0.025 of each distribution
  outliers_indexes <- sample(1:n, size = round(prop_outliers * n), replace = FALSE)
  s_outliers <- s
  s_outliers[outliers_indexes] <- 0

  # choice of enough ranged random distribution
  min_component <- which.min(sapply(1:k, function(j) sn::qsn(0.05, xi = mu[j], omega = sigma[j], alpha = skew[j])))
  max_component <- which.max(sapply(1:k, function(j) sn::qsn(0.95, xi = mu[j], omega = sigma[j], alpha = skew[j])))
  length_interval <- interval * (sn::qsn(0.95, xi = mu[max_component], omega = sigma[max_component], alpha = skew[max_component]) -
    sn::qsn(0.05, xi = mu[min_component], omega = sigma[min_component], alpha = skew[min_component]))
  x[outliers_indexes] <- runif(
    n = round(prop_outliers * n),
    min = sn::qsn(0.05, xi = mu[min_component], omega = sigma[min_component], alpha = skew[min_component]) - length_interval,
    max = sn::qsn(0.95, xi = mu[max_component], omega = sigma[max_component], alpha = skew[max_component]) + length_interval
  )


  return(list(k = k, p = p, mu = mu, sigma = sigma, skew = skew, x = x, s = s, s_outliers = s_outliers))
}


#'  Global function to launch the initialisation step of the EM algorithm
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

initialize_em <- function(x = NULL, k = 2, nstart = 10L, short_iter = 200, short_eps = 10^-2, prior_prob = 0.05,
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
  }

  else if (initialisation_algorithm == "quantiles") {
    fit <- bgmm::init.model.params(X = x, k = k, method = "all")
    estimated_theta <- list(
      p = fit$pi,
      mu = as.vector(fit$mu),
      sigma = sqrt(as.vector(fit$cvar))
    )
  }
  else if (initialisation_algorithm == "random") {
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
  }
  else if (initialisation_algorithm == "hc") {
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
  }
  else if (initialisation_algorithm == "small em") {
    all_logs <- lapply(1:nstart, function(y) {
      # take some random points using EMCluster simple method
      start <- EMCluster::simple.init(as.matrix(x), nclass = k)
      # small runs of EM from mixtools package, as implementing the true EM algorithm
      fit <- em_Rmixmod(
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
  }
  else if (initialisation_algorithm == "rebmix") {
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





#'  Custom R implementation of the EM algorithm
#'
#' @author Bastien CHASSAGNOL
#'
#' @param x the vector of the observations
#' @param k the number of components
#' @param itmax the maximal number of iterations to reach the threshold
#' @param epsilon the criterion threshold considered as the tolerance between two consecutive log-likelihoods
#' @param start list of initial estimates provided by the user
#' @param initialisation_algorithm,nstart hyperparameters, when the user rather uses one of our implemented initialisation algorithms
#' @param skew the initial guess of the user on the skewness of the distribution (only relevent for em_mixsmn function)
#' @param parallel only relevant for GMKMCharlie package which has a native parallell implementation (by default, takes half of the avalaible clusters)
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
    start <- initialize_em(
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
    if (!check_parameters_validity(list(p = p, mu = mu, sigma = sigma), k = k)) {
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
em_Rmixmod <- function(x = x, k = 2, initialisation_algorithm = "kmeans",
                       itmax = 5000, epsilon = 10^-12, start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em(
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


#' @describeIn emnmix EM implementation with EMCluster package
#' @export
em_EMCluster <- function(x = x, k = 2, initialisation_algorithm = "kmeans",
                         itmax = 5000, epsilon = 10^-12, start = NULL, ...) {

  # initialization section
  if (is.null(start)) {
    start <- initialize_em(
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
    start <- initialize_em(
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
    start <- initialize_em(
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
    start <- initialize_em(
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

# em algorithm for mixture models using mclust package
em_mclust <- function(x = x, k = 2, initialisation_algorithm = "hc", start = NULL,
                      itmax = 5000, epsilon = 10^-12, prior = NULL, ...) {

  # set relevant parameters to be done
  control <- mclust::emControl(tol = epsilon, itmax = itmax)
  # initialization section
  if (is.null(start)) {
    start <- initialize_em(
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
    start <- initialize_em(
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
em_GMKMcharlie <- function(x = x, k = 2, initialisation_algorithm = "hc",
                           itmax = 5000, epsilon = 10^-12, start = NULL, parallel = FALSE, ...) {
  # initialization section
  if (is.null(start)) {
    start <- initialize_em(
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
  }
  else {
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


#' @describeIn emnmix EM implementation with mixsmsn package (designed to deal with skewed GMMs especially)
#' @export
em_mixsmsn <- function(x = x, k = 2, initialisation_algorithm = "hc", skew = rep(0, k),
                       itmax = 5000, epsilon = 10^-12, start = NULL, ...) {
  if (is.null(start)) {
    start <- initialize_em(
      x = x, k = k, nstart = 10L, itmax = 200, epsilon = 10^-2,
      initialisation_algorithm = initialisation_algorithm, ...
    )
  }

  # nu is the kurtosis, always equal to 3 for a classic Gaussian distribution(not skewed)
  if (all(skew == 0)) {
    # fit a classic Gaussian distribution
    fit <- mixsmsn::smsn.mix(
      y = x, nu = 3,
      mu = start$mu, sigma2 = start$sigma^2, shape = skew, pii = start$p,
      g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Normal",
      error = 10^-12, iter.max = 5000, calc.im = FALSE
    )
  }
  else {
    # fit a skewed distribution
    fit <- mixsmsn::smsn.mix(
      y = x, nu = 3,
      mu = start$mu, sigma2 = start$sigma, shape = skew, pii = start$p,
      g = 2, get.init = FALSE, criteria = FALSE, group = FALSE, family = "Skew.normal",
      error = 10^-12, iter.max = 5000, calc.im = FALSE
    )
  }


  # return an ordered list by mean values
  ordered_estimated_theta <- list(
    p = fit$pii[order(fit$mu)],
    mu = sort(fit$mu),
    sigma = sqrt(fit$sigma2[order(fit$mu)])
  )

  ordered_estimated_theta <- ordered_estimated_theta %>% purrr::map(unname)
  return(ordered_estimated_theta)
}



#' Estimate the parameters in supervised case (all labels associated to the observations are available)
#'
#' @author Bastien CHASSAGNOL
#'
#' @param simulated_distribution the object returned by generating a sample drawn from a GMM
#'
#' @return observed_estimated_theta a list of the estimated parameters in supervised case
#'
#' @seealso  \code{\link{rnmix_skewed_with_outliers}}
#'
#' @export

# retrieve the parameters, when the latent variable S is observed
compute_parameters_complete_observations <- function(simulated_distribution) {

  ##### retrieve relevant values
  x <- simulated_distribution$x
  s_outliers <- simulated_distribution$s_outliers
  skew <- simulated_distribution$skew
  k <- simulated_distribution$k

  # remove outliers
  x <- x[s_outliers != 0]
  s_outliers <- s_outliers[s_outliers != 0]

  parameters_colnames <- c(
    names(unlist(simulated_distribution[c("p", "mu", "sigma")])),
    paste0("mean", 1:k), paste0("sd", 1:k)
  )

  # assumption of skewness or not for data, depending on skewness parameters
  skewed <- NULL
  if (all(skew == 0)) {
    skewed <- list(alpha = 0)
  }

  ##### retrieve the parameters

  parameters_per_component <- lapply(split(x, f = as.factor(s_outliers)), function(x_subset) {
    # fit model for each component
    model_fitted <- sn::selm(x_subset ~ 1, family = "SN", fixed.param = skewed)
    coefs_dp <- sn::coef(model_fitted, param.type = "dp")
    coefs_cp <- sn::coef(model_fitted, param.type = "cp")
    return(list(
      mu = coefs_dp["xi"], sigma = coefs_dp["omega"],
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




logsumexp <- function(l) {
  i <- which.max(l)
  res <- l[i] + log1p(sum(exp(l[-i] - l[i])))
  if (is.nan(res)) res <- -Inf
  return(res)
}


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


check_parameters_validity <- function(estimated_theta, k = 2) {
  machine_limit <- .Machine$double.eps
  machine_max <- .Machine$double.xmax
  estimated_theta <- estimated_theta[c("p", "mu", "sigma")]
  if (any(is.na(unlist(estimated_theta))) | length(unlist(estimated_theta)) == 0) {
    warning(paste0(
      "Missing data with estimated theta, with ratios: ", paste(estimated_theta$p, collapse = " / "), ", ",
      "mu: ", paste(estimated_theta$mu, collapse = " / "), " and sigma: ", paste(estimated_theta$sigma, collapse = " / ")
    ))
    return(FALSE)
  }
  else if (any(estimated_theta$p < machine_limit | estimated_theta$p > 1 - machine_limit) |
    any(estimated_theta$sigma < machine_limit) | any(estimated_theta$sigma > machine_max) | any(sapply(estimated_theta, length) != k)) {
    warning(paste0(
      "Out of computational data with estimated theta, with ratios: ", paste(estimated_theta$p, collapse = " / "), ", ",
      "mu: ", paste(estimated_theta$mu, collapse = " / "), " and sigma: ", paste(estimated_theta$sigma, collapse = " / ")
    ))
    return(FALSE)
  }
  else {
    return(TRUE)
  }
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


#  function used to compute bias and mse of a distribution
get_local_scores <- function(estimated_theta, true_theta, alpha = 0.2) {
  names_param <- names(unlist(true_theta[c("p", "mu", "sigma")]))
  mean_parameter <- apply(estimated_theta[, names_param], 2, mean) # mean of the distribution
  empiric_sds <- apply(estimated_theta[, names_param], 2, sd) # sd of the distribution
  bias <- apply(estimated_theta, 2, mean) - unlist(true_theta) # bias of the distribution
  mse <- apply(estimated_theta, 2, var) + bias^2 # mse of the distribution
  return(tibble::tibble(
    scores = c("mean", "sd", "bias", "mse"),
    dplyr::bind_rows(mean_parameter, empiric_sds, bias, mse)
  ))
}
