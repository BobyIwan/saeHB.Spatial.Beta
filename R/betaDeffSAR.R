#' @title Small Area Estimation using Hierarchical Bayesian Method under Spatial Beta SAR Model with Design Effect
#'
#' @description
#' \itemize{
#'   \item {This function is implemented to variable of interest \eqn{y} that is assumed to follow a Beta distribution. The range of data is \eqn{0 < y < 1}.}
#'   \item {This function gives estimation of small area means under Spatial Simultaneous Autoregressive (SAR) Model using Hierarchical Bayesian Method with Design Effect (DEFF) adjustment.}
#' }
#'
#' @param formula Formula that describes the fitted model.
#' @param DEFF String specifying the name of the design effect variable in the data frame.
#' @param n_i String specifying the name of the sample size variable in the data frame.
#' @param proxmat \eqn{N \times N} row-standardized proximity matrix with values in the interval \code{[0,1]} representing the spatial proximity between areas. The rows sum to \code{1}.
#' @param data The data frame.
#' @param iter.update Number of updates performed during Gibbs sampling. Default is \code{3}.
#' @param iter.mcmc Total number of MCMC iterations per chain. Default is \code{2000}.
#' @param thin Thinning rate for MCMC sampling. Must be a positive integer. Default is \code{1}.
#' @param burn.in Number of burn-in iterations discarded from each MCMC chain. Default is \code{1000}.
#' @param chains Number of parallel MCMC chains. Default is \code{2}.
#' @param n.adapt Number of iterations used for the adaptation phase in JAGS. Default is \code{1000}.
#' @param coef Optional vector containing the mean of the prior distribution of the regression model coefficients.
#' @param var.coef Optional vector containing the variances of the prior distribution of the regression model coefficients.
#' @param tau.u Initial value or shape for the random effect precision. Default is \code{1}.
#' @param seed An integer seed for the random number generator to ensure reproducibility. Default is \code{123}.
#' @param quiet Logical; if \code{TRUE}, suppresses the JAGS terminal output. Default is \code{FALSE}.
#' @param plot Logical; if \code{TRUE}, generates MCMC diagnostic trace, autocorrelation, and density plots. Default is \code{TRUE}.
#' @param keep.fit Logical; if \code{TRUE}, keeps the raw MCMC \code{coda} samples object in the output list. Default is \code{FALSE}.
#'
#' @return This function returns a list with the following objects:
#' \describe{
#'   \item{Est}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the small area means estimated using the Hierarchical Bayesian method.}
#'   \item{refVar}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the area-specific random effect variances \eqn{(a.var)}.}
#'   \item{randeff}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the area-specific random effects \eqn{(v)}.}
#'   \item{coefficient}{A dataframe containing the posterior mean estimates, posterior standard deviations, 95\% credible intervals, Rhat convergence diagnostics, and Effective Sample Sizes (ESS) for the regression coefficients \eqn{(\beta)} and the spatial autoregressive parameter \eqn{(\rho)}.}
#' }
#'
#' @examples
#' # Load dataset
#' data(dataBeta)
#' data(weight_mat)
#'
#' \donttest{
#' # Fit the Spatial Beta-SAR model with Design Effect
#' result <- betaDeffSAR(
#'   formula = y ~ x1 + x2,
#'   DEFF = "deff",
#'   n_i = "n_i",
#'   proxmat = weight_mat,
#'   data = dataBeta
#' )
#'
#' # View the estimation results
#' result$Est
#' result$refVar
#' result$randeff
#' result$coefficient
#' }
#'
#' @import rjags
#' @import coda
#' @import stats
#' @import grDevices
#' @import graphics
#'
#' @export betaDeffSAR
betaDeffSAR <- function(formula, DEFF, n_i, proxmat, data,
                        iter.update = 3, iter.mcmc = 2000,
                        thin = 1, burn.in = 1000, chains = 2, n.adapt = 1000,
                        coef = NULL, var.coef = NULL, tau.u = 1,
                        seed = 123, quiet = FALSE, plot = TRUE, keep.fit = FALSE) {

  result <- list(Est = NA, refVar = NA, randeff = NA, coefficient = NA)

  formuladata <- stats::model.frame(formula, data, na.action = NULL)
  y <- formuladata[, 1, drop = FALSE]

  if (ncol(formuladata) < 2) stop("Formula must include response and at least 1 predictor.")
  if (any(is.na(formuladata[, -1, drop = FALSE]))) stop("Auxiliary variables contain NA values.")

  y_all <- as.numeric(y[, 1])
  N     <- length(y_all)
  if (N < 2) stop("Need at least 2 areas.")
  if (any(y_all[!is.na(y_all)] <= 0) || any(y_all[!is.na(y_all)] >= 1)) {
    stop("Response variable must satisfy 0 < y < 1.")
  }

  DEFF <- data[, DEFF]
  n_i  <- data[, n_i]

  if (length(DEFF) != N) stop("Length of DEFF must equal number of areas.")
  if (length(n_i) != N) stop("Length of n_i must equal number of areas.")
  if (any(n_i[!is.na(y_all)] <= DEFF[!is.na(y_all)])) {
    stop("There is at least one sampled area where n_i <= DEFF. Effective sample size must be > 1.")
  }
  if (iter.update < 3) stop("The number of iteration updates must be at least 3.")

  xmat <- stats::model.matrix(formula, data = formuladata)
  X    <- as.matrix(xmat[, -1, drop = FALSE])
  P    <- ncol(X)
  nvar <- P + 1

  W <- as.matrix(proxmat)
  if (any(is.na(W))) stop("Proximity matrix contains NA.")
  if (nrow(W) != N || ncol(W) != N) stop("Proximity matrix must be N x N.")

  mu_beta  <- if (!is.null(coef)) coef else rep(0, nvar)
  tau_beta <- if (!is.null(var.coef)) 1/var.coef else rep(1, nvar)
  tau.ua <- 1
  tau.ub <- 1

  I <- diag(1, N)
  O <- rep(0, N)

  inits_list <- lapply(1:chains, function(c) {
    list(v = rep(0, N), beta = mu_beta, tau_u = tau.u, rho = 0,
         .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + c)
  })

  # Model definitions
  # Model 1: Fully Sampled Data (SAR)
  model_sampled <- "model {
    for (i in 1:N) {
      y[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- mu[i] * phi[i]
      shape2[i] <- (1 - mu[i]) * phi[i]
      phi[i] <- (n_i[i] / DEFF[i]) - 1

      logit(mu[i]) <- beta[1] + inprod(beta[2:(P+1)], X[i, ]) + v[i]
      a.var[i] <- sig_v[i, i]
    }

    C <- (I - rho * W)
    tau_v <- tau_u * (t(C) %*% C)
    sig_v <- inverse(tau_v)
    v ~ dmnorm(O, tau_v)

    for (k in 1:(P+1)) {
      beta[k] ~ dnorm(mu_beta[k], tau_beta[k])
    }
    tau_u ~ dgamma(tau.ua, tau.ub)
    rho ~ dunif(-0.9999, 0.9999)
    sigma2_u <- 1 / tau_u
  }"

  # Model 2: Data with NAs (SAR)
  model_nonsampled <- "model {
    for (i in 1:N_samp) {
      y_samp[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- mu[idx_samp[i]] * phi[i]
      shape2[i] <- (1 - mu[idx_samp[i]]) * phi[i]
      phi[i] <- (n_i[idx_samp[i]] / DEFF[idx_samp[i]]) - 1
    }

    for (j in 1:N) {
      logit(mu[j]) <- beta[1] + inprod(beta[2:(P+1)], X[j, ]) + v[j]
      a.var[j] <- sig_v[j, j]
    }

    C <- (I - rho * W)
    tau_v <- tau_u * (t(C) %*% C)
    sig_v <- inverse(tau_v)
    v ~ dmnorm(O, tau_v)
    for (k in 1:(P+1)) {
      beta[k] ~ dnorm(mu_beta[k], tau_beta[k])
    }
    tau_u ~ dgamma(tau.ua, tau.ub)
    rho ~ dunif(-0.9999, 0.9999)
    sigma2_u <- 1 / tau_u
  }"

  params <- c("mu", "a.var", "beta", "rho", "tau_u", "sigma2_u", "v")

  if (!any(is.na(y_all))) {
    for (i in 1:iter.update) {
      dat <- list(N = N, P = P, y = y_all, X = X, DEFF = DEFF, n_i = n_i,
                  W = W, I = I, O = O, mu_beta = mu_beta, tau_beta = tau_beta,
                  tau.ua = tau.ua, tau.ub = tau.ub)

      jags.m <- rjags::jags.model(file = textConnection(model_sampled), data = dat,
                                  inits = inits_list, n.chains = chains, n.adapt = n.adapt, quiet = quiet)

      samps <- rjags::coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin,
                                   progress.bar = if(quiet) "none" else "text")
      samps1 <- stats::window(samps, start = start(samps) + burn.in)

      res_sum = summary(samps1)
      beta_temp = res_sum$statistics[grep("^beta\\[", rownames(res_sum$statistics)), 1:2]
      for (k in 1:nvar) {
        mu_beta[k]  = beta_temp[k, 1]
        tau_beta[k] = 1/(beta_temp[k, 2]^2)
      }
      tau_idx = grep("^tau_u", rownames(res_sum$statistics))
      tau.ua = res_sum$statistics[tau_idx, 1]^2 / res_sum$statistics[tau_idx, 2]^2
      tau.ub = res_sum$statistics[tau_idx, 1] / res_sum$statistics[tau_idx, 2]^2
    }
  } else {
    idx_samp <- which(!is.na(y_all))
    N_samp   <- length(idx_samp)
    y_samp   <- y_all[idx_samp]

    for (i in 1:iter.update) {
      dat <- list(N = N, P = P, N_samp = N_samp, y_samp = y_samp, idx_samp = idx_samp,
                  X = X, DEFF = DEFF, n_i = n_i,
                  W = W, I = I, O = O,
                  mu_beta = mu_beta, tau_beta = tau_beta, tau.ua = tau.ua, tau.ub = tau.ub)

      jags.m <- rjags::jags.model(file = textConnection(model_nonsampled), data = dat,
                                  inits = inits_list, n.chains = chains, n.adapt = n.adapt, quiet = quiet)

      samps <- rjags::coda.samples(jags.m, params, n.iter = iter.mcmc, thin = thin,
                                   progress.bar = if(quiet) "none" else "text")
      samps1 <- stats::window(samps, start = start(samps) + burn.in)

      res_sum = summary(samps1)
      beta_temp = res_sum$statistics[grep("^beta\\[", rownames(res_sum$statistics)), 1:2]
      for (k in 1:nvar) {
        mu_beta[k]  = beta_temp[k, 1]
        tau_beta[k] = 1/(beta_temp[k, 2]^2)
      }
      tau_idx = grep("^tau_u", rownames(res_sum$statistics))
      tau.ua = res_sum$statistics[tau_idx, 1]^2 / res_sum$statistics[tau_idx, 2]^2
      tau.ub = res_sum$statistics[tau_idx, 1] / res_sum$statistics[tau_idx, 2]^2
    }
  }

  #Output
  res_sum <- summary(samps1)
  ESS <- coda::effectiveSize(samps1)
  if (chains > 1) {
    Rhat_raw <- coda::gelman.diag(samps1, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
  } else {
    Rhat_raw <- rep(NA, length(ESS))
  }

  mu_idx <- grep("^mu\\[", rownames(res_sum$statistics))
  Estimation <- data.frame(res_sum$statistics[mu_idx, 1:2], res_sum$quantiles[mu_idx, c(1,5)])
  colnames(Estimation) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  v_idx <- grep("^v\\[", rownames(res_sum$statistics))
  randeff <- data.frame(res_sum$statistics[v_idx, 1:2], res_sum$quantiles[v_idx, c(1,5)])
  colnames(randeff) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  a_var_idx <- grep("^a\\.var\\[", rownames(res_sum$statistics))
  refVar <- data.frame(res_sum$statistics[a_var_idx, 1:2], res_sum$quantiles[a_var_idx, c(1,5)])
  colnames(refVar) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  b_idx <- grep("^beta\\[", rownames(res_sum$statistics))
  rho_idx <- grep("^rho", rownames(res_sum$statistics))

  coef_stats <- rbind(res_sum$statistics[b_idx, 1:2], res_sum$statistics[rho_idx, 1:2])
  coef_quant <- rbind(res_sum$quantiles[b_idx, c(1,5)], res_sum$quantiles[rho_idx, c(1,5)])
  coef_rhat <- c(Rhat_raw[b_idx], Rhat_raw[rho_idx])
  coef_ess  <- c(ESS[b_idx], ESS[rho_idx])
  coefficient <- data.frame(coef_stats, coef_quant, coef_rhat, coef_ess)

  b_varnames <- character(nvar)
  for (i in 1:nvar) {
    b_varnames[i] <- paste0("beta[", i - 1, "]")
  }
  rownames(coefficient) <- c(b_varnames, "rho")
  colnames(coefficient) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "ESS")

  result$Est         <- Estimation
  result$refVar      <- refVar
  result$randeff     <- randeff
  result$coefficient <- coefficient

  if (keep.fit) result$fit <- samps1

  if (plot) {
    plot_idx <- c(b_idx, rho_idx)
    result_mcmc <- samps1[, plot_idx, drop = FALSE]
    coda::varnames(result_mcmc) <- rownames(coefficient)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mar = c(2, 2, 2, 2))

    coda::autocorr.plot(result_mcmc, col = "brown2", lwd = 2)
    plot(result_mcmc, col = "brown2", lwd = 2)
  }

  return(result)
}
