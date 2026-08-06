#' @title Small Area Estimation using Hierarchical Bayesian Method under Spatial Beta SAR Model
#'
#' @description This function estimates small area proportions using a Hierarchical Bayes (HB) method under a Spatial Simultaneous Autoregressive (SAR) Model with a Beta distribution without DEFF adjustments, by estimating the unknown precision parameter. It is designed for a variable of interest (\eqn{y}) that represents proportions, strictly bounded between 0 and 1 (\eqn{0 < y < 1}).
#'
#' @param formula An object of class \code{\link[stats]{formula}} that describes the fitted model.
#' @param proxmat An \eqn{N \times N} row-standardized spatial weights matrix (\code{style = "W"}) representing the spatial proximity between areas. The rows must sum to \code{1} and diagonal elements must be \code{0}.
#' @param data The data frame containing the variables named in \code{formula}.
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
#'   \item{est}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the small area means estimated using the Hierarchical Bayesian method.}
#'   \item{randeff}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the area-specific random effects \eqn{(v)}.}
#'   \item{refvar}{A dataframe containing the posterior mean estimates, posterior standard deviations, and 95\% credible intervals of the area-specific random effect variances \eqn{(a.var)}.}
#'   \item{coefficient}{A dataframe containing the posterior mean estimates, posterior standard deviations, 95\% credible intervals, Rhat convergence diagnostics, and Effective Sample Sizes (ESS) for the regression coefficients \eqn{(\beta)}, the spatial autoregressive parameter \eqn{(\rho)}, and the global precision parameter \eqn{(\phi)}.}
#'   \item{fit}{The raw MCMC \code{coda} object (included only if \code{keep.fit = TRUE}).}
#' }
#'
#' @examples
#' # Load dataset and proximity matrix
#' data(databeta)
#' data(weight_mat)
#'
#' \donttest{
#' # Fit the Spatial Beta-SAR model
#' result <- beta_sar(
#'   formula = y ~ x1 + x2,
#'   proxmat = weight_mat,
#'   data = databeta
#' )
#'
#' # View the estimation results
#' # 1. Small Area Estimates
#' result$est
#' # 2. Estimated area-specific random effects
#' result$randeff
#' # 3. Estimated variance of the random effects
#' result$refvar
#' # 4. Estimated regression coefficients, spatial, and precision parameters
#' result$coefficient
#' }
#'
#' @import rjags
#' @import coda
#' @import stats
#' @import grDevices
#' @import graphics
#'
#' @export beta_sar
beta_sar <- function(formula, proxmat, data,
                     iter.update = 3, iter.mcmc = 2000,
                     thin = 1, burn.in = 1000, chains = 2, n.adapt = 1000,
                     coef = NULL, var.coef = NULL, tau.u = 1,
                     seed = 123, quiet = FALSE, plot = TRUE, keep.fit = FALSE) {

  result <- list(est = NA, randeff = NA, refvar = NA, coefficient = NA)

  if (attr(terms(formula), "intercept") == 0) stop("Model must include an intercept.")

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

  if (iter.mcmc <= burn.in) stop("iter.mcmc must exceed burn.in.")
  if (thin < 1) stop("thin must be >= 1.")
  if (chains < 1) stop("chains must be >= 1.")
  if (iter.update < 3) stop("The number of iteration updates must be at least 3.")
  if (tau.u <= 0) stop("tau.u must be positive.")
  if (seed <= 0) stop("seed must be positive.")

  xmat <- stats::model.matrix(formula, data = formuladata)
  X    <- as.matrix(xmat[, -1, drop = FALSE])
  P    <- ncol(X)
  nvar <- P + 1

  if (!is.null(coef) && length(coef) != nvar) stop("coef must have length equal to the number of regression coefficients (including intercept).")
  if (!is.null(var.coef) && length(var.coef) != nvar) stop("var.coef must have length equal to the number of regression coefficients.")
  if (!is.null(var.coef) && any(var.coef <= 0)) stop("All values in var.coef must be positive.")

  W <- as.matrix(proxmat)
  if (any(is.na(W))) stop("Proximity matrix contains NA.")
  if (nrow(W) != N || ncol(W) != N) stop("Proximity matrix must be N x N.")
  if (any(abs(rowSums(W) - 1) > 1e-8)) stop("Each row of proxmat must sum to 1 (row-standardized).")
  if (any(diag(W) != 0)) stop("Diagonal elements of proxmat must be zero.")

  eig <- eigen(W)$values
  rho.min <- 1 / min(Re(eig))
  rho.max <- 1 / max(Re(eig))

  mu_beta  <- if (!is.null(coef)) coef else rep(0, nvar)
  tau_beta <- if (!is.null(var.coef)) 1/var.coef else rep(1, nvar)
  tau.ua <- 1
  tau.ub <- 1

  I <- diag(1, N)
  O <- rep(0, N)

  inits_list <- lapply(1:chains, function(c) {
    list(v = rep(0, N), beta = mu_beta, tau_u = tau.u, rho = 0, phi = 1,
         .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + c)
  })

  # Model definitions
  # Model 1: Fully Sampled Data (SAR)
  model_sampled <- "model {
    for (i in 1:N) {
      y[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- mu[i] * phi
      shape2[i] <- (1 - mu[i]) * phi

      logit(mu[i]) <- beta[1] + inprod(beta[2:(P+1)], X[i, ]) + v[i]
    }
    C <- (I - rho * W)
    tau_v <- tau_u * (t(C) %*% C)
    v ~ dmnorm(O, tau_v)

    for (k in 1:(P+1)) {
      beta[k] ~ dnorm(mu_beta[k], tau_beta[k])
    }
    tau_u ~ dgamma(tau.ua, tau.ub)
    rho ~ dunif(rho_min, rho_max)
    sigma2_u <- 1 / tau_u
    phi ~ dgamma(4, 0.2)
  }"

  # Model 2: Data with NAs (SAR)
  model_nonsampled <- "model {
    for (i in 1:N_samp) {
      y_samp[i] ~ dbeta(shape1[i], shape2[i])
      shape1[i] <- mu[idx_samp[i]] * phi
      shape2[i] <- (1 - mu[idx_samp[i]]) * phi
    }
    for (j in 1:N) {
      logit(mu[j]) <- beta[1] + inprod(beta[2:(P+1)], X[j, ]) + v[j]
    }
    C <- (I - rho * W)
    tau_v <- tau_u * (t(C) %*% C)
    v ~ dmnorm(O, tau_v)

    for (k in 1:(P+1)) {
      beta[k] ~ dnorm(mu_beta[k], tau_beta[k])
    }
    tau_u ~ dgamma(tau.ua, tau.ub)
    rho ~ dunif(rho_min, rho_max)
    sigma2_u <- 1 / tau_u
    phi ~ dgamma(4, 0.2)
  }"

  params <- c("mu", "beta", "rho", "tau_u", "sigma2_u", "v", "phi")

  if (!any(is.na(y_all))) {
    for (i in 1:iter.update) {
      dat <- list(N = N, P = P, y = y_all, X = X, W = W, I = I, O = O,
                  mu_beta = mu_beta, tau_beta = tau_beta, tau.ua = tau.ua, tau.ub = tau.ub,
                  rho_min = rho.min, rho_max = rho.max)
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
                  X = X, W = W, I = I, O = O, mu_beta = mu_beta, tau_beta = tau_beta,
                  tau.ua = tau.ua, tau.ub = tau.ub, rho_min = rho.min, rho_max = rho.max)
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

  # Output
  res_sum <- summary(samps1)

  ESS <- coda::effectiveSize(samps1)
  if (chains > 1) {
    Rhat_raw <- coda::gelman.diag(samps1, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
  } else {
    Rhat_raw <- rep(NA, length(ESS))
  }

  mu_idx <- grep("^mu\\[", rownames(res_sum$statistics))
  estimation <- data.frame(res_sum$statistics[mu_idx, 1:2], res_sum$quantiles[mu_idx, c(1,5)])
  colnames(estimation) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  v_idx <- grep("^v\\[", rownames(res_sum$statistics))
  randeff <- data.frame(res_sum$statistics[v_idx, 1:2], res_sum$quantiles[v_idx, c(1,5)])
  colnames(randeff) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  coda_mat <- as.matrix(samps1)
  tau_u_samp <- coda_mat[, "tau_u"]
  rho_samp   <- coda_mat[, "rho"]
  n_iter_samp <- nrow(coda_mat)
  a_var_mat <- matrix(NA, nrow = n_iter_samp, ncol = N)

  for(i in 1:n_iter_samp) {
    C_iter <- I - rho_samp[i] * W
    tau_v_iter <- tau_u_samp[i] * (t(C_iter) %*% C_iter)
    sig_v_iter <- solve(tau_v_iter)
    a_var_mat[i, ] <- diag(sig_v_iter)
  }

  refvar <- data.frame(
    Estimate   = colMeans(a_var_mat),
    Est.Error  = apply(a_var_mat, 2, stats::sd),
    `l-95% CI` = apply(a_var_mat, 2, stats::quantile, probs = 0.025),
    `u-95% CI` = apply(a_var_mat, 2, stats::quantile, probs = 0.975),
    check.names = FALSE
  )
  rownames(refvar) <- paste0("a.var[", 1:N, "]")
  colnames(refvar) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")

  b_idx   <- grep("^beta\\[", rownames(res_sum$statistics))
  rho_idx <- grep("^rho", rownames(res_sum$statistics))
  phi_idx <- grep("^phi", rownames(res_sum$statistics))

  coef_stats <- rbind(res_sum$statistics[b_idx, 1:2, drop = FALSE],
                      res_sum$statistics[rho_idx, 1:2, drop = FALSE],
                      res_sum$statistics[phi_idx, 1:2, drop = FALSE])
  coef_quant <- rbind(res_sum$quantiles[b_idx, c(1,5), drop = FALSE],
                      res_sum$quantiles[rho_idx, c(1,5), drop = FALSE],
                      res_sum$quantiles[phi_idx, c(1,5), drop = FALSE])
  coef_rhat  <- c(Rhat_raw[b_idx], Rhat_raw[rho_idx], Rhat_raw[phi_idx])
  coef_ess   <- c(ESS[b_idx], ESS[rho_idx], ESS[phi_idx])
  coefficient <- data.frame(coef_stats, coef_quant, coef_rhat, coef_ess)

  b_varnames <- character(nvar)
  for (i in 1:nvar) {
    b_varnames[i] <- paste0("beta[", i - 1, "]")
  }
  rownames(coefficient) <- c(b_varnames, "rho", "phi")
  colnames(coefficient) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "ESS")

  result$est         <- estimation
  result$randeff     <- randeff
  result$refvar      <- refvar
  result$coefficient <- coefficient

  if (keep.fit) result$fit <- samps1

  if (plot) {
    plot_idx <- c(b_idx, rho_idx, phi_idx)
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
