# ===============================================================================
# mle_functions.R
#
# Maximum likelihood estimation functions for the heterogeneous reduced SEIR model
# with non-pharmaceutical interventions.
#
# This file is the canonical, merged replacement for:
#   - MLE_functions_paper.R  (older version; now retired)
#   - Maxlik_functions_reduced_model.R  (clean version; base for this file)
#
# What is NOT here (lives in utility_functions.R instead):
#   - logit() and expit()
#   - Reduced.m_intervene()  (the ODE)
#   - simulate_cases_reduced_model()
#
# Functions provided:
#   Section 1: Core likelihood
#     poisson.loglik.NPI.reduced()        main Poisson log-likelihood (times via arg or global)
#     poisson.loglik.withNPI.reduced()    alternative (derives times from sim.data$time)
#   Section 2: Single epidemic MLE
#     f4_optim_reducedm.NPI()             objective function (returns neg-loglik)
#     fit4_reducedm_loglik.NPI()          fits model to one epidemic
#   Section 3: Two epidemics MLE
#     f4_2epi_optim_reduced.NPI()         objective function, two epidemics
#     fit4_2epic_reduced.loglik.NPI()     fits model to two epidemics simultaneously
#   Section 4: Profile likelihood
#     profile_likelihood_reducedm()       single epidemic
#     two_epi_profile_obj_fn()            helper wrapper
#     profile_likelihood_two_epidemics()  two epidemics
#   Section 5: Trajectory generation
#     generate_trajectory()               run ODE forward for given parameters
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# Date: August 2025
# ===============================================================================

library(deSolve)
library(tidyverse)

# ===============================================================================
# SECTION 1: CORE LIKELIHOOD FUNCTIONS
# ===============================================================================

#' Poisson log-likelihood for the reduced SEIR model with NPI
#'
#' Integrates the reduced SEIR ODE (Reduced.m_intervene from utility_functions.R)
#' and evaluates the Poisson log-likelihood against observed daily incidence.
#'
#' @param params  Named vector of model parameters (R0, v, t0, c_value2, etc.)
#' @param sim.data  Data frame with columns 'time' and 'reports'
#' @param initial_state  Named vector [S, E, I, R, C]
#' @param times  Time vector for integration; if NULL uses global variable 'times'
#' @return  Scalar log-likelihood value
poisson.loglik.NPI.reduced <- function(params, sim.data, initial_state, times = NULL) {

  if (is.null(times)) {
    if (!exists("times", envir = .GlobalEnv)) {
      stop("'times' must be passed as an argument or defined in the global environment")
    }
    times <- get("times", envir = .GlobalEnv)
  }

  out <- as.data.frame(ode(
    y     = initial_state,
    times = times,
    func  = Reduced.m_intervene,
    parms = params
  ))

  Daily_incidence <- c(0, diff(out[, "C"]))
  df <- out %>% mutate(Inc = Daily_incidence)

  lambda_ <- pmax(df[, "Inc"], 0.0001)

  if (!is.data.frame(sim.data)) sim.data <- as.data.frame(sim.data)
  if (!"reports" %in% names(sim.data)) stop("'reports' column missing from sim.data")
  if (nrow(sim.data) != nrow(df))      stop("sim.data and model output have different lengths")

  loglik <- sum(dpois(x = sim.data[, "reports"], lambda = lambda_, log = TRUE))
  return(loglik)
}


#' Alternative Poisson log-likelihood that derives the time vector from sim.data$time
#'
#' Handles the case where data begin at t = 1 (not t = 0) by prepending t = 0
#' for integration and trimming it before computing likelihood.
#'
#' @param params  Named vector of model parameters
#' @param sim.data  Data frame with columns 'time' and 'reports'
#' @param initial_state  Named vector [S, E, I, R, C]
#' @return  Scalar log-likelihood value
poisson.loglik.withNPI.reduced <- function(params, sim.data, initial_state) {

  times_data <- sim.data$time

  if (min(times_data) > 0) {
    times_integration <- c(0, times_data)
    needs_trimming <- TRUE
  } else {
    times_integration <- times_data
    needs_trimming <- FALSE
  }

  out <- as.data.frame(ode(
    y     = initial_state,
    times = times_integration,
    func  = Reduced.m_intervene,
    parms = params
  ))

  if (needs_trimming) out <- out[-1, ]

  Daily_incidence <- diff(c(0, out[, "C"]))
  df <- out %>% mutate(Inc = Daily_incidence)

  lambda_ <- ifelse(df[, "Inc"] <= 0, 0.0001, df[, "Inc"])

  if (!is.data.frame(sim.data)) sim.data <- as.data.frame(sim.data)
  if (!"reports" %in% names(sim.data)) stop("'reports' column missing from sim.data")
  if (nrow(sim.data) != length(lambda_)) stop("sim.data and model output have different lengths")

  loglik <- sum(dpois(x = sim.data[, "reports"], lambda = lambda_, log = TRUE))
  return(loglik)
}


# ===============================================================================
# SECTION 2: SINGLE EPIDEMIC MLE
# ===============================================================================

#' Objective function for single-epidemic MLE (returns negative log-likelihood)
#'
#' Parameter transformations on the optimisation scale:
#'   par[1] = log(R0)        -> R0 = exp(par[1])
#'   par[2] = log(v)         -> v  = exp(par[2])
#'   par[3] = log(t0)        -> t0 = exp(par[3])
#'   par[4] = logit(c_value2)-> c_value2 = expit(par[4])
#'
#' @param par      Length-4 vector of transformed parameters
#' @param sim.data Data frame with 'time' and 'reports'
#' @return  Negative log-likelihood (for minimisation)
f4_optim_reducedm.NPI <- function(par, sim.data) {

  params <- c(
    R0       = exp(par[1]),
    v        = exp(par[2]),
    t0       = exp(par[3]),
    t1       = t1_spec,
    t2       = t2_spec,
    t3       = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[4]),
    c_value3 = c_value3_spec,
    rho      = rho_spec,
    delta    = delta_spec,
    gamma    = gamma_spec,
    N        = N,
    tfinal   = tfinal_spec
  )

  loglik <- poisson.loglik.NPI.reduced(
    params        = params,
    sim.data      = sim.data,
    initial_state = initial_state
  )

  return(-loglik)
}


#' Fit the heterogeneous reduced SEIR model to a single epidemic
#'
#' Estimates R0, v (CV of susceptibility), t0 (NPI onset), and c_value2
#' (NPI contact-reduction strength) by Nelder-Mead maximisation of the
#' Poisson log-likelihood.
#'
#' @param dat  Data frame with columns 'time' and 'reports'
#' @return  List with elements:
#'   parms         natural-scale estimates plus AIC, negloglik, loglik, convergence
#'   trans_parms   transformed-scale estimates
#'   trans_hessian Hessian matrix (for Wald confidence intervals)
fit4_reducedm_loglik.NPI <- function(dat) {

  start_par <- c(log(2), log(2), log(12), logit(0.4))

  fit1 <- optim(
    par     = start_par,
    fn      = f4_optim_reducedm.NPI,
    sim.data = dat,
    method  = "Nelder-Mead",
    control = list(trace = 0, maxit = 1000),
    hessian = TRUE
  )

  fittedparams <- c(
    R0          = exp(fit1$par[1]),
    v           = exp(fit1$par[2]),
    t0          = exp(fit1$par[3]),
    c_value2    = expit(fit1$par[4]),
    AIC         = 2 * length(fit1$par) + 2 * fit1$value,
    negloglik   = fit1$value,
    loglik      = -fit1$value,
    convergence = fit1$convergence
  )

  return(list(
    parms        = fittedparams,
    trans_parms  = fit1$par,
    trans_hessian = fit1$hessian
  ))
}


# ===============================================================================
# SECTION 3: TWO EPIDEMICS MLE
# ===============================================================================

#' Objective function for two-epidemic MLE
#'
#' Sums log-likelihoods across two independent epidemics sharing the same
#' parameters (R0, v, t0, c_value2) but with different initial conditions.
#'
#' @param par        Length-4 vector of transformed parameters
#' @param sim.data_1 First epidemic data frame
#' @param sim.data_2 Second epidemic data frame
#' @return  Combined negative log-likelihood
f4_2epi_optim_reduced.NPI <- function(par, sim.data_1, sim.data_2) {

  params <- c(
    R0       = exp(par[1]),
    v        = exp(par[2]),
    t0       = exp(par[3]),
    t1       = t1_spec,
    t2       = t2_spec,
    t3       = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[4]),
    c_value3 = c_value3_spec,
    rho      = rho_spec,
    delta    = delta_spec,
    gamma    = gamma_spec,
    N        = N,
    tfinal   = tfinal_spec
  )

  loglik1 <- poisson.loglik.NPI.reduced(
    params        = params,
    sim.data      = sim.data_1,
    initial_state = initial_state_1
  )

  loglik2 <- poisson.loglik.NPI.reduced(
    params        = params,
    sim.data      = sim.data_2,
    initial_state = initial_state_2
  )

  return(-(loglik1 + loglik2))
}


#' Fit the heterogeneous reduced SEIR model to two concurrent epidemics
#'
#' Simultaneously maximises the combined Poisson log-likelihood across two
#' epidemics with different initial conditions, resolving the nu-c* compensation
#' ridge that makes single-epidemic inference unreliable.
#'
#' @param dat1  First epidemic data frame
#' @param dat2  Second epidemic data frame
#' @return  List with parms, trans_parms, trans_hessian (same structure as single-epidemic fit)
fit4_2epic_reduced.loglik.NPI <- function(dat1, dat2) {

  start_par <- c(log(2), log(2), log(12), logit(0.4))

  fit1 <- optim(
    par        = start_par,
    fn         = f4_2epi_optim_reduced.NPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method     = "Nelder-Mead",
    control    = list(trace = 0, maxit = 1500),
    hessian    = TRUE
  )

  fittedparams <- c(
    R0          = exp(fit1$par[1]),
    v           = exp(fit1$par[2]),
    t0          = exp(fit1$par[3]),
    c_value2    = expit(fit1$par[4]),
    AIC         = 2 * length(fit1$par) + 2 * fit1$value,
    negloglik   = fit1$value,
    loglik      = -fit1$value,
    convergence = fit1$convergence
  )

  return(list(
    parms        = fittedparams,
    trans_parms  = fit1$par,
    trans_hessian = fit1$hessian
  ))
}


# ===============================================================================
# SECTION 4: PROFILE LIKELIHOOD FUNCTIONS
# ===============================================================================

#' Profile likelihood for a single epidemic
#'
#' Fixes one parameter at a grid of values and re-maximises over the remaining
#' three, producing a profile likelihood curve and 95% profile confidence
#' interval (chi-squared threshold with 1 df).
#'
#' @param sim_data         Epidemic data frame
#' @param param_to_profile One of "R0", "v", "t0", "c_value2"
#' @param value_range      Grid of values; auto-computed from Hessian if NULL
#' @param n_points         Number of grid points (default 15)
#' @param plot             Whether to produce a ggplot (default TRUE)
#' @param custom_title     Optional plot title override
#' @return  List: profile_results, profile_plot, mle, ci, ci_width
profile_likelihood_reducedm <- function(sim_data, param_to_profile, value_range = NULL,
                                        n_points = 15, plot = TRUE, custom_title = NULL) {

  mle_fit        <- fit4_reducedm_loglik.NPI(dat = sim_data)
  mle_params     <- mle_fit$parms
  mle_trans_params <- mle_fit$trans_parms

  if (!is.null(mle_fit$trans_hessian)) {
    tryCatch({
      z_variance <- solve(mle_fit$trans_hessian)
      z_se       <- sqrt(diag(z_variance))

      if (param_to_profile == "R0")       param_index <- 1
      else if (param_to_profile == "v")   param_index <- 2
      else if (param_to_profile == "t0")  param_index <- 3
      else if (param_to_profile == "c_value2") param_index <- 4

      if (is.null(value_range)) {
        se          <- z_se[param_index]
        trans_lower <- mle_trans_params[param_index] - 2.5 * se
        trans_upper <- mle_trans_params[param_index] + 2.5 * se

        if (param_to_profile %in% c("R0", "v", "t0")) {
          value_range <- seq(exp(trans_lower), exp(trans_upper), length.out = n_points)
        } else if (param_to_profile == "c_value2") {
          value_range <- seq(expit(trans_lower), expit(trans_upper), length.out = n_points)
        }
      }
    }, error = function(e) {
      if (is.null(value_range)) {
        if (param_to_profile == "R0") {
          value_range <<- seq(max(0.5, mle_params["R0"] * 0.7), mle_params["R0"] * 1.3, length.out = n_points)
        } else if (param_to_profile == "v") {
          value_range <<- seq(max(0.1, mle_params["v"] * 0.7), mle_params["v"] * 1.3, length.out = n_points)
        } else if (param_to_profile == "t0") {
          value_range <<- seq(max(1, mle_params["t0"] * 0.8), mle_params["t0"] * 1.2, length.out = n_points)
        } else if (param_to_profile == "c_value2") {
          value_range <<- seq(max(0.1, mle_params["c_value2"] * 0.7), min(0.9, mle_params["c_value2"] * 1.3), length.out = n_points)
        }
      }
    })
  }

  profile_results <- data.frame(
    param_value = value_range, neg_loglik = NA,
    R0 = NA, v = NA, t0 = NA, c_value2 = NA
  )

  true_values <- list(
    R0       = if (exists("R0_spec"))      R0_spec      else mle_params["R0"],
    v        = if (exists("CV_true"))      CV_true      else mle_params["v"],
    t0       = if (exists("t0_spec"))      t0_spec      else mle_params["t0"],
    c_value2 = if (exists("c_value2_spec")) c_value2_spec else mle_params["c_value2"]
  )

  for (i in seq_along(value_range)) {
    profile_value <- value_range[i]

    if (param_to_profile == "R0") {
      profile_trans_value <- log(profile_value); param_index <- 1
    } else if (param_to_profile == "v") {
      profile_trans_value <- log(profile_value); param_index <- 2
    } else if (param_to_profile == "t0") {
      profile_trans_value <- log(profile_value); param_index <- 3
    } else if (param_to_profile == "c_value2") {
      profile_trans_value <- logit(profile_value); param_index <- 4
    }

    other_indices <- setdiff(1:4, param_index)
    start_par     <- mle_trans_params[other_indices]

    profile_obj_fn <- function(par, fixed_param_index, fixed_value) {
      full_par <- numeric(4)
      full_par[fixed_param_index] <- fixed_value
      full_par[other_indices]     <- par
      return(f4_optim_reducedm.NPI(full_par, sim_data))
    }

    opt_result <- tryCatch({
      optim(par = start_par, fn = profile_obj_fn,
            fixed_param_index = param_index, fixed_value = profile_trans_value,
            method = "Nelder-Mead", control = list(maxit = 1000))
    }, error = function(e) list(value = NA, par = rep(NA, length(start_par))))

    profile_results$neg_loglik[i] <- opt_result$value

    opt_full_par <- numeric(4)
    opt_full_par[param_index]   <- profile_trans_value
    opt_full_par[other_indices] <- opt_result$par

    profile_results$R0[i]       <- exp(opt_full_par[1])
    profile_results$v[i]        <- exp(opt_full_par[2])
    profile_results$t0[i]       <- exp(opt_full_par[3])
    profile_results$c_value2[i] <- expit(opt_full_par[4])
  }

  profile_results  <- profile_results[!is.na(profile_results$neg_loglik), ]
  min_neg_loglik   <- min(profile_results$neg_loglik)
  conf_threshold   <- qchisq(0.95, df = 1)
  profile_results$LR_stat <- 2 * (profile_results$neg_loglik - min_neg_loglik)

  ci_data  <- profile_results[profile_results$LR_stat <= conf_threshold, ]
  ci_lower <- if (nrow(ci_data) > 0) min(ci_data$param_value) else min(profile_results$param_value)
  ci_upper <- if (nrow(ci_data) > 0) max(ci_data$param_value) else max(profile_results$param_value)

  mle_x  <- mle_params[param_to_profile]
  true_x <- true_values[[param_to_profile]]

  param_labels <- list(
    R0       = expression(R[0]),
    v        = expression(nu),
    t0       = expression(t[0]),
    c_value2 = expression(c[1])
  )

  p <- NULL
  if (plot) {
    p <- ggplot(profile_results, aes(x = param_value, y = -neg_loglik)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_point(color = "darkblue", size = 2) +
      geom_hline(yintercept = -min_neg_loglik - conf_threshold / 2,
                 linetype = "dashed", color = "red", linewidth = 1.2, alpha = 0.7) +
      geom_vline(xintercept = mle_x,  color = "purple", linetype = "solid",  linewidth = 1.1) +
      geom_vline(xintercept = true_x, color = "green4", linetype = "dotted", linewidth = 1.5) +
      geom_vline(xintercept = ci_lower, color = "blue", linetype = "dotted", linewidth = 1.2) +
      geom_vline(xintercept = ci_upper, color = "blue", linetype = "dotted", linewidth = 1.2) +
      annotate("text", x = mle_x, y = max(-profile_results$neg_loglik),
               label = "MLE", color = "purple", hjust = -0.3, size = 4, fontface = "bold") +
      labs(x = param_labels[[param_to_profile]], y = "Log likelihood") +
      theme_minimal() +
      theme(
        axis.title       = element_text(size = 12, face = "bold"),
        axis.text        = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )
    print(p)
  }

  return(list(
    profile_results = profile_results,
    profile_plot    = p,
    mle             = mle_params,
    ci              = c(ci_lower, ci_upper),
    ci_width        = ci_upper - ci_lower
  ))
}


#' Helper objective function wrapper for two-epidemic profile likelihood
#'
#' @param par               Free parameters (length 3)
#' @param fixed_param_index Index of the fixed parameter (1-4)
#' @param fixed_value       Fixed value on the transformed scale
#' @param sim_data_1        First epidemic data frame
#' @param sim_data_2        Second epidemic data frame
#' @return  Negative log-likelihood with one parameter fixed
two_epi_profile_obj_fn <- function(par, fixed_param_index, fixed_value,
                                   sim_data_1, sim_data_2) {
  full_par <- numeric(4)
  full_par[fixed_param_index]            <- fixed_value
  full_par[setdiff(1:4, fixed_param_index)] <- par
  return(f4_2epi_optim_reduced.NPI(full_par, sim_data_1, sim_data_2))
}


#' Profile likelihood for two concurrent epidemics
#'
#' Same structure as profile_likelihood_reducedm() but optimises over the
#' combined two-epidemic likelihood.
#'
#' @param sim_data_large  First (larger initial condition) epidemic data frame
#' @param sim_data_small  Second (smaller initial condition) epidemic data frame
#' @param param_to_profile  One of "R0", "v", "t0", "c_value2"
#' @param value_range  Grid of values; auto-computed if NULL
#' @param n_points  Number of grid points (default 15)
#' @param plot  Whether to produce a ggplot (default FALSE)
#' @param custom_title  Optional plot title override
#' @return  List: profile_results, profile_plot, mle, ci, ci_width
profile_likelihood_two_epidemics <- function(sim_data_large, sim_data_small,
                                             param_to_profile,
                                             value_range = NULL, n_points = 15,
                                             plot = FALSE, custom_title = NULL) {

  mle_fit          <- fit4_2epic_reduced.loglik.NPI(dat1 = sim_data_large, dat2 = sim_data_small)
  mle_params       <- mle_fit$parms
  mle_trans_params <- mle_fit$trans_parms

  if (!is.null(mle_fit$trans_hessian)) {
    tryCatch({
      z_variance <- solve(mle_fit$trans_hessian)
      z_se       <- sqrt(diag(z_variance))

      if (param_to_profile == "R0")            param_index <- 1
      else if (param_to_profile == "v")        param_index <- 2
      else if (param_to_profile == "t0")       param_index <- 3
      else if (param_to_profile == "c_value2") param_index <- 4

      if (is.null(value_range)) {
        se          <- z_se[param_index]
        trans_lower <- mle_trans_params[param_index] - 2.5 * se
        trans_upper <- mle_trans_params[param_index] + 2.5 * se

        if (param_to_profile %in% c("R0", "v", "t0")) {
          value_range <- seq(exp(trans_lower), exp(trans_upper), length.out = n_points)
        } else if (param_to_profile == "c_value2") {
          value_range <- seq(expit(trans_lower), expit(trans_upper), length.out = n_points)
        }
      }
    }, error = function(e) {
      if (is.null(value_range)) {
        if (param_to_profile == "R0") {
          value_range <<- seq(max(0.5, mle_params["R0"] * 0.7), mle_params["R0"] * 1.3, length.out = n_points)
        } else if (param_to_profile == "v") {
          value_range <<- seq(max(0.1, mle_params["v"] * 0.7), mle_params["v"] * 1.3, length.out = n_points)
        } else if (param_to_profile == "t0") {
          value_range <<- seq(max(1, mle_params["t0"] * 0.8), mle_params["t0"] * 1.2, length.out = n_points)
        } else if (param_to_profile == "c_value2") {
          value_range <<- seq(max(0.1, mle_params["c_value2"] * 0.7), min(0.9, mle_params["c_value2"] * 1.3), length.out = n_points)
        }
      }
    })
  }

  profile_results <- data.frame(
    param_value = value_range, neg_loglik = NA,
    R0 = NA, v = NA, t0 = NA, c_value2 = NA
  )

  true_values <- list(
    R0       = if (exists("R0_spec"))       R0_spec       else mle_params["R0"],
    v        = if (exists("CV_true"))       CV_true       else mle_params["v"],
    t0       = if (exists("t0_spec"))       t0_spec       else mle_params["t0"],
    c_value2 = if (exists("c_value2_spec")) c_value2_spec else mle_params["c_value2"]
  )

  for (i in seq_along(value_range)) {
    profile_value <- value_range[i]

    if (param_to_profile == "R0") {
      profile_trans_value <- log(profile_value); param_index <- 1
    } else if (param_to_profile == "v") {
      profile_trans_value <- log(profile_value); param_index <- 2
    } else if (param_to_profile == "t0") {
      profile_trans_value <- log(profile_value); param_index <- 3
    } else if (param_to_profile == "c_value2") {
      profile_trans_value <- logit(profile_value); param_index <- 4
    }

    other_indices <- setdiff(1:4, param_index)
    start_par     <- mle_trans_params[other_indices]

    opt_result <- tryCatch({
      optim(par = start_par,
            fn  = two_epi_profile_obj_fn,
            fixed_param_index = param_index,
            fixed_value       = profile_trans_value,
            sim_data_1        = sim_data_large,
            sim_data_2        = sim_data_small,
            method  = "Nelder-Mead",
            control = list(maxit = 1000))
    }, error = function(e) list(value = NA, par = rep(NA, length(start_par))))

    profile_results$neg_loglik[i] <- opt_result$value

    opt_full_par <- numeric(4)
    opt_full_par[param_index]   <- profile_trans_value
    opt_full_par[other_indices] <- opt_result$par

    profile_results$R0[i]       <- exp(opt_full_par[1])
    profile_results$v[i]        <- exp(opt_full_par[2])
    profile_results$t0[i]       <- exp(opt_full_par[3])
    profile_results$c_value2[i] <- expit(opt_full_par[4])
  }

  profile_results  <- profile_results[!is.na(profile_results$neg_loglik), ]
  min_neg_loglik   <- min(profile_results$neg_loglik)
  conf_threshold   <- qchisq(0.95, df = 1)
  profile_results$LR_stat <- 2 * (profile_results$neg_loglik - min_neg_loglik)

  ci_data  <- profile_results[profile_results$LR_stat <= conf_threshold, ]
  ci_lower <- if (nrow(ci_data) > 0) min(ci_data$param_value) else min(profile_results$param_value)
  ci_upper <- if (nrow(ci_data) > 0) max(ci_data$param_value) else max(profile_results$param_value)

  mle_x  <- mle_params[param_to_profile]
  true_x <- true_values[[param_to_profile]]

  param_labels <- list(
    R0       = expression(R[0]),
    v        = expression(nu),
    t0       = expression(t[0]),
    c_value2 = expression(c[1])
  )

  p <- NULL
  if (plot) {
    p <- ggplot(profile_results, aes(x = param_value, y = -neg_loglik)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_point(color = "darkblue", size = 2) +
      geom_hline(yintercept = -min_neg_loglik - conf_threshold / 2,
                 linetype = "dashed", color = "red", linewidth = 1.2, alpha = 0.7) +
      geom_vline(xintercept = mle_x,    color = "purple", linetype = "solid",  linewidth = 1.1) +
      geom_vline(xintercept = true_x,   color = "green4", linetype = "dotted", linewidth = 1.5) +
      geom_vline(xintercept = ci_lower, color = "blue",   linetype = "dotted", linewidth = 1.2) +
      geom_vline(xintercept = ci_upper, color = "blue",   linetype = "dotted", linewidth = 1.2) +
      annotate("text", x = mle_x, y = max(-profile_results$neg_loglik),
               label = "MLE", color = "purple", hjust = -0.3, size = 4, fontface = "bold") +
      labs(x = param_labels[[param_to_profile]], y = "Log likelihood") +
      theme_minimal() +
      theme(
        axis.title       = element_text(size = 12, face = "bold"),
        axis.text        = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
      )
    print(p)
  }

  return(list(
    profile_results = profile_results,
    profile_plot    = p,
    mle             = mle_params,
    ci              = c(ci_lower, ci_upper),
    ci_width        = ci_upper - ci_lower
  ))
}


# ===============================================================================
# SECTION 5: TRAJECTORY GENERATION
# ===============================================================================

#' Run the reduced SEIR model forward and return daily incidence
#'
#' Convenience function for plotting fitted or hypothetical trajectories.
#' Initial conditions are looked up from the global environment in order:
#'   E0_1 / I0_1, then E0_fixed / I0_fixed, then E0 / I0.
#'
#' @param R0_val    Basic reproduction number
#' @param v_val     Susceptibility coefficient of variation
#' @param t0_val    NPI onset time
#' @param c_val     NPI contact-reduction strength (c_value2)
#' @param times_vec Integer time vector (t = 1, 2, ..., T)
#' @return  Numeric vector of daily incidence of length equal to times_vec
generate_trajectory <- function(R0_val, v_val, t0_val, c_val, times_vec) {

  params <- c(
    R0       = R0_val,   v       = v_val,
    rho      = rho_spec, delta   = delta_spec, gamma = gamma_spec,
    N        = N,        t0      = t0_val,
    t1       = t1_spec,  t2      = t2_spec,    t3    = t3_spec,
    c_value1 = c_value1_spec, c_value2 = c_val, c_value3 = c_value3_spec,
    tfinal   = max(times_vec)
  )

  if (exists("E0_1") && exists("I0_1")) {
    E0_use <- E0_1; I0_use <- I0_1
  } else if (exists("E0_fixed") && exists("I0_fixed")) {
    E0_use <- E0_fixed; I0_use <- I0_fixed
  } else if (exists("E0") && exists("I0")) {
    E0_use <- E0; I0_use <- I0
  } else {
    stop("Initial conditions (E0, I0) not found in the global environment")
  }

  init <- c(S = N - E0_use - I0_use, E = E0_use, I = I0_use, R = 0, C = 0)

  times_with_zero <- if (min(times_vec) > 0) c(0, times_vec) else times_vec

  out <- as.data.frame(ode(y = init, times = times_with_zero,
                           func = Reduced.m_intervene, parms = params))

  daily_incidence <- c(0, diff(out$C))

  if (min(times_vec) > 0) return(daily_incidence[-1]) else return(daily_incidence)
}


# ===============================================================================
# END
#
# Source order for analysis scripts:
#   source("R/utility_functions.R")   <- ODE, simulate_cases_reduced_model, logit, expit
#   source("R/mle_functions.R")       <- this file
# ===============================================================================
