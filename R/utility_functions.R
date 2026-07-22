# ===============================================================================
# utility_functions.R
# 
# 
# This file contains all mathematical and utility functions needed for the
# simultaneous inference of susceptibility distributions and intervention
# effects from epidemic curves analysis using the reduced model with parameter v.
#
# Functions include:
# - Mathematical transformations (logit, expit)
# - ODE system definition for reduced SEIR model with heterogeneity parameter ν
# - Data simulation function for reduced model
# 
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# Date: August 2025
# ===============================================================================

# Load required libraries (called by scripts that source this file)
library(deSolve)     # For solving differential equations
library(tidyverse)   # For data manipulation and plotting

# ===============================================================================
# SECTION 1: MATHEMATICAL UTILITY FUNCTIONS
# ===============================================================================

#' Logit transformation
#' 
#' Transforms a probability p to the logit scale. This is useful for ensuring
#' parameters stay between 0 and 1 during optimization.
#' 
#' @param p A probability between 0 and 1
#' @return The logit of p: log(p/(1-p))
#' @examples
#' logit(0.5)  # Returns 0
#' logit(0.8)  # Returns approximately 1.39
logit <- function(p) {
  log(p/(1-p))
}

#' Inverse logit (expit) transformation
#' 
#' Transforms a value on the real line back to a probability between 0 and 1.
#' This is the inverse of the logit function.
#' 
#' @param x A real number
#' @return The inverse logit of x: 1/(1+exp(-x))
#' @examples
#' expit(0)    # Returns 0.5
#' expit(1.39) # Returns approximately 0.8
expit <- function(x) {
  1/(1+exp(-x))
}

# ===============================================================================
# SECTION 2: ODE SYSTEM DEFINITION FOR REDUCED MODEL
# ===============================================================================

#' ODE function for reduced SEIR model with intervention and heterogeneity parameter ν
#' 
#' This is the main ODE system for the reduced SEIR model with time-varying 
#' transmission due to non-pharmaceutical interventions (NPIs) and heterogeneity
#' in susceptibility characterized by parameter ν.
#' 
#' The model includes:
#' - S: Susceptible individuals
#' - E: Exposed (infected but not yet infectious) individuals  
#' - I: Infectious individuals
#' - R: Recovered individuals
#' - C: Cumulative cases (for data fitting)
#' 
#' Key features:
#' - Heterogeneity parameter ν: when ν=0 (homogeneous), when ν>0 (heterogeneous)
#' - Time-varying transmission: c(t) function models intervention effects
#' - Force of infection includes both E and I compartments with relative infectiousness ρ
#' 
#' @param t Current time
#' @param y State variables [S, E, I, R, C]  
#' @param parms Named vector of parameters including:
#'   - R0: Basic reproduction number
#'   - gamma: Recovery rate (1/infectious period)
#'   - rho: Relative infectiousness in E compartment
#'   - delta: Rate of transition from E to I (1/incubation period)
#'   - v: Heterogeneity parameter (coefficient of variation)
#'   - N: Total population size
#'   - t0, t1, t2, t3: Intervention timing parameters
#'   - c_value1, c_value2, c_value3: Intervention strength parameters
#' @return List containing derivatives [dS/dt, dE/dt, dI/dt, dR/dt, dC/dt]

#' ODE function for homogeneous SEIR model with time-varying transmission
#' 
#' @param time Current time point
#' @param y Vector of state variables
#' @param parms List of model parameters
#' @return List containing the derivatives of all state variables
seir.ct <- function(time, y, parms) {
  with(as.list(parms), {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]
    
    # Define intervention function
    prox <- 1.0  # Default value
    
    if (time <= t0) {
      prox <- c_value1
    } else if (time <= t1) {
      prox <- c_value1 + (c_value2 - c_value1) * (time - t0) / (t1 - t0)
    } else if (time <= t2) {
      prox <- c_value2
    } else if (time <= t3) {
      prox <- c_value3 + (c_value3 - c_value2) * (time - t3) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    # SEIR model equations
    dS <- -Beta*prox*(rho*E+I)*(S/N)^(1+v^2)    
    dE <- Beta*prox*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI <- delta*E-gamma*I
    dR <- gamma*I
    dC <- delta*E
    return(list(c(dS, dE, dI, dR, dC)))
  })
}

#' ODE function for reduced model with intervention (homogeneous when v=0)
#' 
#' @param t Current time point
#' @param y Vector of state variables
#' @param parms List of model parameters
#' @return List containing the derivatives of all state variables
Reduced.m_intervene <- function(t, y, parms) {
  with(as.list(c(t, y, parms)), {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]
    
    # Define intervention function
    prox <- 1.0  # Default value
    
    if (t <= t0) {
      prox <- c_value1
    } else if (t <= t1) {
      prox <- c_value1 - (c_value1 - c_value2) * (t - t0) / (t1 - t0)
    } else if (t <= t2) {
      prox <- c_value2
    } else if (t <= t3) {
      prox <- c_value3 + (c_value3 - c_value2) * (t - t3) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    # Calculate transmission rate with intervention effect
    
    Beta <- R0*prox/(rho / delta + 1 / gamma)
    
    
    # SEIR model equations
    dS<- -Beta*(rho*E+I)*(S/N)^(1+v^2)    
    dE<- Beta*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI<- delta*E-gamma*I
    dR<- gamma*I
    dC<-delta*E
    return(list(c(dS, dE, dI, dR,dC)))
  })
  
}
# ===============================================================================
# section 3: Data simulation functions
# ===============================================================================

#' Simulate cases using a reduced SEIR model with intervention
#' 
#' This function simulates epidemic data using a reduced SEIR model with
#' intervention (non-pharmaceutical interventions). The model can be homogeneous
#' (v=0) or heterogeneous (v>0).
#' 
#' @param R0 Basic reproduction number
#' @param delta Rate of leaving exposed compartment
#' @param rho Relative infectiousness in E compartment
#' @param gamma Recovery rate
#' @param v Coefficient of variation of susceptibility (0 for homogeneous)
#' @param N Total population size
#' @param E0 Initial number of exposed individuals
#' @param I0 Initial number of infectious individuals
#' @param t0 Time when adaptive behavior begins
#' @param t1 Time when lockdown begins
#' @param t2 Time when lockdown ends
#' @param t3 Time when transmission returns to baseline
#' @param c_value1 Initial transmission factor
#' @param c_value2 Intervention strength
#' @param c_value3 Final transmission factor
#' @param tfinal Final simulation time
#' @return List containing simulated data and plots
simulate_cases_reduced_model <- function(R0, delta, rho, gamma, v, N, E0, I0, 
                                         t0, t1, t2, t3, c_value1, c_value2, c_value3, tfinal) {
  
  params <- list(R0, delta, rho, gamma, v, N, E0, I0, t0, t1, t2, t3, c_value1, c_value2, c_value3, tfinal)
  names(params) <- c("R0", "delta", "rho", "gamma", "v", "N", "E0", "I0", "t0", 
                     "t1", "t2", "t3", "c_value1", "c_value2", "c_value3", "tfinal")
  
  parms <- c(R0, gamma, rho, delta, v, N, t0, t1, t2, t3, c_value1, c_value2, c_value3)
  names(parms) <- c("R0", "gamma", "rho", "delta", "v", "N", "t0", "t1", "t2", "t3", 
                    "c_value1", "c_value2", "c_value3")
  
  # Initial state
  initial_state <- c(S = N-E0-I0, E = E0, I = I0, R = 0, C = 0)
  
  # Time vector for simulation
  times <- seq(0, params$tfinal, by = 1)
  
  # Perform numerical integration
  result <- as.data.frame(ode(y = initial_state, times = times, func = Reduced.m_intervene, parms = parms))
  Daily_incidence <- c(0, diff(result[, "C"]))
  df <- result %>% mutate(Inc = Daily_incidence)
  
  # Generate Poisson-distributed cases
  cases <- rpois(length(df[, "Inc"]), lambda = df$Inc)
  
  # Create data frame of simulated data
  sim_data <- data.frame(time = times, reports = cases)
  
  
  
  return(list(sim_data = sim_data))
}

# ===============================================================================
# END OF UTILITY FUNCTIONS
# ===============================================================================
#
# Summary of functions provided:
# 1. logit() - Transform probability to logit scale
# 2. expit() - Transform from logit scale back to probability  
# 3. Reduced.m_intervene() - ODE system for reduced SEIR model with parameter v
# 4. simulate_cases_reduced_model() - Simulate epidemic data with Poisson noise
#
# These functions support both homogeneous (v=0) and heterogeneous (v>0) 
# epidemic modeling with time-varying transmission due to interventions.
#
# ===============================================================================
