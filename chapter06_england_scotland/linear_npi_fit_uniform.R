
# =============================================================================
# linear_npi_fit_uniform.R
#
# Chapter 6: England and Scotland application
# Joint fit with a PIECEWISE-LINEAR NPI ramp, UNIFORM (flat) priors,
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu).
#
# Uniform priors are specified via parameter bounds in the Stan data block.
# This is the sensitivity companion to linear_npi_fit.R (hierarchical priors).
#
# Data: chapter06_england_scotland/data/GB_data.csv
#
# Outputs:
#   outputs/stan_fits/fit_het_lin_unif.rds     fitted object
#   
#   outputs/figures/fit_het_lin_unif_fitted.png       fitted vs observed
#   outputs/figures/fit_het_lin_unif_residuals.png    residuals vs time
#   outputs/figures/fit_het_lin_unif_trace.png        MCMC trace plots
#   outputs/figures/fit_het_lin_unif_densities.png    posterior densities
#  
#
# The fitted object is named fit_het_lin_unif. Required by
# compare_npi_specifications_uniform.R. Do not rename it.
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

library(rstan)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (!dir.exists("outputs/figures"))   dir.create("outputs/figures",   recursive = TRUE)
if (!dir.exists("outputs/stan_fits")) dir.create("outputs/stan_fits", recursive = TRUE)


# ----------------------------
# 0) Stan code: JOINT LINEAR (HET) — common CV_shared, IFR fixed, I0 estimated
# ----------------------------
stan_code_joint_lin_het <- "
functions {
  vector delay_probs(int L, real lmu, real lsd) {
    vector[L] p;
    for (l in 1:L) {
      real F_hi = lognormal_cdf(l, lmu, lsd);
      real F_lo = (l == 1) ? 0 : lognormal_cdf(l-1, lmu, lsd);
      p[l] = fmax(F_hi - F_lo, 1e-12);
    }
    p /= sum(p);
    return p;
  }
}

data {
  int<lower=1> C;
  int<lower=2> T;
  int<lower=0> y[C, T];
  vector<lower=0>[C] N;

  real<lower=0> rho;
  real<lower=0> delta;
  real<lower=0> gamma;

  int<lower=1> L;
  real lmu_id_fix;
  real<lower=0> lsd_id_fix;

  int<lower=1, upper=T> t1[C];
  real<lower=0, upper=1> IFR_fixed;

  // Uniform bounds (same style as your single-country template)
  real<lower=0> R0_min;    real<lower=R0_min> R0_max;
  real<lower=0, upper=1> cstar_min;  real<lower=cstar_min, upper=1> cstar_max;
  real<lower=0> w_min;     real<lower=w_min>  w_max;

  real<lower=0> CV_min;    real<lower=CV_min> CV_max;
  real<lower=0> phi_min;   real<lower=phi_min> phi_max;

  real<lower=0> I0_min;    real<lower=I0_min> I0_max;

  // Fixed multiplier E0 = m * I0 (m fixed at 2.5 in your thesis)
  real<lower=1> E0_multiplier_min;
  real<upper=10> E0_multiplier_max;
  real<lower=E0_multiplier_min, upper=E0_multiplier_max> E0_multiplier;
}

parameters {
  vector<lower=R0_min, upper=R0_max>[C] R0;
  vector<lower=cstar_min, upper=cstar_max>[C] cstar;
  vector<lower=w_min, upper=w_max>[C] w;

  real<lower=CV_min, upper=CV_max> CV_shared;

  vector<lower=phi_min, upper=phi_max>[C] phi;

  vector<lower=I0_min, upper=I0_max>[C] I0;
}

transformed parameters {
  real CV2 = square(CV_shared);
  vector[L] p_delay = delay_probs(L, lmu_id_fix, lsd_id_fix);

  vector[C] t0;
  vector[C] E0 = E0_multiplier * I0;

  matrix[C, T] mu;
  matrix<lower=0>[C, T] inc_EI;
  inc_EI = rep_matrix(0.0, C, T);

  for (c in 1:C) t0[c] = t1[c] - w[c];

  for (c in 1:C) {
    real beta = R0[c] / (rho / delta + 1.0 / gamma);

    vector[T] S;
    vector[T] E;
    vector[T] I;
    vector[T] R;
    vector[T] inc_new;

    S[1] = fmax(N[c] - E0[c] - I0[c], 1.0);
    E[1] = E0[c];
    I[1] = I0[c];
    R[1] = 0.0;

    inc_new = rep_vector(0.0, T);
    inc_EI[c, 1] = delta * E[1];

    for (t in 1:(T-1)) {
      real c_t;

      if (t < t0[c]) {
        c_t = 1.0;
      } else if (t < t1[c]) {
        real denom = fmax(t1[c] - t0[c], 1e-6);
        real frac  = (t - t0[c]) / denom;
        c_t = 1.0 - (1.0 - cstar[c]) * fmin(fmax(frac, 0.0), 1.0);
      } else {
        c_t = cstar[c];
      }

      real lambda = c_t * beta * (rho * E[t] + I[t]);
      real sus_pw = pow(S[t] / N[c], 1.0 + CV2);

      real dS = - lambda * sus_pw;
      real dE =   lambda * sus_pw - delta * E[t];
      real dI =   delta * E[t] - gamma * I[t];
      real dR =   (1.0 - IFR_fixed) * gamma * I[t];

      S[t+1] = fmax(S[t] + dS, 1e-9);
      E[t+1] = fmax(E[t] + dE, 0.0);
      I[t+1] = fmax(I[t] + dI, 0.0);
      R[t+1] = fmax(R[t] + dR, 0.0);

      inc_new[t]     = fmax(-dS, 0.0);
      inc_EI[c, t+1] = delta * E[t+1];
    }

    // infections -> deaths (convolution)
    for (t in 1:T) {
      real m = 0.0;
      int mlag = (t-1 < L) ? (t-1) : L;
      if (mlag > 0)
        for (l in 1:mlag) m += inc_new[t - l] * p_delay[l];
      mu[c, t] = IFR_fixed * fmax(m, 1e-12);
    }
  }
}

model {
  // implicit uniform via bounds
  
  // Priors
  R0 ~ uniform(R0_min,R0_max);
  w ~ uniform(w_min, w_max);
 CV_shared ~uniform (CV_min,CV_max);
  cstar ~ uniform(cstar_min, cstar_max);
  phi ~ uniform(phi_min, phi_max);
  I0 ~ uniform(I0_min, I0_max);
  
  for (c in 1:C)
    y[c] ~ neg_binomial_2(mu[c] + 1e-12, phi[c]);
}

generated quantities {
  int y_rep[C, T];
  matrix[C, T] log_lik;

  // store the piecewise-linear contact modifier c_t for each country and time
  matrix[C, T] ct;

  for (c in 1:C)
    for (t in 1:T) {
      y_rep[c, t] = neg_binomial_2_rng(mu[c, t] + 1e-12, phi[c]);
      log_lik[c, t] = neg_binomial_2_lpmf(y[c, t] | mu[c, t] + 1e-12, phi[c]);
    }

  // recompute c_t on the same piecewise-linear schedule used in the SEIR loop
  for (c in 1:C) {
    for (t in 1:T) {
      if (t < t0[c]) {
        ct[c, t] = 1.0;
      } else if (t < t1[c]) {
        real denom = fmax(t1[c] - t0[c], 1e-6);
        real frac  = (t - t0[c]) / denom;
        ct[c, t] = 1.0 - (1.0 - cstar[c]) * fmin(fmax(frac, 0.0), 1.0);
      } else {
        ct[c, t] = cstar[c];
      }
    }
  }
}

"

writeLines(stan_code_joint_lin_het, "outputs/stan_fits/seir_nb_joint_linear_uniform_het_CVcommon_IFRfixed_I0est.stan")


# ----------------------------
# 2) Compile
# ----------------------------
mod_het_lin_unif <- stan_model("outputs/stan_fits/seir_nb_joint_linear_uniform_het_CVcommon_IFRfixed_I0est.stan")

# ----------------------------
# 3) Data prep 
# ----------------------------
prepare_joint_data <- function(filepath = "chapter06_england_scotland/data/GB_data.csv") {
  cat("Preparing joint data for England and Scotland...\n")
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data$Date <- dmy(data$Date)
  data <- data[order(data$Date), ]
  
  lockdown_start_EN <- as.Date("2020-03-26")
  lockdown_end_EN   <- as.Date("2020-05-12")
  lockdown_start_SC <- as.Date("2020-03-24")
  lockdown_end_SC   <- as.Date("2020-05-28")
  
  analysis_start <- as.Date("2020-01-31")
  analysis_end   <- as.Date("2020-06-01")
  idx <- which(data$Date >= analysis_start & data$Date <= analysis_end)
  
  deaths_EN <- data$Deaths_EN[idx]; deaths_EN[is.na(deaths_EN)] <- 0L
  deaths_SC <- data$Deaths_SC[idx]; deaths_SC[is.na(deaths_SC)] <- 0L
  dates     <- data$Date[idx]
  
  t1_EN <- which.min(abs(dates - lockdown_start_EN))
  t2_EN <- which.min(abs(dates - lockdown_end_EN))
  t1_SC <- which.min(abs(dates - lockdown_start_SC))
  t2_SC <- which.min(abs(dates - lockdown_end_SC))
  
  cat(sprintf("Using %s to %s; %d days\n", min(dates), max(dates), length(dates)))
  
  list(
    england = list(deaths = as.integer(deaths_EN), dates = dates, n_steps = length(dates),
                   t1 = t1_EN, t2 = t2_EN, population = 56e6),
    scotland = list(deaths = as.integer(deaths_SC), dates = dates, n_steps = length(dates),
                    t1 = t1_SC, t2 = t2_SC, population = 5.5e6)
  )
}

joint_data <- prepare_joint_data("chapter06_england_scotland/data/GB_data.csv")


joint_data <- prepare_joint_data("GB_data.csv")

# Fixed process params 
rho   <- 0.5
delta <- 1/5.5
gamma <- 1/4

# Delay distribution params 
L     <- 60L
lmu_id_fix <- 3.151
lsd_id_fix <- 0.469

# IFR fixed (your baseline choice)
IFR_fixed <- 0.009

# Bounds 
bounds <- list(
  R0_min = 2,       R0_max = 8.0,
  cstar_min = 0.05, cstar_max = 0.9,
  w_min  = 0.01,    w_max  = 60.0,
  CV_min = 0.01,    CV_max = 3.0,
  phi_min = 10.0,   phi_max = 1000.0,
  I0_min = 1,      I0_max = 10000.0,
  E0_multiplier_min = 1,
  E0_multiplier_max = 4,
  E0_multiplier     = 2.5
)

# Build joint stan_data
dates  <- joint_data$england$dates
T_all  <- joint_data$england$n_steps

y_mat <- rbind(joint_data$england$deaths, joint_data$scotland$deaths)
storage.mode(y_mat) <- "integer"

stan_data_joint <- list(
  C = 2L,
  T = as.integer(T_all),
  y = y_mat,
  N = as.numeric(c(joint_data$england$population, joint_data$scotland$population)),
  
  rho = as.numeric(rho),
  delta = as.numeric(delta),
  gamma = as.numeric(gamma),
  
  L = as.integer(L),
  lmu_id_fix = as.numeric(lmu_id_fix),
  lsd_id_fix = as.numeric(lsd_id_fix),
  
  t1 = as.integer(c(joint_data$england$t1, joint_data$scotland$t1)),
  IFR_fixed = as.numeric(IFR_fixed),
  
  R0_min = bounds$R0_min, R0_max = bounds$R0_max,
  cstar_min = bounds$cstar_min, cstar_max = bounds$cstar_max,
  w_min = bounds$w_min, w_max = bounds$w_max,
  CV_min = bounds$CV_min, CV_max = bounds$CV_max,
  phi_min = bounds$phi_min, phi_max = bounds$phi_max,
  I0_min = bounds$I0_min, I0_max = bounds$I0_max,
  
  E0_multiplier_min = bounds$E0_multiplier_min,
  E0_multiplier_max = bounds$E0_multiplier_max,
  E0_multiplier     = bounds$E0_multiplier
)


# ----------------------------
# 4) Fit
# ----------------------------
fit_one <- function(mod, stan_data, seed) {
  sampling(
    object  = mod,
    data    = stan_data,
    chains  = 4,
    iter    = 2500,
    warmup  = 1250,
    seed    = seed,
    control = list(adapt_delta = 0.98, max_treedepth = 15)
  )
}

fit_het_lin_unif <- fit_one(mod_het_lin_unif, stan_data_joint, seed = 20250120)

rstan::check_hmc_diagnostics(fit_het_lin_unif)

#saveRDS(fit_het_lin_unif, "outputs/stan_fits/fit_het_lin_unif.rds")


create_custom_theme <- function(
    title_size = 16,
    subtitle_size = 14,
    axis_title_size = 14,
    axis_text_size = 13,
    legend_title_size = 13,
    legend_text_size = 12,
    legend_position = "right",
    grid_color_major = "grey90",
    grid_color_minor = "grey95",
    base_size = 12,
    base_family = "",
    remove_minor_grid = TRUE
) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = subtitle_size, hjust = 0.5),
      
      axis.title      = element_text(size = axis_title_size),
      axis.text       = element_text(size = axis_text_size),
      
      legend.title    = element_text(size = legend_title_size),
      legend.text     = element_text(size = legend_text_size),
      legend.position = legend_position,
      
      panel.grid.major = element_line(colour = grid_color_major),
      panel.grid.minor = if (remove_minor_grid) element_blank() else element_line(colour = grid_color_minor),
      
      # helps plots look “finished” (still minimal)
      axis.line        = element_line(colour = "grey60", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      
      strip.text       = element_text(size = axis_text_size, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

# =========================================================
# 6) Choose model for plotting
# =========================================================
fit_use     <- fit_het_lin_unif
model_label <- "Joint Linear ramp:uniform"
include_CV  <- TRUE

post <- rstan::extract(fit_use)

country_colors <- c("England" = "#1f77b4", "Scotland" = "#ff7f0e")

# =========================================================
# 7) Fitted trajectories
# =========================================================
n_days <- length(dates)
y_EN <- as.numeric(joint_data$england$deaths)
y_SC <- as.numeric(joint_data$scotland$deaths)

mu_EN <- apply(post$mu[,1,], 2, median)
mu_SC <- apply(post$mu[,2,], 2, median)

phi_EN <- median(post$phi[,1])
phi_SC <- median(post$phi[,2])

band_en <- list(
  lo = qnbinom(0.025, size = phi_EN, mu = mu_EN),
  hi = qnbinom(0.975, size = phi_EN, mu = mu_EN)
)
band_sc <- list(
  lo = qnbinom(0.025, size = phi_SC, mu = mu_SC),
  hi = qnbinom(0.975, size = phi_SC, mu = mu_SC)
)

fitted_data <- data.frame(
  Date = rep(dates, 2),
  Country = rep(c("England", "Scotland"), each = n_days),
  Observed = c(y_EN, y_SC),
  Fitted = c(mu_EN, mu_SC),
  Lower = c(band_en$lo, band_sc$lo),
  Upper = c(band_en$hi, band_sc$hi)
)

p1_fitted <- ggplot(fitted_data, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Country), alpha = 0.3) +
  geom_line(aes(y = Fitted, color = Country), linewidth = 1.2) +
  geom_point(aes(y = Observed, color = Country), size = 1, alpha = 0.7) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors) +
  
  create_custom_theme() +
  theme(legend.position = "bottom") +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p1_fitted)
ggsave("outputs/figures/fit_het_lin_unif_fitted.png", p1_fitted, width = 10, height = 7, dpi = 300)

# =========================================================
# 8) Residuals
# =========================================================
residuals_EN <- y_EN - mu_EN
residuals_SC <- y_SC - mu_SC
rmse_EN <- sqrt(mean(residuals_EN^2))
rmse_SC <- sqrt(mean(residuals_SC^2))

cat(sprintf("\n=== RMSE (%s) ===\n", model_label))
cat(sprintf("England RMSE: %.2f\n", rmse_EN))
cat(sprintf("Scotland RMSE: %.2f\n", rmse_SC))

resid_data <- data.frame(
  Date = rep(dates, 2),
  Country = rep(c("England","Scotland"), each = n_days),
  Residual = c(residuals_EN, residuals_SC),
  Fitted = c(mu_EN, mu_SC)
)

p2_resid <- ggplot(resid_data, aes(x = Date, y = Residual, color = Country)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = country_colors) +
  labs(title = paste0("Residuals vs time - ", model_label),
       subtitle = paste0("RMSE: EN=", round(rmse_EN,2), ", SC=", round(rmse_SC,2)),
       x = "Date", y = "Residual") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y")

print(p2_resid)
ggsave("outputs/figures/fit_het_lin_unif_residuals.png", p2_resid, width = 10, height = 6, dpi = 300)

# =========================================================
# 9) Trace plots
# =========================================================
params_to_trace <- c(
  "R0[1]","R0[2]",
  "cstar[1]","cstar[2]",
  "w[1]","w[2]",
  "t0[1]","t0[2]",
  "I0[1]","I0[2]",
  "phi[1]","phi[2]"
)
if (include_CV && ("CV_shared" %in% names(post))) params_to_trace <- c(params_to_trace, "CV_shared")

arr <- as.array(fit_use, pars = params_to_trace)
it <- dim(arr)[1]; ch <- dim(arr)[2]; pn <- dim(arr)[3]
par_names <- dimnames(arr)[[3]]

trace_df <- data.frame(
  Iteration = rep(1:it, times = ch*pn),
  Chain = factor(rep(rep(1:ch, each = it), times = pn)),
  Parameter = factor(rep(par_names, each = it*ch), levels = par_names),
  Value = as.vector(arr)
) %>%
  mutate(Parameter = ifelse(as.character(Parameter) == "CV_shared", "CV", as.character(Parameter))) %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

p_trace <- ggplot(trace_df, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.75, linewidth = 0.35) +
  labs(title = paste0("Trace plots - ", model_label),
       x = "Iteration", y = "Value") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)
ggsave("outputs/figures/fit_het_lin_unif_trace.png", p_trace, width = 12, height = 9, dpi = 300)

# =========================================================
# 10) Posterior densities
# =========================================================
R0_EN <- as.numeric(post$R0[,1]);  R0_SC <- as.numeric(post$R0[,2])
cs_EN <- as.numeric(post$cstar[,1]); cs_SC <- as.numeric(post$cstar[,2])
I0_EN <- as.numeric(post$I0[,1]);   I0_SC <- as.numeric(post$I0[,2])
phi_EN_d <- as.numeric(post$phi[,1]); phi_SC_d <- as.numeric(post$phi[,2])
w_EN  <- as.numeric(post$w[,1]);    w_SC  <- as.numeric(post$w[,2])
t0_EN <- as.numeric(post$t0[,1]);   t0_SC <- as.numeric(post$t0[,2])
CV_draw <- if (include_CV && ("CV_shared" %in% names(post))) as.numeric(post$CV_shared) else NULL

draw_list <- list(
  "R0[1]"    = R0_EN,    "R0[2]"    = R0_SC,
  "cstar[1]" = cs_EN,    "cstar[2]" = cs_SC,
  "w[1]"     = w_EN,     "w[2]"     = w_SC,
  "t0[1]"    = t0_EN,    "t0[2]"    = t0_SC,
  "I0[1]"    = I0_EN,    "I0[2]"    = I0_SC,
  "phi[1]"   = phi_EN_d, "phi[2]"   = phi_SC_d
)
if (include_CV && !is.null(CV_draw)) draw_list[["CV"]] <- CV_draw

post_dist <- bind_rows(lapply(names(draw_list), function(nm) {
  tibble(Parameter = nm, Value = as.numeric(draw_list[[nm]]))
})) %>%
  mutate(Parameter = factor(Parameter, levels = names(draw_list)))

meds <- post_dist %>% group_by(Parameter) %>% summarise(med = median(Value), .groups="drop")

p_dens <- ggplot(post_dist, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.6) +
  geom_vline(data = meds, aes(xintercept = med), linetype = "dashed", linewidth = 1) +
  labs(title = paste0("Posterior densities - ", model_label),
       subtitle = "Dashed = posterior median",
       x = "Value", y = "Density") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_dens)
ggsave("outputs/figures/fit_het_lin_unif_densities.png", p_dens, width = 12, height = 9, dpi = 300)

