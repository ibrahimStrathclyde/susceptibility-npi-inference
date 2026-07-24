
# =============================================================================
# logistic_npi_fit.R
#
# Chapter 6: England and Scotland application
# Joint hierarchical fit with a LOGISTIC NPI profile,
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu).
#
# The contact modifier follows a smooth logistic decline:
#   c(t) = 1 - (1 - c*) * inv_logit(k * (t - tm))
# where tm is the midpoint of the decline and k its steepness.
#
# England and Scotland are fitted jointly. R0, c*, tm and k are partially
# pooled across countries; nu is shared.
#
# Data: chapter06_england_scotland/data/GB_data.csv
#       31 Jan 2020 to 01 Jun 2020 (first wave)
#
# Outputs:
#   outputs/stan_fits/seir_nb_joint_logistic_uniform_HET_HIER.stan  Stan model
#   outputs/stan_fits/fit_het_log_hier.rds                          fitted object
#   outputs/stan_fits/loo_het_log_hier.rds                          LOO object
#   outputs/figures/                                                 diagnostics
#
# The fitted object is named fit_het_log_hier. This name is required by
# compare_npi_specifications.R. Do not rename it.
#
# Sampler: 4 chains, 2500 iterations, 1250 warmup,
#          adapt_delta = 0.98, max_treedepth = 15.
# Run time: roughly 2 to 4 hours depending on hardware.
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

library(ggplot2)

create_custom_theme <- function(
    title_size = 16,
    subtitle_size = 14,
    axis_title_size = 14,
    axis_text_size = 13,
    legend_title_size = 13,
    legend_text_size = 12,
    legend_position = "right",   # e.g., "right", "bottom", or c(0.8, 0.2)
    grid_color_major = "grey90",
    grid_color_minor = "grey95",
    base_size = 11,
    base_family = ""
) {
  # Use a real base theme here (was custom_theme())
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
      panel.grid.minor = element_line(colour = grid_color_minor),
      panel.border     = element_blank(),
      strip.text       = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}


# =============================
# seir_nb_joint_hier_logistic_HET.stan
# =============================
stan_code_joint_log_het <- "
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

  real<lower=0, upper=1> IFR;

  // country-specific prior for log(I0)
  vector[C] log_I0_mu;
  vector<lower=0>[C] log_I0_sd;

  // retained for interface compatibility; unused in the heterogeneous model
  real<lower=0> CV_shared_fixed;
}

parameters {
  // population (hyper) parameters
  real  mu_log_R0;
  real<lower=0> sigma_log_R0;

  real  mu_tm;
  real<lower=0> sigma_tm;

  real  mu_log_k;
  real<lower=0> sigma_log_k;

  real  mu_logit_cstar;
  real<lower=0> sigma_logit_cstar;

  // non-centred country effects
  vector[C] z_log_R0;
  vector[C] z_tm;
  vector[C] z_log_k;
  vector[C] z_logit_cstar;

  // initial condition (country-specific)
  vector[C] log_I0;

  // heterogeneity (shared)
  real<lower=0> CV_shared;

  // NB dispersion (country-specific)
  vector<lower=1e-3>[C] phi;
}


transformed parameters {
  vector<lower=0>[C] R0 = exp(mu_log_R0 + sigma_log_R0 * z_log_R0);
  vector[C]           tm = mu_tm + sigma_tm * z_tm;
  vector<lower=0>[C]  k  = exp(mu_log_k + sigma_log_k * z_log_k);
  vector<lower=0,upper=1>[C] cstar =
    inv_logit(mu_logit_cstar + sigma_logit_cstar * z_logit_cstar);

  vector<lower=0>[C] I0 = exp(log_I0);
  vector<lower=0>[C] E0 = 2.5 * I0;

  real CV2 = square(CV_shared);

  vector[L] p_delay = delay_probs(L, lmu_id_fix, lsd_id_fix);

  matrix[C, T] mu;
  matrix<lower=0>[C, T] inc_EI;

  //  infections (S -> E), stored on time index t (1..T)
  matrix<lower=0>[C, T] inc_SE;
  inc_SE = rep_matrix(0.0, C, T);

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
      real s_t    = inv_logit(k[c] * (t - tm[c]));
      real c_t    = 1.0 - (1.0 - cstar[c]) * s_t;

      real lambda = c_t * beta * (rho * E[t] + I[t]);
      real sus_pw = pow(S[t] / N[c], 1.0 + CV2);

      real dS = -lambda * sus_pw;
      real dE =  lambda * sus_pw - delta * E[t];
      real dI =  delta * E[t] - gamma * I[t];
      real dR =  (1.0 - IFR) * gamma * I[t];

      S[t+1] = fmax(S[t] + dS, 1e-9);
      E[t+1] = fmax(E[t] + dE, 0.0);
      I[t+1] = fmax(I[t] + dI, 0.0);
      R[t+1] = fmax(R[t] + dR, 0.0);

      inc_new[t]    = fmax(-dS, 0.0);
      inc_SE[c, t]  = inc_new[t];        
      inc_EI[c,t+1] = delta * E[t+1];
    }

    for (t in 1:T) {
      real m = 0.0;
      int mlag = (t-1 < L) ? (t-1) : L;
      if (mlag > 0) {
        for (l in 1:mlag) m += inc_new[t - l] * p_delay[l];
      }
      mu[c,t] = IFR * fmax(m, 1e-12);
    }
  }
}



model {
  // hierarchical priors
  mu_log_R0    ~ normal(log(3.0), 0.5);
  sigma_log_R0 ~ normal(0, 0.5);

  mu_tm        ~ normal(40, 5);
  sigma_tm     ~ normal(0, 5);

  mu_log_k     ~ normal(log(0.35), 0.5);
  sigma_log_k  ~ normal(0, 0.5);

  mu_logit_cstar    ~ normal(logit(0.33), 0.7);
  sigma_logit_cstar ~ normal(0, 0.7);

  z_log_R0 ~ normal(0, 1);
  z_tm     ~ normal(0, 1);
  z_log_k  ~ normal(0, 1);
  z_logit_cstar ~ normal(0, 1);

  // initial conditions
  for (c in 1:C)
    log_I0[c] ~ normal(log_I0_mu[c], log_I0_sd[c]);

  // heterogeneity
  CV_shared ~ normal(1.6, 0.35);

  // dispersion
  phi ~ lognormal(log(15), 0.6);

  // likelihood
  for (c in 1:C)
    y[c] ~ neg_binomial_2(mu[c] + 1e-12, phi[c]);
}

generated quantities {
  int y_rep[C, T];
  matrix[C, T] log_lik;

  // store contact modifier c(t) for each country and time
  matrix[C, T] ct;

  real R0_pop_median = exp(mu_log_R0);
  real R0_pop_mean   = exp(mu_log_R0 + 0.5 * square(sigma_log_R0));

  real k_pop_median  = exp(mu_log_k);
  real k_pop_mean    = exp(mu_log_k + 0.5 * square(sigma_log_k));

  real tm_pop_median = mu_tm;
  real tm_pop_mean   = mu_tm;

  real cstar_pop_median = inv_logit(mu_logit_cstar);

  for (c in 1:C) {

    // ct[c,1] = 1; and ct[c,t+1] uses c_t computed at time t
    ct[c,1] = 1.0;
    for (t in 1:(T-1)) {
      real s_t = inv_logit(k[c] * (t - tm[c]));
      ct[c,t+1] = 1.0 - (1.0 - cstar[c]) * s_t;
    }

    for (t in 1:T) {
      y_rep[c,t]   = neg_binomial_2_rng(mu[c,t] + 1e-12, phi[c]);
      log_lik[c,t] = neg_binomial_2_lpmf(y[c,t] | mu[c,t] + 1e-12, phi[c]);
    }
  }
}

"
writeLines(stan_code_joint_log_het, "outputs/stan_fits/seir_nb_joint_logistic_uniform_HET_HIER.stan")

# ============================================================
# Joint hierarchical pooling (England+Scotland)
# Logistic NPI
# Heterogeneous susceptibility (CV estimated)
# LOO/WAIC + thesis-style diagnostics
# ============================================================

library(rstan)
library(loo)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Create output directories if they do not already exist
if (!dir.exists("outputs/figures")) {
  dir.create("outputs/figures", recursive = TRUE)
}
if (!dir.exists("outputs/stan_fits")) {
  dir.create("outputs/stan_fits", recursive = TRUE)
}

# -----------------------
# Data prep
# -----------------------
prepare_joint_data <- function(filepath = "chapter06_england_scotland/data/GB_data.csv") {
  d <- read.csv(filepath, stringsAsFactors = FALSE)
  d$Date <- dmy(d$Date)
  d <- d[order(d$Date), ]
  
  analysis_start <- as.Date("2020-01-31")
  analysis_end   <- as.Date("2020-06-01")
  idx <- which(d$Date >= analysis_start & d$Date <= analysis_end)
  
  dates <- d$Date[idx]
  deaths_EN <- d$Deaths_EN[idx]; deaths_EN[is.na(deaths_EN)] <- 0L
  deaths_SC <- d$Deaths_SC[idx]; deaths_SC[is.na(deaths_SC)] <- 0L
  
  list(
    england = list(deaths = as.integer(deaths_EN),
                   dates  = dates,
                   n_steps = length(dates),
                   population = 56e6),
    scotland = list(deaths = as.integer(deaths_SC),
                    dates  = dates,
                    n_steps = length(dates),
                    population = 5.5e6)
  )
}

joint_data <- prepare_joint_data("chapter06_england_scotland/data/GB_data.csv")

stopifnot(joint_data$england$n_steps == joint_data$scotland$n_steps)
T_len <- joint_data$england$n_steps
C <- 2L

y_mat <- rbind(
  as.integer(joint_data$england$deaths),
  as.integer(joint_data$scotland$deaths)
)
N_vec <- c(joint_data$england$population, joint_data$scotland$population)

# -----------------------
# Fixed process params
# -----------------------
rho   <- 0.5
delta <- 1/5.5
gamma <- 1/4
L     <- 60L

IFR <- 0.009

# Wood delay
lmu_id_fix <- 3.151
lsd_id_fix <- 0.469

# E0 priors (you can change)
log_E0_mu <- log(c(100, 10))
log_E0_sd <- c(0.1, 0.1)

log_I0_mu <- log_E0_mu - log(2.5)
log_I0_sd <- log_E0_sd


stan_data_common <- list(
  C = C,
  T = as.integer(T_len),
  y = y_mat,
  N = as.numeric(N_vec),
  rho = as.numeric(rho),
  delta = as.numeric(delta),
  gamma = as.numeric(gamma),
  L = as.integer(L),
  lmu_id_fix = as.numeric(lmu_id_fix),
  lsd_id_fix = as.numeric(lsd_id_fix),
  IFR = as.numeric(IFR),
  log_I0_mu  = as.numeric(log_I0_mu),
  log_I0_sd = as.numeric(log_I0_sd),
  CV_shared_fixed = 0.0
)

dates <- joint_data$england$dates
n_days <- length(dates)

# -----------------------
# Compile
# -----------------------
mod_het_log_hier <- stan_model("outputs/stan_fits/seir_nb_joint_logistic_uniform_HET_HIER.stan")

# -----------------------
# Sample
# -----------------------
fit_het_log_hier <- sampling(
  object  = mod_het_log_hier,
  data    = stan_data_common,
  seed    = 15012026,
  chains  = 4,
  iter    = 2500,
  warmup  = 1250,
  control = list(adapt_delta = 0.98, max_treedepth = 15)
)


rstan::check_hmc_diagnostics(fit_het_log_hier)

saveRDS(fit_het_log_hier, "outputs/stan_fits/fit_het_log_hier.rds")

# 1. Force the explicit conversion to a data frame
fit_df <- as.data.frame(fit_het_log_hier)

# 2. Save the draws to CSV
write.csv(fit_df, "fit_joint_hier_logistic_het_draws.csv", row.names = FALSE)
post_log_het <- rstan::extract(fit_het_log_hier
)
stopifnot(!is.null(post_log_het$inc_SE))

# draws x C x T
dim(post_log_het$inc_SE)

# save arrays for later overlay plotting
saveRDS(
  list(inc_SE = post_log_het$inc_SE, inc_EI = post_log_het$inc_EI),
  "logistic_hier_HET_arrays.rds",
  compress = "xz"
)


# -----------------------
# Choose fit to diagnose
# -----------------------
fit_log     <- fit_het_log_hier
model_label <- "Hierarchical Models"
include_CV  <- TRUE

# include_CV  <- FALSE

post <- rstan::extract(fit_log)

# -----------------------
# Quick parameter tables
# -----------------------
q <- function(x) quantile(x, c(0.025, 0.5, 0.975))
include_CV <- TRUE
post <- rstan::extract(fit_log)

# helper: coerce quantiles to plain numeric scalars
q025 <- function(x) as.numeric(stats::quantile(x, 0.025))
q050 <- function(x) as.numeric(stats::quantile(x, 0.50))
q975 <- function(x) as.numeric(stats::quantile(x, 0.975))

# ----------------------------
# Derived population midpoint: average country midpoint per draw
# ----------------------------
tm_bar <- rowMeans(post$tm)

# ----------------------------
# POPULATION table (MEDIAN everywhere)
# ----------------------------
tab_pop <- tibble::tibble(
  Parameter = c("R0_pop", "k_pop", "tm_bar (avg midpoint)", "cstar_pop",
                if (include_CV) "CV_shared" else NULL),
  Median = c(
    stats::median(post$R0_pop_median),
    stats::median(post$k_pop_median),
    stats::median(tm_bar),
    stats::median(post$cstar_pop_median),
    if (include_CV) stats::median(post$CV_shared) else NULL
  ),
  SD = c(
    stats::sd(post$R0_pop_median),
    stats::sd(post$k_pop_median),
    stats::sd(tm_bar),
    stats::sd(post$cstar_pop_median),
    if (include_CV) stats::sd(post$CV_shared) else NULL
  ),
  Q2.5 = c(
    q025(post$R0_pop_median),
    q025(post$k_pop_median),
    q025(tm_bar),
    q025(post$cstar_pop_median),
    if (include_CV) q025(post$CV_shared) else NULL
  ),
  Q50 = c(
    q050(post$R0_pop_median),
    q050(post$k_pop_median),
    q050(tm_bar),
    q050(post$cstar_pop_median),
    if (include_CV) q050(post$CV_shared) else NULL
  ),
  Q97.5 = c(
    q975(post$R0_pop_median),
    q975(post$k_pop_median),
    q975(tm_bar),
    q975(post$cstar_pop_median),
    if (include_CV) q975(post$CV_shared) else NULL
  )
)

# ----------------------------
# COUNTRY table (MEDIAN everywhere)
# ----------------------------
R0_EN <- post$R0[,1]; R0_SC <- post$R0[,2]
k_EN  <- post$k[,1];  k_SC  <- post$k[,2]
tm_EN <- post$tm[,1]; tm_SC <- post$tm[,2]
cs_EN <- post$cstar[,1]; cs_SC <- post$cstar[,2]
phi_EN <- post$phi[,1]; phi_SC <- post$phi[,2]
I0_EN<- exp(post$log_I0[,1]); I0_SC<-exp(post$log_I0[,2])
tab_cty <- tibble::tibble(
  Parameter = c("R0_EN","R0_SC","k_EN","k_SC","tm_EN","tm_SC",
                "cstar_EN","cstar_SC","phi_EN","phi_SC",
                if (include_CV) "CV_shared" else NULL),
  Median = c(
    stats::median(R0_EN), stats::median(R0_SC),
    stats::median(k_EN),  stats::median(k_SC),
    stats::median(tm_EN), stats::median(tm_SC),
    stats::median(cs_EN), stats::median(cs_SC),
    stats::median(phi_EN),stats::median(phi_SC),
    if (include_CV) stats::median(post$CV_shared) else NULL
  ),
  SD = c(
    stats::sd(R0_EN), stats::sd(R0_SC),
    stats::sd(k_EN),  stats::sd(k_SC),
    stats::sd(tm_EN), stats::sd(tm_SC),
    stats::sd(cs_EN), stats::sd(cs_SC),
    stats::sd(phi_EN),stats::sd(phi_SC),
    if (include_CV) stats::sd(post$CV_shared) else NULL
  ),
  Q2.5 = c(
    q025(R0_EN), q025(R0_SC),
    q025(k_EN),  q025(k_SC),
    q025(tm_EN), q025(tm_SC),
    q025(cs_EN), q025(cs_SC),
    q025(phi_EN),q025(phi_SC),
    if (include_CV) q025(post$CV_shared) else NULL
  ),
  Q50 = c(
    q050(R0_EN), q050(R0_SC),
    q050(k_EN),  q050(k_SC),
    q050(tm_EN), q050(tm_SC),
    q050(cs_EN), q050(cs_SC),
    q050(phi_EN),q050(phi_SC),
    if (include_CV) q050(post$CV_shared) else NULL
  ),
  Q97.5 = c(
    q975(R0_EN), q975(R0_SC),
    q975(k_EN),  q975(k_SC),
    q975(tm_EN), q975(tm_SC),
    q975(cs_EN), q975(cs_SC),
    q975(phi_EN),q975(phi_SC),
    if (include_CV) q975(post$CV_shared) else NULL
  )
)

cat("\n== POPULATION (medians + 95% CrI) ==\n")
print(tab_pop %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4))))

cat("\n== COUNTRY-SPECIFIC (medians + 95% CrI) ==\n")
print(tab_cty %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4))))







cat("\n== POPULATION (hierarchical centers) ==\n")
print(tab_pop %>% mutate(across(where(is.numeric), ~round(., 4))))

cat("\n== COUNTRY-SPECIFIC ==\n")
print(tab_cty %>% mutate(across(where(is.numeric), ~round(., 4))))






# -----------------------
# Fitted trajectories + NB predictive bands
# -----------------------
mu_EN_med <- apply(post$mu[,1,], 2, median)
mu_SC_med <- apply(post$mu[,2,], 2, median)

phi_EN_med <- median(phi_EN)
phi_SC_med <- median(phi_SC)

band_en <- list(
  mu = mu_EN_med,
  lo = qnbinom(0.025, size = phi_EN_med, mu = mu_EN_med),
  hi = qnbinom(0.975, size = phi_EN_med, mu = mu_EN_med)
)
band_sc <- list(
  mu = mu_SC_med,
  lo = qnbinom(0.025, size = phi_SC_med, mu = mu_SC_med),
  hi = qnbinom(0.975, size = phi_SC_med, mu = mu_SC_med)
)

fitted_data <- data.frame(
  Date = rep(dates, 2),
  Country = rep(c("England","Scotland"), each = n_days),
  Observed = c(joint_data$england$deaths, joint_data$scotland$deaths),
  Fitted = c(mu_EN_med, mu_SC_med),
  Lower  = c(band_en$lo, band_sc$lo),
  Upper  = c(band_en$hi, band_sc$hi)
)

country_colors <- c("England"="#1f77b4","Scotland"="#ff7f0e")


#country_colors <- c("England"="grey","Scotland"="skyblue")


p_fit <- ggplot(fitted_data, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Country), alpha = 0.3) +
  geom_line(aes(y = Fitted, color = Country), linewidth = 1.1) +
  geom_point(aes(y = Observed, color = Country), size = 1, alpha = 0.7) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors) +
  labs(title = paste0("Fitted vs observed : ", model_label),
       #subtitle = "Line = posterior median mu; ribbon = NB( mu, phi_med ) 95% interval",
       x = "Date", y = "Daily deaths") +
  create_custom_theme() +
  theme(legend.position = "bottom") +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_fit)

# -----------------------
# Residual diagnostics
# -----------------------
res_EN <- joint_data$england$deaths - mu_EN_med
res_SC <- joint_data$scotland$deaths - mu_SC_med

rmse_EN <- sqrt(mean(res_EN^2))
rmse_SC <- sqrt(mean(res_SC^2))

resid_data <- data.frame(
  Date = rep(dates, 2),
  Country = rep(c("England","Scotland"), each = n_days),
  Residual = c(res_EN, res_SC),
  Fitted = c(mu_EN_med, mu_SC_med)
)

p_res_time <- ggplot(resid_data, aes(x = Date, y = Residual, color = Country)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = country_colors) +
  labs(title = paste0("Residuals vs time : ", model_label),
       subtitle = paste0("RMSE: EN=", round(rmse_EN,2), ", SC=", round(rmse_SC,2)),
       x = "Date", y = "Observed - fitted") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y")

p_res_fit <- ggplot(resid_data, aes(x = Fitted, y = Residual, color = Country)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = country_colors) +
  labs(title = paste0("Residuals vs fitted: ", model_label),
       x = "Fitted", y = "Observed - fitted") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free")

print(p_res_time)
print(p_res_fit)

# -----------------------
# Trace plots (key parameters)
# -----------------------
params_to_trace <- c("k[1]","k[2]","tm[1]","tm[2]")#,"cstar[1]","cstar[2]","phi[1]","phi[2]")#,"I0[1]"  ,"I0[2]"   )
if (include_CV) params_to_trace <- c(params_to_trace, "CV_shared")

arr <- as.array(fit_log, pars = params_to_trace)
it <- dim(arr)[1]; ch <- dim(arr)[2]; pn <- dim(arr)[3]
par_names <- dimnames(arr)[[3]]

trace_df <- data.frame(
  Iteration = rep(1:it, times = ch*pn),
  Chain = factor(rep(rep(1:ch, each = it), times = pn)),
  Parameter = factor(rep(par_names, each = it*ch), levels = par_names),
  Value = as.vector(arr)
)

p_trace <- ggplot(trace_df, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.75) +
  labs(title = paste0("Trace plots - ", model_label),
       x = "Iteration", y = "Value") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)




# --- Trace plots (key parameters) ---
params_to_trace <- c("R0[1]","R0[2]","k[1]","k[2]","tm[1]","tm[2]","cstar[1]","cstar[2]","phi[1]","phi[2]")
if (include_CV) params_to_trace <- c(params_to_trace, "CV_shared")

arr <- as.array(fit_log, pars = params_to_trace)
it <- dim(arr)[1]; ch <- dim(arr)[2]; pn <- dim(arr)[3]
par_names <- dimnames(arr)[[3]]

trace_df <- data.frame(
  Iteration = rep(1:it, times = ch*pn),
  Chain = factor(rep(rep(1:ch, each = it), times = pn)),
  Parameter = factor(rep(par_names, each = it*ch), levels = par_names),
  Value = as.vector(arr)
) %>%
  mutate(Parameter = ifelse(Parameter == "CV_shared", "CV", as.character(Parameter))) %>%
  mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

p_trace <- ggplot(trace_df, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.75) +
  labs(title = paste0("Trace plots — ", model_label),
       x = "Iteration", y = "Value") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)



# -----------------------
# Posterior histograms  JOINT
# -----------------------
draw_list <- list(
  "R0[1]"=R0_EN, "R0[2]"=R0_SC,
  "k[1]"=k_EN,   "k[2]"=k_SC,
  "tm[1]"=tm_EN, "tm[2]"=tm_SC,
  "cstar[1]"=cs_EN, "cstar[2]"=cs_SC,
  "phi[1]"=phi_EN, "phi[2]"=phi_SC,
  "I0[1]"=I0_EN, "I0[2]"=I0_SC
)
if (include_CV) draw_list[["CV"]] <- as.numeric(post$CV_shared)   # <- rename here

post_dist <- bind_rows(lapply(names(draw_list), function(nm) {
  tibble(Parameter = nm, Value = as.numeric(draw_list[[nm]]))
}))

meds <- post_dist %>% group_by(Parameter) %>% summarise(med = median(Value), .groups="drop")

p_hist <- ggplot(post_dist, aes(x = Value)) +
  geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue", color = "black") +
  geom_vline(data = meds, aes(xintercept = med), linetype = "dashed", linewidth = 1) +
  labs(title = paste0("Posterior distributions — ", model_label),
       subtitle = "Dashed = posterior median",
       x = "Value", y = "Count") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_hist)

# -----------------------
# Posterior histograms 
# -----------------------
draw_list <- list(
  "R0[1]"=R0_EN, "R0[2]"=R0_SC,
  "k[1]"=k_EN,   "k[2]"=k_SC,
  "tm[1]"=tm_EN, "tm[2]"=tm_SC,
  "cstar[1]"=cs_EN, "cstar[2]"=cs_SC,
  "phi[1]"=phi_EN, "phi[2]"=phi_SC,
  "I0[1]"=I0_EN, "I0[2]"=I0_SC
)
if (include_CV) draw_list[["CV_shared"]] <- post$CV_shared

post_dist <- bind_rows(lapply(names(draw_list), function(nm) {
  tibble(Parameter = nm, Value = as.numeric(draw_list[[nm]]))
}))

p_hist <- ggplot(post_dist, aes(x = Value)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "steelblue") +
  geom_vline(data = post_dist %>% group_by(Parameter) %>% summarise(med = median(Value), .groups="drop"),
             aes(xintercept = med), linetype = "dashed", linewidth = 1) +
  labs(title = paste0("Posterior distributions — ", model_label),
       subtitle = "Dashed = posterior median",
       x = "Value", y = "Count") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_hist)


# -----------------------
# Incidence (EI flow)
# -----------------------
inc_EN <- post$inc_EI[,1,]
inc_SC <- post$inc_EI[,2,]

inc_sum <- bind_rows(
  as.data.frame(inc_EN) %>%
    pivot_longer(cols = everything()) %>%
    mutate(Country="England", time = as.integer(gsub("V","", name))) %>%
    transmute(Country, time, inc = value),
  as.data.frame(inc_SC) %>%
    pivot_longer(cols = everything()) %>%
    mutate(Country="Scotland", time = as.integer(gsub("V","", name))) %>%
    transmute(Country, time, inc = value)
) %>%
  group_by(Country, time) %>%
  summarise(Med = median(inc),
            Lo  = quantile(inc, 0.025),
            Hi  = quantile(inc, 0.975),
            .groups="drop") %>%
  mutate(Date = min(dates) + (time - 1))

lines_df <- tibble(
  Country = c("England", "Scotland"),   # <- must match facet variable name
  line    = as.Date(c("2020-03-26", "2020-03-24"))
)

p_inc <- ggplot(inc_sum, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.3) +
  geom_line(aes(y = Med), linewidth = 1.1) +
  geom_vline(data = lines_df, aes(xintercept = line), linewidth = 0.5) +
  labs(title = paste0("Simulated incidence — ", model_label),
       x = "Date", y = "delta * E(t)") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_inc)
