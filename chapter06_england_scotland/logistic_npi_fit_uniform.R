
# =============================================================================
# logistic_npi_fit_uniform.R
#
# Chapter 6: England and Scotland application
# Joint fit with a LOGISTIC NPI profile, UNIFORM (flat) priors,
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu).
#
# Sensitivity companion to logistic_npi_fit.R (hierarchical priors).
#
# Outputs:
#   outputs/stan_fits/fit_het_log_unif.rds
#   outputs/stan_fits/loo_het_log_unif.rds
#   outputs/figures/
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

library(rstan)

library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(bayesplot)
library(gridExtra)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (!dir.exists("outputs/stan_fits")) dir.create("outputs/stan_fits",recursive=TRUE)


if (!dir.exists("outputs/figures")) dir.create("outputs/figures", recursive=TRUE)

# ============================================================
# 1) Data prep 
# ============================================================
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
                   n_steps= length(dates),
                   population = 56e6),
    scotland= list(deaths = as.integer(deaths_SC),
                   dates  = dates,
                   n_steps= length(dates),
                   population = 5.5e6)
  )
}

joint_data <- prepare_joint_data("chapter06_england_scotland/data/GB_data.csv")
stopifnot(joint_data$england$n_steps == joint_data$scotland$n_steps)

dates <- joint_data$england$dates
T_len <- joint_data$england$n_steps
C <- 2L

y_mat <- rbind(
  as.integer(joint_data$england$deaths),
  as.integer(joint_data$scotland$deaths)
)
N_vec <- c(joint_data$england$population, joint_data$scotland$population)

# ============================================================
#  Fixed process + observation inputs
# ============================================================
rho   <- 0.5
delta <- 1/5.5
gamma <- 1/4

L     <- 60L
lmu_id_fix <- 3.151
lsd_id_fix <- 0.469

IFR_fixed  <- 0.009

# ============================================================
#  Uniform bounds (priors via parameter bounds)
# ============================================================
bounds <- list(
  R0_min = 2,     R0_max = 8.0,
  k_min  = 0.01,  k_max  = 2.0,
  tm_lower = 1.0, tm_upper = as.numeric(T_len),
  cstar_min = 0.05, cstar_max = 0.9,
  CV_min = 0.01,  CV_max = 3.0,          # only used in HET
  phi_min = 1.0, phi_max = 1000.0,
  I0_min = 1,    I0_max = 10000.0,
  E0_multiplier_min = 1,
  E0_multiplier_max = 4,
  E0_multiplier = 2.5
)

# ============================================================
# Stan code 
# ============================================================


stan_code_joint_log_uniform_het <- "
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

  real<lower=0, upper=1> IFR_fixed;

  // Uniform bounds
  real<lower=0> R0_min;  real<lower=R0_min> R0_max;
  real<lower=0> k_min;   real<lower=k_min>  k_max;
  real tm_lower;         real tm_upper;
  real<lower=0, upper=1> cstar_min;  real<lower=cstar_min, upper=1> cstar_max;

  real<lower=0> CV_min;  real<lower=CV_min> CV_max;
  real<lower=0> phi_min; real<lower=phi_min> phi_max;

  real<lower=0> I0_min;  real<lower=I0_min> I0_max;
  real<lower=1> E0_multiplier_min;
  real<upper=10> E0_multiplier_max;
  real<lower=E0_multiplier_min, upper=E0_multiplier_max> E0_multiplier;
}

parameters {
  vector<lower=R0_min,  upper=R0_max>[C] R0;
  vector<lower=k_min,   upper=k_max>[C]  k;
  vector<lower=tm_lower, upper=tm_upper>[C] tm;
  vector<lower=cstar_min, upper=cstar_max>[C] cstar;

  real<lower=CV_min, upper=CV_max> CV_shared;

  vector<lower=phi_min, upper=phi_max>[C] phi;

  vector<lower=I0_min, upper=I0_max>[C] I0;
}

transformed parameters {
  real CV2 = square(CV_shared);
  vector[L] p_delay = delay_probs(L, lmu_id_fix, lsd_id_fix);

  matrix[C, T] mu;
  matrix<lower=0>[C, T] inc_EI;
  inc_EI = rep_matrix(0.0, C, T);

  vector[C] E0 = E0_multiplier * I0;

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
      real s_t = inv_logit(k[c] * (t - tm[c]));
      real c_t = 1.0 - (1.0 - cstar[c]) * s_t;

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

    for (t in 1:T) {
      real m = 0.0;
      int mlag = (t-1 < L) ? t-1 : L;
      if (mlag > 0)
        for (l in 1:mlag) m += inc_new[t - l] * p_delay[l];
      mu[c, t] = IFR_fixed * fmax(m, 1e-12);
    }
  }
}

model {
  // Uniform priors via bounds
  for (c in 1:C)
    y[c] ~ neg_binomial_2(mu[c] + 1e-12, phi[c]);
}

generated quantities {
  int y_rep[C, T];
  matrix[C, T] log_lik;
  matrix[C, T] ct;

  for (c in 1:C) {
    ct[c,1] = 1.0;
    for (t in 1:(T-1)) {
      real s_t = inv_logit(k[c] * (t - tm[c]));
      ct[c,t+1] = 1.0 - (1.0 - cstar[c]) * s_t;
    }

    for (t in 1:T) {
      y_rep[c, t] = neg_binomial_2_rng(mu[c, t] + 1e-12, phi[c]);
      log_lik[c, t] = neg_binomial_2_lpmf(y[c, t] | mu[c, t] + 1e-12, phi[c]);
    }
  }
}
"


writeLines(stan_code_joint_log_uniform_het, "outputs/stan_fits/seir_nb_joint_linear_uniform_het_CVcommon_IFRfixed_I0est.stan")




# ============================================================
# Stan data
# ============================================================
stan_data_common <- list(
  C = as.integer(C),
  T = as.integer(T_len),
  y = y_mat,
  N = as.numeric(N_vec),
  
  rho = as.numeric(rho),
  delta = as.numeric(delta),
  gamma = as.numeric(gamma),
  
  L = as.integer(L),
  lmu_id_fix = as.numeric(lmu_id_fix),
  lsd_id_fix = as.numeric(lsd_id_fix),
  
  IFR_fixed = as.numeric(IFR_fixed),
  
  R0_min = as.numeric(bounds$R0_min), R0_max = as.numeric(bounds$R0_max),
  k_min  = as.numeric(bounds$k_min),  k_max  = as.numeric(bounds$k_max),
  tm_lower = as.numeric(bounds$tm_lower), tm_upper = as.numeric(bounds$tm_upper),
  cstar_min = as.numeric(bounds$cstar_min), cstar_max = as.numeric(bounds$cstar_max),
  
  CV_min = as.numeric(bounds$CV_min), CV_max = as.numeric(bounds$CV_max),
  phi_min = as.numeric(bounds$phi_min), phi_max = as.numeric(bounds$phi_max),
  
  I0_min = as.numeric(bounds$I0_min), I0_max = as.numeric(bounds$I0_max),
  E0_multiplier_min = as.numeric(bounds$E0_multiplier_min),
  E0_multiplier_max = as.numeric(bounds$E0_multiplier_max),
  E0_multiplier = as.numeric(bounds$E0_multiplier)
)


# ============================================================
# Compile + sample
# ============================================================
mod_het_log_unif <- stan_model("outputs/stan_fits/seir_nb_joint_logistic_uniform_HET_BASE.stan")
#mod_het_log_unif <- stan_model("seir_nb_joint_logistic_uniform_HET_BASE.stan")

fit_het_log_unif <- sampling(
  object  = mod_het_log_unif,
  data    = stan_data_common,
  seed    = 15012026,
  chains  = 4,
  iter    = 2500,
  warmup  = 1250,
  control = list(adapt_delta = 0.98, max_treedepth = 15)
)


rstan::check_hmc_diagnostics(fit_het_log_unif)

saveRDS(fit_het_log_unif, "outputs/stan_fits/fit_het_log_unif.rds")





# ============================================================
#  Thesis plotting theme 
# ============================================================
create_custom_theme <- function(
    base_size = 12,
    title_size = 14,
    subtitle_size = 11,
    axis_title_size = 12,
    axis_text_size = 10,
    legend_title_size = 11,
    legend_text_size = 10,
    legend_position = "bottom",
    grid_color_major = "grey85",
    grid_color_minor = "grey92",
    remove_minor_grid = TRUE,
    base_family = ""
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
      
      axis.line        = element_line(colour = "grey60", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      
      strip.text       = element_text(size = axis_text_size, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

country_colors <- c("England"="#1f77b4","Scotland"="#ff7f0e")

# ============================================================
#  Diagnostics
# ============================================================
fit_use     <- fit_het_log_unif
model_label <- "Joint logistic:Uniform priors"
include_CV  <- TRUE

post <- rstan::extract(fit_use)
q <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

# ============================================================
# 9) Posterior summaries (quick tables)
# ============================================================
tab_country <- tibble(
  Parameter = c("R0[1]","R0[2]","k[1]","k[2]","tm[1]","tm[2]","cstar[1]","cstar[2]","phi[1]","phi[2]",
                if (include_CV) "CV_shared" else NULL),
  Mean  = c(colMeans(post$R0), colMeans(post$k), colMeans(post$tm), colMeans(post$cstar), colMeans(post$phi),
            if (include_CV) mean(post$CV_shared) else NULL),
  SD    = c(apply(post$R0,2,sd), apply(post$k,2,sd), apply(post$tm,2,sd), apply(post$cstar,2,sd), apply(post$phi,2,sd),
            if (include_CV) sd(post$CV_shared) else NULL),
  Q2.5  = c(apply(post$R0,2,function(x) q(x)[1]),
            apply(post$k,2,function(x) q(x)[1]),
            apply(post$tm,2,function(x) q(x)[1]),
            apply(post$cstar,2,function(x) q(x)[1]),
            apply(post$phi,2,function(x) q(x)[1]),
            if (include_CV) q(post$CV_shared)[1] else NULL),
  Q50   = c(apply(post$R0,2,function(x) q(x)[2]),
            apply(post$k,2,function(x) q(x)[2]),
            apply(post$tm,2,function(x) q(x)[2]),
            apply(post$cstar,2,function(x) q(x)[2]),
            apply(post$phi,2,function(x) q(x)[2]),
            if (include_CV) q(post$CV_shared)[2] else NULL),
  Q97.5 = c(apply(post$R0,2,function(x) q(x)[3]),
            apply(post$k,2,function(x) q(x)[3]),
            apply(post$tm,2,function(x) q(x)[3]),
            apply(post$cstar,2,function(x) q(x)[3]),
            apply(post$phi,2,function(x) q(x)[3]),
            if (include_CV) q(post$CV_shared)[3] else NULL)
)


print(tab_country %>% mutate(across(where(is.numeric), ~round(. , 4))), n = Inf)
#write.csv(tab_country, "posterior_summaries_joint_logistic_uniform.csv", row.names = FALSE)


post <- rstan::extract(fit_use)

qfun <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

# ============================================================
# Posterior summaries: country-level parameters + I0
# ============================================================

tab_country <- tibble(
  Parameter = c(
    "R0[1]", "R0[2]",
    "k[1]", "k[2]",
    "tm[1]", "tm[2]",
    "cstar[1]", "cstar[2]",
    "I0[1]", "I0[2]",
    "phi[1]", "phi[2]",
    if (include_CV) "CV_shared" else NULL
  ),
  
  Mean = c(
    colMeans(post$R0),
    colMeans(post$k),
    colMeans(post$tm),
    colMeans(post$cstar),
    colMeans(post$I0),
    colMeans(post$phi),
    if (include_CV) mean(post$CV_shared) else NULL
  ),
  
  SD = c(
    apply(post$R0, 2, sd),
    apply(post$k, 2, sd),
    apply(post$tm, 2, sd),
    apply(post$cstar, 2, sd),
    apply(post$I0, 2, sd),
    apply(post$phi, 2, sd),
    if (include_CV) sd(post$CV_shared) else NULL
  ),
  
  Q2.5 = c(
    apply(post$R0, 2, function(x) qfun(x)[1]),
    apply(post$k, 2, function(x) qfun(x)[1]),
    apply(post$tm, 2, function(x) qfun(x)[1]),
    apply(post$cstar, 2, function(x) qfun(x)[1]),
    apply(post$I0, 2, function(x) qfun(x)[1]),
    apply(post$phi, 2, function(x) qfun(x)[1]),
    if (include_CV) qfun(post$CV_shared)[1] else NULL
  ),
  
  Q50 = c(
    apply(post$R0, 2, function(x) qfun(x)[2]),
    apply(post$k, 2, function(x) qfun(x)[2]),
    apply(post$tm, 2, function(x) qfun(x)[2]),
    apply(post$cstar, 2, function(x) qfun(x)[2]),
    apply(post$I0, 2, function(x) qfun(x)[2]),
    apply(post$phi, 2, function(x) qfun(x)[2]),
    if (include_CV) qfun(post$CV_shared)[2] else NULL
  ),
  
  Q97.5 = c(
    apply(post$R0, 2, function(x) qfun(x)[3]),
    apply(post$k, 2, function(x) qfun(x)[3]),
    apply(post$tm, 2, function(x) qfun(x)[3]),
    apply(post$cstar, 2, function(x) qfun(x)[3]),
    apply(post$I0, 2, function(x) qfun(x)[3]),
    apply(post$phi, 2, function(x) qfun(x)[3]),
    if (include_CV) qfun(post$CV_shared)[3] else NULL
  )
)


print(tab_country %>% mutate(across(where(is.numeric), ~ round(.x, 4))), n = Inf)


post <- rstan::extract(fit_use)

summ_vec <- function(x, labels) {
  tibble(
    Parameter = labels,
    Mean  = colMeans(x),
    SD    = apply(x, 2, sd),
    Q2.5  = apply(x, 2, quantile, probs = 0.025, na.rm = TRUE),
    Q50   = apply(x, 2, quantile, probs = 0.500, na.rm = TRUE),
    Q97.5 = apply(x, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

summ_scalar <- function(x, label) {
  tibble(
    Parameter = label,
    Mean  = mean(x, na.rm = TRUE),
    SD    = sd(x, na.rm = TRUE),
    Q2.5  = quantile(x, 0.025, na.rm = TRUE),
    Q50   = quantile(x, 0.500, na.rm = TRUE),
    Q97.5 = quantile(x, 0.975, na.rm = TRUE)
  )
}

tab_country <- bind_rows(
  summ_vec(post$R0,    c("R0[1]", "R0[2]")),
  summ_vec(post$k,     c("k[1]", "k[2]")),
  summ_vec(post$tm,    c("tm[1]", "tm[2]")),
  summ_vec(post$cstar, c("cstar[1]", "cstar[2]")),
  summ_vec(post$I0,    c("I0[1]", "I0[2]")),
  summ_vec(post$phi,   c("phi[1]", "phi[2]")),
  if (include_CV) summ_scalar(post$CV_shared, "CV_shared") else NULL
)


print(tab_country %>% mutate(across(where(is.numeric), ~ round(.x, 4))), n = Inf)

# ============================================================
#  Fitted vs observed + NB predictive bands
# ============================================================
n_days <- length(dates)

mu_EN_med <- apply(post$mu[,1,], 2, median)
mu_SC_med <- apply(post$mu[,2,], 2, median)

phi_EN_med <- median(post$phi[,1])
phi_SC_med <- median(post$phi[,2])

band_en <- list(
  lo = qnbinom(0.025, size = phi_EN_med, mu = mu_EN_med),
  hi = qnbinom(0.975, size = phi_EN_med, mu = mu_EN_med)
)
band_sc <- list(
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

p_fit <- ggplot(fitted_data, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Country), alpha = 0.3) +
  geom_line(aes(y = Fitted, color = Country), linewidth = 1.1) +
  geom_point(aes(y = Observed, color = Country), size = 1, alpha = 0.7) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors) +
  labs(title = paste0("Fitted vs observed: ", model_label),
       x = "Date", y = "Daily deaths") +
  create_custom_theme() +
  theme(legend.position = "bottom") +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_fit)
#ggsave("outputs/figures/fitted_vs_observed_joint_logistic_uniform.png", p_fit, width = 10, height = 7, dpi = 300)

# ============================================================
# 11) Residual diagnostics
# ============================================================
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
  labs(title = paste0("Residuals vs time — ", model_label),
       subtitle = paste0("RMSE: EN=", round(rmse_EN,2), ", SC=", round(rmse_SC,2)),
       x = "Date", y = "Observed - fitted") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y")


print(p_res_time)

# ============================================================
#  Trace plots 
# ============================================================
params_to_trace <- c(
  "R0[1]","R0[2]",
  "k[1]","k[2]",
  "tm[1]","tm[2]",
  "cstar[1]","cstar[2]",
  "phi[1]","phi[2]"
)
if (include_CV) params_to_trace <- c(params_to_trace, "CV_shared")

arr <- as.array(fit_use, pars = params_to_trace)
it <- dim(arr)[1]; ch <- dim(arr)[2]; pn <- dim(arr)[3]
par_names <- dimnames(arr)[[3]]

trace_df <- data.frame(
  Iteration = rep(1:it, times = ch*pn),
  Chain = factor(rep(rep(1:ch, each = it), times = pn)),
  Parameter = factor(rep(par_names, each = it*ch), levels = par_names),
  Value = as.vector(arr)
)

p_trace <- ggplot(trace_df, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.75, linewidth = 0.35) +
  labs(title = paste0("Trace plots — ", model_label),
       x = "Iteration", y = "Value") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)

# ============================================================
#  Posterior densities 
# ============================================================
draw_list <- list(
  "R0[1]"    = as.numeric(post$R0[,1]),
  "R0[2]"    = as.numeric(post$R0[,2]),
  "k[1]"     = as.numeric(post$k[,1]),
  "k[2]"     = as.numeric(post$k[,2]),
  "tm[1]"    = as.numeric(post$tm[,1]),
  "tm[2]"    = as.numeric(post$tm[,2]),
  "cstar[1]" = as.numeric(post$cstar[,1]),
  "cstar[2]" = as.numeric(post$cstar[,2]),
  "phi[1]"   = as.numeric(post$phi[,1]),
  "phi[2]"   = as.numeric(post$phi[,2])
)
if (include_CV) draw_list[["CV"]] <- as.numeric(post$CV_shared)

post_dist <- bind_rows(lapply(names(draw_list), function(nm) {
  tibble(Parameter = nm, Value = as.numeric(draw_list[[nm]]))
})) %>%
  mutate(Parameter = factor(Parameter, levels = names(draw_list)))

meds <- post_dist %>%
  group_by(Parameter) %>%
  summarise(med = median(Value), .groups="drop")

p_dens <- ggplot(post_dist, aes(x = Value)) +
  geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.6) +
  geom_vline(data = meds, aes(xintercept = med), linetype = "dashed", linewidth = 1) +
  labs(title = paste0("Posterior densities — ", model_label),
       subtitle = "Dashed = posterior median",
       x = "Value", y = "Density") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_dens)
ggsave("outputs/figures/density_joint_logistic_uniform.png", p_dens, width = 12, height = 8, dpi = 300)
