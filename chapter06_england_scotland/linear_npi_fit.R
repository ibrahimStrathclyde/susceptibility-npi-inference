# =============================================================================
# linear_npi_fit.R
#
# Chapter 6: England and Scotland application
# Joint hierarchical fit with a PIECEWISE-LINEAR NPI ramp,
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu).
#
# Fits England and Scotland jointly. R0 and c* are partially pooled across the
# two countries; nu is shared. The contact modifier c(t) declines linearly from
# 1 at t0 to c* at t1 (the lockdown date), then holds at c*.
#
# Data: chapter06_england_scotland/data/GB_data.csv
#       31 Jan 2020 to 01 Jun 2020 (first wave)
#
# Outputs:
#   outputs/stan_fits/seir_nb_joint_linear_hier_HET.stan   Stan model
#   outputs/stan_fits/fit_het_lin_hier.rds                 fitted object
#   outputs/stan_fits/loo_het_lin_hier.rds                 LOO object
#   outputs/stan_fits/linear_hier_HET_arrays.rds           incidence draws
#   outputs/figures/                                        diagnostic plots
#
# The fitted object is named fit_het_lin_hier. This name is required by
# compare_npi_specifications.R. Do not rename it.
#
# Sampler: 4 chains, 2500 iterations, 1250 warmup,
#          adapt_delta = 0.98, max_treedepth = 15.
# Run time: roughly 30 to 45 minutes  depending on hardware.
#
# Authors: Ibrahim Mohammed, Chris Robertson
# =============================================================================

library(rstan)
library(loo)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)

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
# Theme
# -----------------------
create_custom_theme <- function(
    title_size = 16,
    subtitle_size = 13,
    axis_title_size = 13,
    axis_text_size = 12,
    legend_title_size = 12,
    legend_text_size = 11,
    legend_position = "right",
    base_size = 12,
    base_family = "",
    grid_color_major = "grey90",
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
      panel.grid.minor = if (remove_minor_grid) element_blank() else element_line(colour = "grey95"),
      axis.line        = element_line(colour = "grey60", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(size = axis_text_size, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

# =========================================================
# 0) Stan code — HIERARCHICAL LINEAR RAMP (HET)
# =========================================================
stan_code_joint_lin_hier_het <- "
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

  // lockdown start index per country
  int<lower=1, upper=T> t1[C];

  // country-specific prior for log(I0)
  vector[C] log_I0_mu;
  vector<lower=0>[C] log_I0_sd;
}

parameters {
  // hierarchical hyperparameters
  real  mu_log_R0;
  real<lower=0> sigma_log_R0;

  real  mu_log_w;
  real<lower=0> sigma_log_w;

  real  mu_logit_cstar;
  real<lower=0> sigma_logit_cstar;

  // non-centred country effects
  vector[C] z_log_R0;
  vector[C] z_log_w;
  vector[C] z_logit_cstar;

  // initial condition (country-specific)
  vector[C] log_I0;

  // heterogeneity (shared)
  real<lower=0> CV_shared;

  // NB dispersion (country-specific)
  vector<lower=1e-3>[C] phi;
}
transformed parameters {
  vector<lower=0>[C] R0    = exp(mu_log_R0 + sigma_log_R0 * z_log_R0);
  vector<lower=0>[C] w     = exp(mu_log_w  + sigma_log_w  * z_log_w);
  vector<lower=0,upper=1>[C] cstar =
    inv_logit(mu_logit_cstar + sigma_logit_cstar * z_logit_cstar);

  vector<lower=0>[C] I0 = exp(log_I0);
  vector<lower=0>[C] E0 = 2.5 * I0;

  vector[C] t0;
  for (c in 1:C) t0[c] = t1[c] - w[c];

  real CV2 = square(CV_shared);

  vector[L] p_delay = delay_probs(L, lmu_id_fix, lsd_id_fix);

  matrix[C, T] mu;

  // EI flow (delta * E) already present
  matrix<lower=0>[C, T] inc_EI;
  inc_EI = rep_matrix(0.0, C, T);

  // NEW: infections flow S->E (same as inc_new)
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
    // inc_SE[c, 1] stays 0.0 by initialization

    for (t in 1:(T-1)) {
      real c_t;

      if (t < t0[c]) {
        c_t = 1.0;
      } else if (t < t1[c]) {
        real denom = fmax(t1[c] - t0[c], 1e-6);  // ~w[c], safer
        real frac  = (t - t0[c]) / denom;
        c_t = 1.0 - (1.0 - cstar[c]) * fmin(fmax(frac, 0.0), 1.0);
      } else {
        c_t = cstar[c];
      }

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

      inc_new[t]     = fmax(-dS, 0.0);
      inc_SE[c, t]   = inc_new[t];     // NEW: store infections S->E
      inc_EI[c, t+1] = delta * E[t+1];
    }

    // infections -> deaths (convolution)
    for (t in 1:T) {
      real m = 0.0;
      int mlag = (t-1 < L) ? (t-1) : L;
      if (mlag > 0) {
        for (l in 1:mlag) m += inc_new[t - l] * p_delay[l];
      }
      mu[c, t] = IFR * fmax(m, 1e-12);
    }
  }
}


model {
  // priors (same style as your stage-2 logistic hier file)
  mu_log_R0    ~ normal(log(3.0), 0.5);
  sigma_log_R0 ~ normal(0, 0.5);

  // NOTE: prior for w is a modelling choice — adjust if you want it tighter/looser.
  // This default is weakly-informative on a positive scale.
  mu_log_w     ~ normal(log(20.0), 1.0);
  sigma_log_w  ~ normal(0, 1.0);

  mu_logit_cstar    ~ normal(logit(0.33), 0.7);
  sigma_logit_cstar ~ normal(0, 0.7);

  z_log_R0       ~ normal(0, 1);
  z_log_w        ~ normal(0, 1);
  z_logit_cstar  ~ normal(0, 1);

  for (c in 1:C)
    log_I0[c] ~ normal(log_I0_mu[c], log_I0_sd[c]);

  CV_shared ~ normal(1.6, 0.35);

  phi ~ lognormal(log(15), 0.6);

  for (c in 1:C)
    y[c] ~ neg_binomial_2(mu[c] + 1e-12, phi[c]);
}

generated quantities {
  int y_rep[C, T];
  matrix[C, T] log_lik;

  // store contact modifier c_t for each country and time
  matrix[C, T] ct;

  // population summaries
  real R0_pop_median = exp(mu_log_R0);
  real R0_pop_mean   = exp(mu_log_R0 + 0.5 * square(sigma_log_R0));

  real w_pop_median  = exp(mu_log_w);
  real w_pop_mean    = exp(mu_log_w + 0.5 * square(sigma_log_w));

  real cstar_pop_median = inv_logit(mu_logit_cstar);

  for (c in 1:C) {
    ct[c,1] = 1.0;
    for (t in 1:(T-1)) {
      real t0_c = t1[c] - exp(mu_log_w + sigma_log_w * z_log_w[c]);
      real w_c  = exp(mu_log_w + sigma_log_w * z_log_w[c]);
      real c_t;

      if (t < t0_c) {
        c_t = 1.0;
      } else if (t < t1[c]) {
        real denom = fmax(t1[c] - t0_c, 1e-6);
        real frac  = (t - t0_c) / denom;
        c_t = 1.0 - (1.0 - inv_logit(mu_logit_cstar + sigma_logit_cstar * z_logit_cstar[c]))
                    * fmin(fmax(frac, 0.0), 1.0);
      } else {
        c_t = inv_logit(mu_logit_cstar + sigma_logit_cstar * z_logit_cstar[c]);
      }

      ct[c, t+1] = c_t;
    }
  }

  for (c in 1:C) {
    for (t in 1:T) {
      y_rep[c, t]   = neg_binomial_2_rng(mu[c, t] + 1e-12, phi[c]);
      log_lik[c, t] = neg_binomial_2_lpmf(y[c, t] | mu[c, t] + 1e-12, phi[c]);
    }
  }
}
"

writeLines(stan_code_joint_lin_hier_het, "outputs/stan_fits/seir_nb_joint_linear_hier_HET.stan")


#writeLines(stan_code_joint_lin_hier_het, "seir_nb_joint_linear_hier_HET.stan")



# =========================================================
# 2) Compile
# =========================================================
mod_lin_hier_het <- stan_model("outputs/stan_fits/seir_nb_joint_linear_hier_HET.stan")

# =========================================================
# 3) Data prep
# =========================================================
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
  
  lockdown_start_EN <- as.Date("2020-03-26")
  lockdown_start_SC <- as.Date("2020-03-24")
  t1_EN <- which.min(abs(dates - lockdown_start_EN))
  t1_SC <- which.min(abs(dates - lockdown_start_SC))
  
  list(
    england = list(
      deaths = as.integer(deaths_EN),
      dates  = dates,
      n_steps = length(dates),
      population = 56e6,
      t1 = t1_EN
    ),
    scotland = list(
      deaths = as.integer(deaths_SC),
      dates  = dates,
      n_steps = length(dates),
      population = 5.5e6,
      t1 = t1_SC
    )
  )
}

joint_data <- prepare_joint_data("chapter06_england_scotland/data/GB_data.csv")
stopifnot(joint_data$england$n_steps == joint_data$scotland$n_steps)

dates <- joint_data$england$dates
T_len <- joint_data$england$n_steps

# fixed process params
rho   <- 0.5
delta <- 1/5.5
gamma <- 1/4

# delay (Wood) params
L <- 60L
lmu_id_fix <- 3.151
lsd_id_fix <- 0.469

IFR <- 0.009

y_mat <- rbind(
  as.integer(joint_data$england$deaths),
  as.integer(joint_data$scotland$deaths)
)
N_vec <- c(joint_data$england$population, joint_data$scotland$population)

# initial condition priors (same pattern you used)
log_E0_mu <- log(c(100, 10))
log_E0_sd <- c(0.1, 0.1)
log_I0_mu <- log_E0_mu - log(2.5)
log_I0_sd <- log_E0_sd

stan_data_common <- list(
  C = 2L,
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
  t1 = as.integer(c(joint_data$england$t1, joint_data$scotland$t1)),
  log_I0_mu = as.numeric(log_I0_mu),
  log_I0_sd = as.numeric(log_I0_sd)
)

# =========================================================
# 4) Fit
# =========================================================
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

fit_het_lin_hier <- fit_one(mod_lin_hier_het, stan_data_common, seed = 20260122)

   rstan::check_hmc_diagnostics(fit_het_lin_hier)

saveRDS(fit_het_lin_hier, "outputs/stan_fits/fit_het_lin_hier.rds")

#
#rstan::check_hmc_diagnostics(fit_joint_lin_hier_het)
#
post_lin_het <- rstan::extract(fit_het_lin_hier)
stopifnot(!is.null(post_lin_het$inc_SE))

# draws x C x T
dim(post_lin_het$inc_SE)

# save arrays for later overlay plotting
saveRDS(
  list(inc_SE = post_lin_het$inc_SE, inc_EI = post_lin_het$inc_EI),
  "outputs/stan_fits/linear_hier_HET_arrays.rds",
  compress = "xz"
)




# =========================================================
# 5) LOO / WAIC
# =========================================================
ll_het <- loo::extract_log_lik(fit_het_lin_hier, parameter_name = "log_lik", merge_chains = FALSE)

r_eff_het <- loo::relative_eff(ll_het)

loo_het  <- loo::loo(ll_het, r_eff = r_eff_het)
waic_het <- loo::waic(ll_het)

cat("\n===== JOINT HIER LINEAR (HET): PSIS-LOO =====\n")
print(loo_het)

cat("\n===== JOINT HIER LINEAR (HET): WAIC =====\n")
print(waic_het)

cat("\n===== Pareto-k summary =====\n")
print(summary(loo_het$diagnostics$pareto_k))

# Save LOO object so compare_npi_specifications.R can build the LOOIC table
saveRDS(loo_het, "outputs/stan_fits/loo_het_lin_hier.rds")




# =========================================================
# POSTERIOR SUMMARY TABLES (FUNCTION-FREE)
#   (A) Population-level summaries (from generated quantities)
#   (B) Country-level summaries (from draws; shared CV repeated)
# =========================================================



post <- rstan::extract(fit_het_lin_hier)

q <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

summ_row <- function(x) {
  qq <- q(x)
  tibble(
    Mean  = mean(x, na.rm = TRUE),
    SD    = sd(x, na.rm = TRUE),
    Q2.5  = qq[[1]],
    Q50   = qq[[2]],
    Q97.5 = qq[[3]]
  )
}

# -----------------------------
# (A) Population-level table
# -----------------------------
# These exist because we included them in generated quantities in the hierarchical linear Stan code:
#   R0_pop_median, R0_pop_mean, w_pop_median, w_pop_mean, cstar_pop_median
# (Plus: CV_shared is a top-level parameter, so we can summarise it too)

pop_list <- list(
  R0_pop_median   = as.numeric(post$R0_pop_median),
  R0_pop_mean     = as.numeric(post$R0_pop_mean),
  w_pop_median    = as.numeric(post$w_pop_median),
  w_pop_mean      = as.numeric(post$w_pop_mean),
  cstar_pop_median= as.numeric(post$cstar_pop_median)
)

# include CV if present (HET model)
if ("CV_shared" %in% names(post)) {
  pop_list$CV_shared <- as.numeric(post$CV_shared)
}

tab_pop <- bind_rows(lapply(names(pop_list), function(nm) {
  summ_row(pop_list[[nm]]) %>% mutate(Parameter = nm)
})) %>%
 dplyr:: select(Parameter, Mean, SD, Q2.5, Q50, Q97.5)

cat("\n== POPULATION-LEVEL POSTERIOR SUMMARIES (HIER LINEAR) ==\n")
print(tab_pop %>% mutate(across(where(is.numeric), ~round(., 4))), n = Inf)

write.csv(tab_pop, "posterior_summaries_pop_joint_linear_hier.csv", row.names = FALSE)

# -----------------------------
# (B) Country-level table
# -----------------------------
# Country parameters available as transformed parameters:
#   R0[c], w[c], t0[c], cstar[c], I0[c], phi[c]
# Shared: CV_shared (if HET)

countries <- c("England", "Scotland")

draws_country <- list(
  "R0"    = post$R0,      # draws x C
  "w"     = post$w,       # draws x C
  "t0"    = post$t0,      # draws x C
  "c*"    = post$cstar,   # draws x C
  "I0"    = post$I0,      # draws x C
  "phi"   = post$phi      # draws x C
)

tab_country <- bind_rows(lapply(1:2, function(c_idx) {
  bind_rows(lapply(names(draws_country), function(par_nm) {
    x <- as.numeric(draws_country[[par_nm]][, c_idx])
    summ_row(x) %>%
      mutate(
        Country = countries[c_idx],
        Parameter = par_nm
      )
  })) %>%
    # add shared CV as an extra "parameter" row for each country (if HET)
    {
      if ("CV_shared" %in% names(post)) {
        bind_rows(
          .,
          summ_row(as.numeric(post$CV_shared)) %>%
            mutate(Country = countries[c_idx], Parameter = "CV (shared)")
        )
      } else .
    }
})) %>%
  select(Country, Parameter, Mean, SD, Q2.5, Q50, Q97.5)

cat("\n== COUNTRY-LEVEL POSTERIOR SUMMARIES (HIER LINEAR) ==\n")
print(tab_country %>% mutate(across(where(is.numeric), ~round(., 4))), n = Inf)

write.csv(tab_country, "posterior_summaries_country_joint_linear_hier.csv", row.names = FALSE)


# =========================================================
# 6) Choose model for plotting
# =========================================================
fit_use     <- fit_het_lin_hier
model_label <- "Joint Linear ramp — Hierarchical — HET"
include_CV  <- TRUE

# include_CV  <- FALSE

post <- rstan::extract(fit_use)

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

country_colors <- c("England" = "#1f77b4", "Scotland" = "#ff7f0e")

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
  labs(title = paste0("Fitted vs Observed — ", model_label),
       subtitle = "Points = observed, line = fitted median, ribbon = 95% NB2 predictive interval",
       x = "Date", y = "Daily deaths") +
  create_custom_theme() +
  theme(legend.position = "bottom") +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p1_fitted)

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

p2a_resid_time <- ggplot(resid_data, aes(x = Date, y = Residual, color = Country)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = country_colors) +
  labs(title = paste0("Residuals vs time — ", model_label),
       subtitle = paste0("RMSE: EN=", round(rmse_EN,2), ", SC=", round(rmse_SC,2)),
       x = "Date", y = "Residual") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y")

print(p2a_resid_time)

# =========================================================
# 9) Trace plots — use your direct array->df method
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
  labs(title = paste0("Trace plots — ", model_label),
       x = "Iteration", y = "Value") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)

# =========================================================
# 10) Posterior densities 
# =========================================================
R0_EN <- as.numeric(post$R0[,1]);  R0_SC <- as.numeric(post$R0[,2])
cs_EN <- as.numeric(post$cstar[,1]); cs_SC <- as.numeric(post$cstar[,2])
w_EN  <- as.numeric(post$w[,1]);    w_SC  <- as.numeric(post$w[,2])
t0_EN <- as.numeric(post$t0[,1]);   t0_SC <- as.numeric(post$t0[,2])
I0_EN <- as.numeric(post$I0[,1]);   I0_SC <- as.numeric(post$I0[,2])
phi_EN_d <- as.numeric(post$phi[,1]); phi_SC_d <- as.numeric(post$phi[,2])
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
  labs(title = paste0("Posterior densities — ", model_label),
       subtitle = "Dashed = posterior median",
       x = "Value", y = "Density") +
  create_custom_theme() +
  facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_dens)

# =========================================================
# 11) Correlation heatmap — joint
# =========================================================
key_mat <- cbind(
  R0_EN = R0_EN, R0_SC = R0_SC,
  cstar_EN = cs_EN, cstar_SC = cs_SC,
  w_EN = w_EN, w_SC = w_SC,
  t0_EN = t0_EN, t0_SC = t0_SC,
  I0_EN = I0_EN, I0_SC = I0_SC,
  phi_EN = phi_EN_d, phi_SC = phi_SC_d
)
if (include_CV && !is.null(CV_draw)) key_mat <- cbind(key_mat, CV = CV_draw)

cor_m <- cor(key_mat, use = "pairwise.complete.obs")
cor_long <- expand.grid(Var1 = colnames(cor_m), Var2 = rownames(cor_m))
cor_long$Correlation <- as.vector(cor_m)

p_cor <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, limits = c(-1,1), name = "Corr") +
  geom_text(aes(label = round(Correlation, 2)), size = 3.3, fontface="bold") +
  labs(title = paste0("Posterior correlation — ", model_label), x = "", y = "") +
  create_custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_cor)

# =========================================================
# 12) Incidence (EI flow): inc_EI = delta * E(t)
# =========================================================
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

# FIX for double vlines: facet variable must match layer variable name/levels
lines_df <- tibble(
  Country = c("England", "Scotland"),
  line    = as.Date(c("2020-03-26", "2020-03-24"))
)

p_inc <- ggplot(inc_sum, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.3) +
  geom_line(aes(y = Med), linewidth = 1.1) +
  geom_vline(data = lines_df, aes(xintercept = line), inherit.aes = FALSE) +
  labs(title = paste0("Simulated incidence — ", model_label),
       x = "Date", y = "delta * E(t)") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_inc)

# =========================================================
# 13) C_t trajectories
# (Preferred: directly from generated quantities ct[c,t])
# =========================================================
ct_EN <- post$ct[,1,]
ct_SC <- post$ct[,2,]

ct_sum <- bind_rows(
  as.data.frame(ct_EN) %>%
    pivot_longer(cols = everything()) %>%
    mutate(Country="England", time = as.integer(gsub("V","", name))) %>%
    transmute(Country, time, ct = value),
  as.data.frame(ct_SC) %>%
    pivot_longer(cols = everything()) %>%
    mutate(Country="Scotland", time = as.integer(gsub("V","", name))) %>%
    transmute(Country, time, ct = value)
) %>%
  group_by(Country, time) %>%
  summarise(Med = median(ct),
            Lo  = quantile(ct, 0.025),
            Hi  = quantile(ct, 0.975),
            .groups="drop") %>%
  mutate(Date = min(dates) + (time - 1))

p_ct <- ggplot(ct_sum, aes(x = Date)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.3) +
  geom_line(aes(y = Med), linewidth = 1.1) +
  geom_vline(data = lines_df, aes(xintercept = line), inherit.aes = FALSE) +
  labs(title = paste0("Simulated NPI profile (C_t) — ", model_label),
       subtitle = "Line = median, ribbon = 95% pointwise interval, vertical = lockdown start",
       x = "Date", y = "C_t") +
  create_custom_theme() +
  facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_ct)

# =========================================================
# OPTIONAL: compute Ct outside Stan (if you want a drop-in fun_get_ct analogue)
# =========================================================
fun_get_ct_linear <- function(post, stan_data) {
  # returns array: draws x C x T
  C <- stan_data$C
  T <- stan_data$T
  t1 <- stan_data$t1
  
  # use transformed parameters if available
  w_draw  <- post$w          # draws x C
  cs_draw <- post$cstar      # draws x C
  t0_draw <- post$t0         # draws x C (already computed as t1-w in Stan)
  
  n_draws <- nrow(w_draw)
  Ct <- array(NA_real_, dim = c(n_draws, C, T))
  
  for (c in 1:C) {
    t1c <- t1[c]
    for (t in 1:T) {
      Ct[,c,t] <- ifelse(
        t < t0_draw[,c], 1.0,
        ifelse(
          t < t1c,
          1.0 - (1.0 - cs_draw[,c]) * pmin(pmax((t - t0_draw[,c]) / pmax(t1c - t0_draw[,c], 1e-6), 0.0), 1.0),
          cs_draw[,c]
        )
      )
    }
  }
  Ct
}

# Example:
# Ct2 <- fun_get_ct_linear(post, stan_data_common)
# all.equal(Ct2[1:10,,], post$ct[1:10,,], tolerance = 1e-6)

