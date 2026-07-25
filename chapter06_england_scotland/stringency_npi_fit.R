
# =============================================================================
# stringency_npi_fit.R
#
# Chapter 6: England and Scotland application
# Joint hierarchical fit with a STRINGENCY-INDEX driven NPI,
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu).
#
# The contact modifier is driven directly by the OxCGRT stringency index:
#   c(t) = c* + (1 - c*) * u(t)
# where u(t) is the country-specific stringency index rescaled to [0, 1] and
# flipped, so u = 1 at minimum stringency and u = 0 at maximum stringency.
# Unlike the linear and logistic specifications, the shape of c(t) is fixed by
# the data; only its depth (c*) is estimated.
#
# England and Scotland are fitted jointly. R0 and c* are partially pooled;
# nu is shared across countries.
#
# Data: chapter06_england_scotland/data/GB_data.csv
#       31 Jan 2020 to 01 Jun 2020 (first wave)
#
# Outputs:
#   outputs/stan_fits/seir_nb_joint_hier_stringency_HET.stan  Stan model
#   outputs/stan_fits/fit_het_str_hier.rds                    fitted object
#   outputs/stan_fits/loo_het_str_hier.rds                    LOO object
#   outputs/figures/                                           diagnostic plots
#
# The fitted object is named fit_het_str_hier. This name is required by
# compare_npi_specifications.R. Do not rename it.
#
# Sampler: 4 chains, 2500 iterations, 1250 warmup,
#          adapt_delta = 0.98, max_treedepth = 15.
# Run time: roughly 1 to 3 hours depending on hardware.
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

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

# ------------------------------------------------------------
# Thesis theme 
# ------------------------------------------------------------
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
      axis.line        = element_line(colour = "grey60", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(size = axis_text_size, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

# ============================================================
# ===================== STAN: HET (HIER) =====================
# ============================================================
stan_code_joint_str_het <- "
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

  // stringency-based NPI input (already flipped+normalised to [0,1])
  matrix<lower=0, upper=1>[C, T] u;

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

  real  mu_logit_cstar;
  real<lower=0> sigma_logit_cstar;

  // non-centred country effects
  vector[C] z_log_R0;
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
  vector<lower=0,upper=1>[C] cstar =
    inv_logit(mu_logit_cstar + sigma_logit_cstar * z_logit_cstar);

  vector<lower=0>[C] I0 = exp(log_I0);
  vector<lower=0>[C] E0 = 2.5 * I0;

  real CV2 = square(CV_shared);

  vector[L] p_delay = delay_probs(L, lmu_id_fix, lsd_id_fix);

  matrix[C, T] mu;
  matrix<lower=0>[C, T] inc_EI;

  // infections (S -> E), stored at time index t (1..T)
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
      // stringency mapping: c_t in [c*,1]
      real c_t = cstar[c] + (1.0 - cstar[c]) * u[c, t];

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

  mu_logit_cstar    ~ normal(logit(0.33), 0.7);
  sigma_logit_cstar ~ normal(0, 0.7);

  z_log_R0 ~ normal(0, 1);
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

  real cstar_pop_median = inv_logit(mu_logit_cstar);

  for (c in 1:C) {
    for (t in 1:T) {
      ct[c,t] = cstar[c] + (1.0 - cstar[c]) * u[c, t];
      y_rep[c,t]   = neg_binomial_2_rng(mu[c,t] + 1e-12, phi[c]);
      log_lik[c,t] = neg_binomial_2_lpmf(y[c,t] | mu[c,t] + 1e-12, phi[c]);
    }
  }
}
"
writeLines(stan_code_joint_str_het, "outputs/stan_fits/seir_nb_joint_hier_stringency_HET.stan")



# ============================================================
# Data prep (stringency-based u_t) 
# ============================================================
prepare_joint_data <- function(filepath = "chapter06_england_scotland/data/GB_data.csv") {
  cat("Preparing joint data for England and Scotland (with stringency-based u_t)...\n")
  
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data$Date <- dmy(data$Date)
  data <- data[order(data$Date), ]
  
  analysis_start <- as.Date("2020-01-31")
  analysis_end   <- as.Date("2020-06-01")
  idx <- which(data$Date >= analysis_start & data$Date <= analysis_end)
  
  dates <- data$Date[idx]
  
  deaths_EN <- data$Deaths_EN[idx]; deaths_EN[is.na(deaths_EN)] <- 0L
  deaths_SC <- data$Deaths_SC[idx]; deaths_SC[is.na(deaths_SC)] <- 0L
  
  str_EN <- data$Stringency_EN[idx]
  str_SC <- data$Stringency_SC[idx]
  
  str_EN[is.na(str_EN)] <- 0
  str_SC[is.na(str_SC)] <- 0
  
  norm_flip <- function(x) {
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    if (isTRUE(all.equal(xmax, xmin))) {
      # constant series: set u_t = 1 (no variation) after flip/normalise
      return(rep(1.0, length(x)))
    }
    r <- (x - xmin) / (xmax - xmin)
    u <- 1 - r
    pmin(pmax(u, 0), 1)
  }
  
  u_EN <- norm_flip(str_EN)
  u_SC <- norm_flip(str_SC)
  
  lockdown_start_EN <- as.Date("2020-03-26")
  lockdown_end_EN   <- as.Date("2020-05-12")
  lockdown_start_SC <- as.Date("2020-03-24")
  lockdown_end_SC   <- as.Date("2020-05-28")
  
  t1_EN <- which.min(abs(dates - lockdown_start_EN))
  t2_EN <- which.min(abs(dates - lockdown_end_EN))
  t1_SC <- which.min(abs(dates - lockdown_start_SC))
  t2_SC <- which.min(abs(dates - lockdown_end_SC))
  
  cat(sprintf("Using %s to %s; %d days\n", min(dates), max(dates), length(dates)))
  
  list(
    england = list(
      deaths = as.integer(deaths_EN),
      dates  = dates,
      n_steps = length(dates),
      t1 = t1_EN,
      t2 = t2_EN,
      population = 56e6,
      u = as.numeric(u_EN)
    ),
    scotland = list(
      deaths = as.integer(deaths_SC),
      dates  = dates,
      n_steps = length(dates),
      t1 = t1_SC,
      t2 = t2_SC,
      population = 5.5e6,
      u = as.numeric(u_SC)
    )
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

u_mat <- rbind(
  as.numeric(joint_data$england$u),
  as.numeric(joint_data$scotland$u)
)

# ============================================================
# Fixed process params 
# ============================================================
rho   <- 0.5
delta <- 1/5.5
gamma <- 1/4
L     <- 60L

IFR <- 0.009

# Wood delay
lmu_id_fix <- 3.151
lsd_id_fix <- 0.469

# I0 prior construction 
log_E0_mu <- log(c(100, 10))
log_E0_sd <- c(0.1, 0.1)

log_I0_mu <- log_E0_mu - log(2.5)
log_I0_sd <- log_E0_sd

stan_data_common <- list(
  C = C,
  T = as.integer(T_len),
  y = y_mat,
  N = as.numeric(N_vec),
  u = u_mat,
  rho = as.numeric(rho),
  delta = as.numeric(delta),
  gamma = as.numeric(gamma),
  L = as.integer(L),
  lmu_id_fix = as.numeric(lmu_id_fix),
  lsd_id_fix = as.numeric(lsd_id_fix),
  IFR = as.numeric(IFR),
  log_I0_mu  = as.numeric(log_I0_mu),
  log_I0_sd  = as.numeric(log_I0_sd),
  CV_shared_fixed = 0.0
)

dates <- joint_data$england$dates
n_days <- length(dates)

# ============================================================
# Compile
# ============================================================
mod_het_str_hier <- stan_model("outputs/stan_fits/seir_nb_joint_hier_stringency_HET.stan")


# ---- Delay sensitivity ----

lmu_id_fix <- 3.014
lsd_id_fix <- 0.469

stan_data_common$lmu_id_fix <- as.numeric(lmu_id_fix)
stan_data_common$lsd_id_fix <- as.numeric(lsd_id_fix)


fit_het_str_hier <- sampling(
  object  = mod_het_str_hier,
  data    = stan_data_common,
  seed    = 16012026,
  chains  = 4,
  iter    = 2500,
  warmup  = 1250,
  control = list(adapt_delta = 0.98, max_treedepth = 15)
)






rstan::check_hmc_diagnostics(fit_het_str_hier)

saveRDS(fit_het_str_hier, "outputs/stan_fits/fit_het_str_hier.rds")


# =========================================================
# DIAGNOSTICS #   - Joint plots, 
#   
# =========================================================

# -----------------------
# Choose fit to diagnose
# -----------------------
fit_use     <- fit_het_str_hier
model_label <- "Joint hierarchical stringency (HET)"
include_CV  <- F

# include_CV  <- FALSE

post <- rstan::extract(fit_use)

# -----------------------
# Helpers (template style)
# -----------------------
q <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

summ_row <- function(x) {
  qq <- q(x)
  tibble::tibble(
    Mean  = mean(x, na.rm = TRUE),
    SD    = stats::sd(x, na.rm = TRUE),
    Q2.5  = qq[[1]],
    Q50   = qq[[2]],
    Q97.5 = qq[[3]]
  )
}

# -----------------------------
# (A) Population-level table
# -----------------------------
# 
#   R0_pop_median, R0_pop_mean, cstar_pop_median (+ CV_shared if HET)

pop_list <- list(
  R0_pop_median    = as.numeric(post$R0_pop_median),
  R0_pop_mean      = as.numeric(post$R0_pop_mean),
  cstar_pop_median = as.numeric(post$cstar_pop_median)
)

if (include_CV && ("CV_shared" %in% names(post))) {
  pop_list$CV_shared <- as.numeric(post$CV_shared)
}

tab_pop <- dplyr::bind_rows(lapply(names(pop_list), function(nm) {
  summ_row(pop_list[[nm]]) %>% dplyr::mutate(Parameter = nm)
})) %>%
  dplyr::select(Parameter, Mean, SD, Q2.5, Q50, Q97.5)

cat("\n== POPULATION-LEVEL POSTERIOR SUMMARIES (STRINGENCY HIER, WEAK PRIORS) ==\n")
print(tab_pop %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4))), n = Inf)

# Optional CSV
# write.csv(tab_pop, "posterior_summaries_pop_joint_stringency_hier_weak.csv", row.names = FALSE)

# -----------------------------
# (B) Country-level table
# -----------------------------
countries <- c("England", "Scotland")

# Prefer I0 from transformed parameters; fallback: exp(log_I0)
I0_draws <- if ("I0" %in% names(post)) post$I0 else exp(post$log_I0)

draws_country <- list(
  "R0"    = post$R0,      # draws x C
  "c*"    = post$cstar,   # draws x C
  "I0"    = I0_draws,     # draws x C
  "phi"   = post$phi      # draws x C
)

tab_country <- dplyr::bind_rows(lapply(1:2, function(c_idx) {
  dplyr::bind_rows(lapply(names(draws_country), function(par_nm) {
    x <- as.numeric(draws_country[[par_nm]][, c_idx])
    summ_row(x) %>%
      dplyr::mutate(Country = countries[c_idx], Parameter = par_nm)
  })) %>%
    {
      if (include_CV && ("CV_shared" %in% names(post))) {
        dplyr::bind_rows(
          .,
          summ_row(as.numeric(post$CV_shared)) %>%
            dplyr::mutate(Country = countries[c_idx], Parameter = "CV (shared)")
        )
      } else .
    }
})) %>%
  dplyr::select(Country, Parameter, Mean, SD, Q2.5, Q50, Q97.5)

cat("\n== COUNTRY-LEVEL POSTERIOR SUMMARIES (STRINGENCY HIER, WEAK PRIORS) ==\n")
print(tab_country %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 4))), n = Inf)

# Optional CSV
# write.csv(tab_country, "posterior_summaries_country_joint_stringency_hier_weak.csv", row.names = FALSE)

# =========================================================
# Fitted vs observed — JOINT (facetted) + NB predictive bands
# =========================================================
dates  <- joint_data$england$dates
n_days <- length(dates)

y_EN <- as.numeric(joint_data$england$deaths)
y_SC <- as.numeric(joint_data$scotland$deaths)

mu_EN <- apply(post$mu[,1,], 2, stats::median)
mu_SC <- apply(post$mu[,2,], 2, stats::median)

phi_EN <- stats::median(post$phi[,1])
phi_SC <- stats::median(post$phi[,2])

band_en <- list(
  lo = stats::qnbinom(0.025, size = phi_EN, mu = mu_EN),
  hi = stats::qnbinom(0.975, size = phi_EN, mu = mu_EN)
)
band_sc <- list(
  lo = stats::qnbinom(0.025, size = phi_SC, mu = mu_SC),
  hi = stats::qnbinom(0.975, size = phi_SC, mu = mu_SC)
)

fitted_data <- data.frame(
  Date     = rep(dates, 2),
  Country  = rep(c("England","Scotland"), each = n_days),
  Observed = c(y_EN, y_SC),
  Fitted   = c(mu_EN, mu_SC),
  Lower    = c(band_en$lo, band_sc$lo),
  Upper    = c(band_en$hi, band_sc$hi)
)

country_colors <- c("England"="#1f77b4","Scotland"="#ff7f0e")

lines_df <- tibble::tibble(
  Country = c("England","Scotland"),
  line    = as.Date(c("2020-03-26","2020-03-24"))
)

p_fitted <- ggplot2::ggplot(fitted_data, ggplot2::aes(x = Date)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper, fill = Country), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = Fitted, color = Country), linewidth = 1.2) +
  ggplot2::geom_point(ggplot2::aes(y = Observed, color = Country), size = 1, alpha = 0.7) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = line),
                      inherit.aes = FALSE, linewidth = 0.5) +
  ggplot2::scale_color_manual(values = country_colors) +
  ggplot2::scale_fill_manual(values = country_colors) +
  ggplot2::labs(
    title = paste0("Fitted vs Observed: ", model_label),   
    x = "Date", y = "Daily deaths"
  ) +
  create_custom_theme() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_fitted)

# =========================================================
# Residuals — JOINT (facetted): vs time
# =========================================================
res_EN <- y_EN - mu_EN
res_SC <- y_SC - mu_SC

rmse_EN <- sqrt(mean(res_EN^2))
rmse_SC <- sqrt(mean(res_SC^2))

cat(sprintf("\n=== RMSE (%s) ===\n", model_label))
cat(sprintf("England RMSE: %.2f\n", rmse_EN))
cat(sprintf("Scotland RMSE: %.2f\n", rmse_SC))

resid_data <- data.frame(
  Date     = rep(dates, 2),
  Country  = rep(c("England","Scotland"), each = n_days),
  Residual = c(res_EN, res_SC),
  Fitted   = c(mu_EN, mu_SC)
)

p_res_time <- ggplot2::ggplot(resid_data, ggplot2::aes(x = Date, y = Residual, color = Country)) +
  ggplot2::geom_point(alpha = 0.7, size = 1) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = line),
                      inherit.aes = FALSE, linewidth = 0.5) +
  ggplot2::scale_color_manual(values = country_colors) +
  ggplot2::labs(
    title = paste0("Residuals vs time : ", model_label),
    subtitle = paste0("RMSE: EN=", round(rmse_EN,2), ", SC=", round(rmse_SC,2)),
    x = "Date", y = "Observed − fitted"
  ) +
  create_custom_theme() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y")

print(p_res_time)


# =========================================================
# Trace plots 
# =========================================================
params_to_trace <- c(
  "R0[1]","R0[2]",
  "cstar[1]","cstar[2]",
  "log_I0[1]","log_I0[2]",
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
) %>%
  dplyr::mutate(Parameter = ifelse(as.character(Parameter) == "CV_shared", "CV", as.character(Parameter))) %>%
  dplyr::mutate(Parameter = factor(Parameter, levels = unique(Parameter)))

p_trace <- ggplot2::ggplot(trace_df, ggplot2::aes(x = Iteration, y = Value, color = Chain)) +
  ggplot2::geom_line(alpha = 0.75, linewidth = 0.35) +
  ggplot2::labs(title = paste0("Trace plots — ", model_label),
                x = "Iteration", y = "Value") +
  create_custom_theme() +
  ggplot2::facet_wrap(~Parameter, scales = "free_y", ncol = 3)

print(p_trace)

# =========================================================
# Posterior densities — same set as trace (joint)
# =========================================================
R0_EN <- as.numeric(post$R0[,1]);    R0_SC <- as.numeric(post$R0[,2])
cs_EN <- as.numeric(post$cstar[,1]); cs_SC <- as.numeric(post$cstar[,2])
I0_EN <- as.numeric(I0_draws[,1]);   I0_SC <- as.numeric(I0_draws[,2])
ph_EN <- as.numeric(post$phi[,1]);   ph_SC <- as.numeric(post$phi[,2])

draw_list <- list(
  "R0[1]" = R0_EN,  "R0[2]" = R0_SC,
  "cstar[1]" = cs_EN, "cstar[2]" = cs_SC,
  "I0[1]" = I0_EN,  "I0[2]" = I0_SC,
  "phi[1]" = ph_EN, "phi[2]" = ph_SC
)
if (include_CV && ("CV_shared" %in% names(post))) draw_list[["CV"]] <- as.numeric(post$CV_shared)

post_dist <- dplyr::bind_rows(lapply(names(draw_list), function(nm) {
  tibble::tibble(Parameter = nm, Value = as.numeric(draw_list[[nm]]))
})) %>%
  dplyr::mutate(Parameter = factor(Parameter, levels = names(draw_list)))

meds <- post_dist %>% dplyr::group_by(Parameter) %>% dplyr::summarise(med = median(Value), .groups="drop")

p_dens <- ggplot2::ggplot(post_dist, ggplot2::aes(x = Value)) +
  ggplot2::geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.6) +
  ggplot2::geom_vline(data = meds, ggplot2::aes(xintercept = med),
                      linetype = "dashed", linewidth = 1) +
  ggplot2::labs(title = paste0("Posterior densities — ", model_label),
                subtitle = "Dashed = posterior median",
                x = "Value", y = "Density") +
  create_custom_theme() +
  ggplot2::facet_wrap(~Parameter, scales = "free", ncol = 3)

print(p_dens)

# =========================================================
# Incidence (EI flow): inc_EI = delta * E(t)  (saved in TP)
# =========================================================
inc_EN <- post$inc_EI[,1,]
inc_SC <- post$inc_EI[,2,]

inc_sum <- dplyr::bind_rows(
  as.data.frame(inc_EN) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    dplyr::mutate(Country="England", time = as.integer(gsub("V","", name))) %>%
    dplyr::transmute(Country, time, inc = value),
  as.data.frame(inc_SC) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    dplyr::mutate(Country="Scotland", time = as.integer(gsub("V","", name))) %>%
    dplyr::transmute(Country, time, inc = value)
) %>%
  dplyr::group_by(Country, time) %>%
  dplyr::summarise(
    Med = median(inc),
    Lo  = stats::quantile(inc, 0.025),
    Hi  = stats::quantile(inc, 0.975),
    .groups="drop"
  ) %>%
  dplyr::mutate(Date = min(dates) + (time - 1))

p_inc <- ggplot2::ggplot(inc_sum, ggplot2::aes(x = Date)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lo, ymax = Hi), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = Med), linewidth = 1.1) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = line),
                      inherit.aes = FALSE, linewidth = 0.5) +
  ggplot2::labs(title = paste0("Simulated incidence  ", model_label),
                subtitle = "Line = median, ribbon = 95% pointwise interval; vertical = lockdown",
                x = "Date", y = "reconstructed incidence") +
  create_custom_theme() +
  ggplot2::facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_inc)

# =========================================================
# c(t) trajectories — USE generated quantities ct
# =========================================================
ct_EN <- post$ct[,1,]
ct_SC <- post$ct[,2,]

ct_sum <- dplyr::bind_rows(
  as.data.frame(ct_EN) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    dplyr::mutate(Country="England", time = as.integer(gsub("V","", name))) %>%
    dplyr::transmute(Country, time, ct = value),
  as.data.frame(ct_SC) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    dplyr::mutate(Country="Scotland", time = as.integer(gsub("V","", name))) %>%
    dplyr::transmute(Country, time, ct = value)
) %>%
  dplyr::group_by(Country, time) %>%
  dplyr::summarise(
    Med = median(ct),
    Lo  = stats::quantile(ct, 0.025),
    Hi  = stats::quantile(ct, 0.975),
    .groups="drop"
  ) %>%
  dplyr::mutate(Date = min(dates) + (time - 1))

p_ct <- ggplot2::ggplot(ct_sum, ggplot2::aes(x = Date)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lo, ymax = Hi), alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(y = Med), linewidth = 1.1) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = line),
                      inherit.aes = FALSE, linewidth = 0.5) +
  ggplot2::labs(title = paste0("Contact modifier profile c(t) — ", model_label),
                subtitle = "Line = median, ribbon = 95% pointwise interval; vertical = lockdown",
                x = "Date", y = "c(t)") +
  create_custom_theme() +
  ggplot2::facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_ct)

















