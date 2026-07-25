# =============================================================================
# compare_npi_specifications.R
#
# Chapter 6: England and Scotland application
# Comparison and diagnostics across the three NPI functional forms, all with
# heterogeneous susceptibility (Gamma, shared coefficient of variation nu):
#
#   (i)   Linear ramp      fit_het_lin_hier
#   (ii)  Logistic         fit_het_log_hier
#   (iii) Stringency index fit_het_str_hier
#
# Run the three fitting scripts before this one. This script loads their saved
# fits from outputs/stan_fits/ if they are not already in the workspace, so it
# can be run in a fresh R session.
#
# Produces:
#   Fitted vs observed (combined and per country)
#   Residuals vs time
#  
#   Contact modifier c(t) 
#   Stacked two-panel figures per country
#   
#   Posterior summary table across all three models
#
# All outputs go to outputs/figures/ and outputs/results/.
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

# -----------------------
# Libraries
# -----------------------
library(rstan)
library(loo)
library(lubridate)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)

# -----------------------
# Stan + output directories
# -----------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fig_dir <- "outputs/figures"
res_dir <- "outputs/results"
fit_dir <- "outputs/stan_fits"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}
if (!dir.exists(res_dir)) {
  dir.create(res_dir, recursive = TRUE)
}

# -----------------------
# Thesis theme
# -----------------------
create_custom_theme <- function(
    title_size = 18,
    subtitle_size = 14,
    axis_title_size = 14,
    axis_text_size = 14,
    legend_title_size = 14,
    legend_text_size = 14,
    legend_position = "bottom",
    grid_color_major = "grey90",
    base_size = 13
) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = subtitle_size, hjust = 0.5),
      axis.title      = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text       = ggplot2::element_text(size = axis_text_size, face = "bold"),
      legend.title    = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.text     = ggplot2::element_text(size = legend_text_size, face = "bold"),
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_line(colour = grid_color_major),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "grey60", linewidth = 0.3),
      strip.text       = ggplot2::element_text(size = axis_text_size, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
    )
}

# ============================================================
# 0) Load the three heterogeneous fits
#
# Uses the object already in the workspace if present, otherwise reads the
# saved fit from outputs/stan_fits/. This lets the script run either straight
# after the fitting scripts or in a fresh session.
# ============================================================
fit_files <- c(
  fit_het_lin_hier = file.path(fit_dir, "fit_het_lin_hier.rds"),
  fit_het_log_hier = file.path(fit_dir, "fit_het_log_hier.rds"),
  fit_het_str_hier = file.path(fit_dir, "fit_het_str_hier.rds")
)

for (obj_name in names(fit_files)) {
  if (exists(obj_name)) {
    cat("Using", obj_name, "from the workspace\n")
  } else if (file.exists(fit_files[[obj_name]])) {
    cat("Loading", obj_name, "from", fit_files[[obj_name]], "\n")
    assign(obj_name, readRDS(fit_files[[obj_name]]))
  } else {
    stop(
      "Cannot find ", obj_name,
      ". It is neither in the workspace nor saved at ", fit_files[[obj_name]],
      ". Run the corresponding fitting script first."
    )
  }
}

# ============================================================
# 1) Data prep (stringency default)
# ============================================================
prepare_joint_data <- function(filepath = "chapter06_england_scotland/data/GB_data.csv") {
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data$Date <- lubridate::dmy(data$Date)
  data <- data[order(data$Date), ]

  analysis_start <- as.Date("2020-01-31")
  analysis_end   <- as.Date("2020-06-01")
  idx <- which(data$Date >= analysis_start & data$Date <= analysis_end)

  dates <- data$Date[idx]

  deaths_EN <- data$Deaths_EN[idx]; deaths_EN[is.na(deaths_EN)] <- 0L
  deaths_SC <- data$Deaths_SC[idx]; deaths_SC[is.na(deaths_SC)] <- 0L

  str_EN <- data$Stringency_EN[idx]; str_SC <- data$Stringency_SC[idx]
  str_EN[is.na(str_EN)] <- 0
  str_SC[is.na(str_SC)] <- 0

  norm_flip <- function(x) {
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    if (isTRUE(all.equal(xmax, xmin))) return(rep(1.0, length(x)))
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
dates <- joint_data$england$dates
n_days <- length(dates)

lock_EN <- as.Date("2020-03-26")
lock_SC <- as.Date("2020-03-24")

lines_df <- tibble::tibble(
  Country = c("England", "Scotland"),
  Lock    = c(lock_EN, lock_SC)
)

u_mat  <- rbind(joint_data$england$u, joint_data$scotland$u)
t1_idx <- c(joint_data$england$t1, joint_data$scotland$t1)

obs_df <- tibble::tibble(
  Date     = rep(dates, 2),
  Country  = rep(c("England", "Scotland"), each = n_days),
  Observed = c(joint_data$england$deaths, joint_data$scotland$deaths)
)

# ============================================================
# 2) Helpers
# ============================================================
inv_logit <- function(x) plogis(x)

phi_median_by_country <- function(post) {
  if (!is.null(post$phi)) {
    ph <- post$phi
    if (is.null(dim(ph))) return(rep(stats::median(ph), 2))
    if (length(dim(ph)) == 2) return(apply(ph, 2, stats::median))
  }
  rep(NA_real_, 2)
}

summarise_draw_array <- function(arr_draw_C_T, dates, model_name,
                                 country_names = c("England","Scotland"),
                                 value_name = "Value") {
  stopifnot(length(dim(arr_draw_C_T)) == 3)
  stopifnot(dim(arr_draw_C_T)[2] == 2)
  stopifnot(dim(arr_draw_C_T)[3] == length(dates))

  out <- vector("list", 2)
  for (c in 1:2) {
    mat <- arr_draw_C_T[, c, , drop = TRUE]
    out[[c]] <- tibble::tibble(
      Country = country_names[c],
      Date    = dates,
      Med     = apply(mat, 2, stats::median),
      Lo      = apply(mat, 2, stats::quantile, probs = 0.025),
      Hi      = apply(mat, 2, stats::quantile, probs = 0.975),
      Model   = model_name
    )
  }
  dplyr::bind_rows(out) %>%
    dplyr::rename(!!value_name := Med)
}

ct_from_fit <- function(fit, model_type, dates, u_mat = NULL, t1_idx = NULL) {
  post <- rstan::extract(fit)

  if (!is.null(post$ct) && length(dim(post$ct)) == 3) return(post$ct)

  if (model_type == "stringency") {
    if (is.null(u_mat)) stop("ct reconstruction needs u_mat for stringency.")
    if (is.null(post$cstar)) stop("ct reconstruction needs cstar for stringency.")

    cstar <- post$cstar
    if (is.null(dim(cstar))) cstar <- cbind(cstar, cstar)
    if (ncol(cstar) == 1) cstar <- cbind(cstar, cstar)

    n_draw <- nrow(cstar)
    C <- 2
    T <- ncol(u_mat)

    ct <- array(NA_real_, dim = c(n_draw, C, T))
    for (c in 1:C) {
      for (t in 1:T) {
        ct[, c, t] <- cstar[, c] + (1 - cstar[, c]) * u_mat[c, t]
      }
    }
    return(ct)
  }

  if (model_type == "logistic") {
    need <- c("cstar", "tm", "k")
    if (!all(need %in% names(post))) stop("ct reconstruction needs: ", paste(need, collapse=", "))

    cstar <- post$cstar; tm <- post$tm; k <- post$k
    if (is.null(dim(cstar))) cstar <- cbind(cstar, cstar)
    if (is.null(dim(tm)))    tm    <- cbind(tm, tm)
    if (is.null(dim(k)))     k     <- cbind(k, k)
    if (ncol(cstar) == 1) cstar <- cbind(cstar, cstar)
    if (ncol(tm) == 1)    tm    <- cbind(tm, tm)
    if (ncol(k) == 1)     k     <- cbind(k, k)

    n_draw <- nrow(cstar)
    C <- 2
    T <- length(dates)

    ct <- array(NA_real_, dim = c(n_draw, C, T))
    for (c in 1:C) {
      for (t in 1:T) {
        s_t <- inv_logit(k[, c] * (t - tm[, c]))
        ct[, c, t] <- 1 - (1 - cstar[, c]) * s_t
      }
    }
    return(ct)
  }

  if (model_type == "linear") {
    if (is.null(post$cstar)) stop("ct reconstruction needs cstar for linear.")
    if (is.null(t1_idx)) stop("ct reconstruction needs t1_idx for linear.")
    cstar <- post$cstar
    if (is.null(dim(cstar))) cstar <- cbind(cstar, cstar)
    if (ncol(cstar) == 1) cstar <- cbind(cstar, cstar)

    if (!is.null(post$t0)) {
      t0 <- post$t0
      if (is.null(dim(t0))) t0 <- cbind(t0, t0)
      if (ncol(t0) == 1) t0 <- cbind(t0, t0)
    } else if (!is.null(post$w)) {
      w <- post$w
      if (is.null(dim(w))) w <- cbind(w, w)
      if (ncol(w) == 1) w <- cbind(w, w)
      t0 <- sweep(w, 2, t1_idx, FUN = function(wc, t1c) t1c - wc)
    } else {
      stop("ct reconstruction for linear needs t0 or w in posterior.")
    }

    n_draw <- nrow(cstar)
    C <- 2
    T <- length(dates)

    ct <- array(NA_real_, dim = c(n_draw, C, T))
    for (c in 1:C) {
      t1c <- t1_idx[c]
      for (t in 1:T) {
        ct[, c, t] <- ifelse(
          t < t0[, c], 1.0,
          ifelse(
            t < t1c,
            1 - (1 - cstar[, c]) * pmin(pmax((t - t0[, c]) / pmax(t1c - t0[, c], 1e-6), 0), 1),
            cstar[, c]
          )
        )
      }
    }
    return(ct)
  }

  stop("Unknown model_type: ", model_type)
}

peak_prob_tbl <- function(inc_draws, model_name, dates, lock_EN, lock_SC) {
  stopifnot(length(dim(inc_draws)) == 3)
  stopifnot(dim(inc_draws)[2] == 2)
  stopifnot(dim(inc_draws)[3] == length(dates))

  one_country <- function(mat_draw_T, lockdown_date) {
    peak_idx  <- apply(mat_draw_T, 1, which.max)
    peak_date <- dates[peak_idx]

    q_date <- function(p) as.Date(
      stats::quantile(as.numeric(peak_date), p),
      origin = "1970-01-01"
    )

    tibble::tibble(
      peak_median = as.Date(stats::median(as.numeric(peak_date)), origin = "1970-01-01"),
      peak_q025   = q_date(0.025),
      peak_q975   = q_date(0.975),
      p_before    = mean(peak_date < lockdown_date),
      p_on        = mean(peak_date == lockdown_date),
      p_after     = mean(peak_date > lockdown_date)
    )
  }

  dplyr::bind_rows(
    one_country(inc_draws[, 1, ], lock_EN) %>% dplyr::mutate(Model = model_name, Country = "England"),
    one_country(inc_draws[, 2, ], lock_SC) %>% dplyr::mutate(Model = model_name, Country = "Scotland")
  ) %>%
    dplyr::relocate(Model, Country)
}

# ============================================================
# 3) Model registry
# ============================================================
models <- list(
  list(name = "Linear",     type = "linear",     fit = fit_het_lin_hier),
  list(name = "Logistic",   type = "logistic",   fit = fit_het_log_hier),
  list(name = "Stringency", type = "stringency", fit = fit_het_str_hier)
)

npi_cols <- c("Linear" = "#F8766D", "Logistic" = "#00BA38", "Stringency" = "#619CFF")

# ============================================================
# 4) FITTED vs OBSERVED 
# ============================================================
fitted_all <- dplyr::bind_rows(lapply(models, function(m) {
  post <- rstan::extract(m$fit)

  if (is.null(post$mu)) stop("Missing mu in fit: ", m$name)
  mu <- post$mu

  mu_EN_med <- apply(mu[, 1, ], 2, stats::median)
  mu_SC_med <- apply(mu[, 2, ], 2, stats::median)

  phi_med <- phi_median_by_country(post)
  if (anyNA(phi_med)) stop("Missing phi in fit: ", m$name)

  band_EN_lo <- stats::qnbinom(0.025, size = phi_med[1], mu = mu_EN_med)
  band_EN_hi <- stats::qnbinom(0.975, size = phi_med[1], mu = mu_EN_med)

  band_SC_lo <- stats::qnbinom(0.025, size = phi_med[2], mu = mu_SC_med)
  band_SC_hi <- stats::qnbinom(0.975, size = phi_med[2], mu = mu_SC_med)

  tibble::tibble(
    Date    = rep(dates, 2),
    Country = rep(c("England", "Scotland"), each = n_days),
    Fitted  = c(mu_EN_med, mu_SC_med),
    Lower   = c(band_EN_lo, band_SC_lo),
    Upper   = c(band_EN_hi, band_SC_hi),
    Model   = m$name
  )
}))

p_fit_compare <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = fitted_all,
    ggplot2::aes(x = Date, ymin = Lower, ymax = Upper, fill = Model),
    alpha = 0.15, colour = NA
  ) +
  ggplot2::geom_line(
    data = fitted_all,
    ggplot2::aes(x = Date, y = Fitted, colour = Model),
    linewidth = 1.1
  ) +
  ggplot2::geom_point(
    data = obs_df,
    ggplot2::aes(x = Date, y = Observed),
    size = 1, alpha = 0.7
  ) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = Lock), linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = npi_cols) +
  ggplot2::scale_fill_manual(values = npi_cols) +
  ggplot2::labs(
    title = "Fitted vs observed - joint hierarchical models",
    x = "Date", y = "Daily deaths", colour = "Model", fill = "Model"
  ) +
  create_custom_theme(legend_position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_fit_compare)
ggplot2::ggsave(
  filename = file.path(fig_dir, "compare3_fitted_vs_observed.png"),
  plot = p_fit_compare, width = 10, height = 7, dpi = 300
)

plot_one_country <- function(ctry, lock_date, file_name) {

  fitted_cty <- fitted_all %>% dplyr::filter(Country == ctry)
  obs_cty    <- obs_df    %>% dplyr::filter(Country == ctry)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = fitted_cty,
      ggplot2::aes(x = Date, ymin = Lower, ymax = Upper, fill = Model),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(
      data = fitted_cty,
      ggplot2::aes(x = Date, y = Fitted, colour = Model),
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = obs_cty,
      ggplot2::aes(x = Date, y = Observed),
      size = 1, alpha = 0.7
    ) +
    ggplot2::geom_vline(xintercept = lock_date, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = npi_cols) +
    ggplot2::scale_fill_manual(values = npi_cols) +
    ggplot2::labs(
      title = paste0("Fitted vs observed - ", ctry),
      x = "Date", y = "Daily deaths", colour = "Model", fill = "Model"
    ) +
    create_custom_theme(legend_position = "bottom")

  print(p)
  ggplot2::ggsave(file.path(fig_dir, file_name), p, width = 10, height = 4.5, dpi = 300)

  invisible(p)
}

plot_one_country(
  ctry = "England",
  lock_date = as.Date(lines_df$Lock[lines_df$Country == "England"]),
  file_name = "compare3_fitted_vs_observed_EN.png"
)

plot_one_country(
  ctry = "Scotland",
  lock_date = as.Date(lines_df$Lock[lines_df$Country == "Scotland"]),
  file_name = "compare3_fitted_vs_observed_SC.png"
)

# ============================================================
# 5) Residuals: vs time  
# ============================================================
resid_all <- fitted_all %>%
  dplyr::left_join(obs_df, by = c("Date", "Country")) %>%
  dplyr::mutate(Residual = Observed - Fitted)

rmse_tbl <- resid_all %>%
  dplyr::group_by(Country, Model) %>%
  dplyr::summarise(
    RMSE = sqrt(mean(Residual^2)),
    MAE  = mean(abs(Residual)),
    Bias = mean(Residual),
    .groups = "drop"
  )



p_res_time <- ggplot2::ggplot(resid_all, ggplot2::aes(x = Date, y = Residual, colour = Model)) +
  ggplot2::geom_point(alpha = 0.55, size = 1) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = 0.15) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = Lock), linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = npi_cols) +
  ggplot2::labs(
    title = "Residuals vs time - joint hierarchical models",   
    x = "Date", y = "Observed - fitted", colour = "Model"
  ) +
  create_custom_theme(legend_position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y")


print(p_res_time)


# ============================================================
#  Contact modifier overlay: c(t)
# ============================================================
ct_all <- dplyr::bind_rows(lapply(models, function(m) {
  ct_arr <- ct_from_fit(m$fit, m$type, dates = dates, u_mat = u_mat, t1_idx = t1_idx)
  summarise_draw_array(ct_arr, dates, model_name = m$name, value_name = "Ct")
}))

p_ct_compare <- ggplot2::ggplot(ct_all,
                                ggplot2::aes(x = Date, y = Ct, colour = Model, fill = Model)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lo, ymax = Hi), alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = Lock), linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = npi_cols) +
  ggplot2::scale_fill_manual(values = npi_cols) +
  ggplot2::labs(
    title = "Contact modifier c(t) - joint hierarchical models",
    subtitle = "Lines = posterior medians; ribbons = 95% pointwise intervals; vertical = lockdown date",
    x = "Date", y = "c(t)", colour = "Model", fill = "Model"
  ) +
  create_custom_theme(legend_position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "fixed", ncol = 1) +
  ggplot2::coord_cartesian(ylim = c(0, 1.05))

print(p_ct_compare)
ggplot2::ggsave(file.path(fig_dir, "compare3_ct.png"),
                p_ct_compare, width = 10, height = 7, dpi = 300)


# ============================================================
# Stacked panels per country
# ============================================================


# ---  fitted over c(t) ---
plot_stacked_fit_ct <- function(country, file_stub) {
  fit_c  <- fitted_all %>% dplyr::filter(Country == country)
  obs_c  <- obs_df     %>% dplyr::filter(Country == country)
  ct_c   <- ct_all     %>% dplyr::filter(Country == country)
  line_c <- lines_df   %>% dplyr::filter(Country == country)

  p_fit <- ggplot() +
    geom_ribbon(data = fit_c, aes(x = Date, ymin = Lower, ymax = Upper, fill = Model),
                alpha = 0.15, colour = NA) +
    geom_line(data = fit_c, aes(x = Date, y = Fitted, colour = Model), linewidth = 1.1) +
    geom_point(data = obs_c, aes(x = Date, y = Observed), size = 1, alpha = 0.7) +
    geom_vline(data = line_c, aes(xintercept = Lock), linewidth = 0.5) +
    scale_colour_manual(values = npi_cols) +
    scale_fill_manual(values = npi_cols) +
    labs(title = paste0("Fitted vs observed - ", country),
         x = NULL, y = "Daily deaths", colour = "Model", fill = "Model") +
    create_custom_theme(legend_position = "bottom") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_ct <- ggplot(ct_c, aes(x = Date, y = Ct, colour = Model, fill = Model)) +
    geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_vline(data = line_c, aes(xintercept = Lock), linewidth = 0.5) +
    coord_cartesian(ylim = c(0, 1.05)) +
    scale_colour_manual(values = npi_cols) +
    scale_fill_manual(values = npi_cols) +
    labs(title = paste0("Contact modifier c(t) - ", country),
         subtitle = "Line = posterior median; ribbon = 95% pointwise interval; vertical = lockdown date",
         x = "Date", y = "c(t)", colour = "Model", fill = "Model") +
    create_custom_theme(legend_position = "bottom")

  p_stack <- (p_fit / p_ct) +
    plot_layout(heights = c(1, 1.05), guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, paste0(file_stub, "_",
                                   ifelse(country == "England", "EN", "SC"), ".png")),
         p_stack, width = 10, height = 9, dpi = 300)
  p_stack
}

p_EN_fit_ct <- plot_stacked_fit_ct("England",  "stack_fit_ct")
p_SC_fit_ct <- plot_stacked_fit_ct("Scotland", "stack_fit_ct")
print(p_EN_fit_ct); print(p_SC_fit_ct)


# ============================================================
#  Posterior summaries across all three models
# ============================================================
posterior_tbl <- function(fit, model_label, pars) {
  avail <- fit@sim$pars_oi
  keep  <- intersect(pars, avail)
  if (length(keep) == 0) {
    return(tibble::tibble(Model = model_label, Parameter = character(),
                          Median = numeric(), SD = numeric(),
                          Q2.5 = numeric(), Q97.5 = numeric()))
  }

  sm <- rstan::summary(fit, pars = keep, probs = c(0.025, 0.5, 0.975))$summary

  tibble::tibble(
    Model     = model_label,
    Parameter = rownames(sm),
    Median    = sm[, "50%"],
    SD        = sm[, "sd"],
    Q2.5      = sm[, "2.5%"],
    Q97.5     = sm[, "97.5%"]
  )
}

pars_hier_common <- c(
  "mu_log_R0", "sigma_log_R0", "z_log_R0", "R0",
  "mu_logit_cstar", "sigma_logit_cstar", "z_logit_cstar", "logit_cstar", "cstar",
  "nu", "phi", "log_I0", "I0",
  "R0_pop_median", "R0_pop_mean", "cstar_pop_median"
)

pars_linear_add   <- c("t0", "w")
pars_logistic_add <- c("k", "tm", "mu_tm", "sigma_tm", "z_tm",
                       "mu_log_k", "sigma_log_k", "z_log_k",
                       "k_pop_median", "k_pop_mean",
                       "tm_pop_median", "tm_pop_mean")
pars_string_add   <- character(0)

post_summ_all <- dplyr::bind_rows(
  posterior_tbl(fit_het_lin_hier, "Linear",     c(pars_hier_common, pars_linear_add)),
  posterior_tbl(fit_het_log_hier, "Logistic",   c(pars_hier_common, pars_logistic_add)),
  posterior_tbl(fit_het_str_hier, "Stringency", c(pars_hier_common, pars_string_add))
) %>%
  dplyr::arrange(Parameter, Model)

print(post_summ_all, n = Inf)
write.csv(post_summ_all, file.path(res_dir, "posterior_summaries_3NPI.csv"), row.names = FALSE)
saveRDS(post_summ_all, file.path(res_dir, "posterior_summaries_3NPI.rds"))

      
