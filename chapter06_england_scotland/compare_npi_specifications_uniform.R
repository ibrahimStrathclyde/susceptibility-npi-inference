# =============================================================================
# compare_npi_specifications_uniform.R
#
# Chapter 6: England and Scotland application
# Comparison and diagnostics across the three NPI specifications with UNIFORM
# (flat) priors and heterogeneous susceptibility.
#
# Parallel to compare_npi_specifications.R (hierarchical priors). Produces the
# same set of figures and tables for the uniform-prior sensitivity analysis.
#
# The three fitted objects are loaded from outputs/stan_fits/ if not already
# in the workspace, so this script can run in a fresh session.
#
# Outputs:
#   outputs/figures/unif_compare3_fitted_vs_observed.png
#   outputs/figures/unif_compare3_fitted_vs_observed_EN.png
#   outputs/figures/unif_compare3_fitted_vs_observed_SC.png
#   outputs/figures/unif_compare3_residuals_vs_time.png
#  
#   outputs/figures/unif_stack_fit_ct_EN.png / _SC.png
#   
#   outputs/results/unif_posterior_summaries_3NPI.csv
#
# Authors: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# =============================================================================

library(rstan)
library(lubridate)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fig_dir <- "outputs/figures"
res_dir <- "outputs/results"
fit_dir <- "outputs/stan_fits"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ============================================================
# Thesis theme
# ============================================================
create_custom_theme <- function(
    title_size = 18, subtitle_size = 14, axis_title_size = 14,
    axis_text_size = 14, legend_title_size = 14, legend_text_size = 14,
    legend_position = "bottom", grid_color_major = "grey90", base_size = 13
) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = title_size,      face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = subtitle_size,   hjust = 0.5),
      axis.title      = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text       = ggplot2::element_text(size = axis_text_size,  face = "bold"),
      legend.title    = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.text     = ggplot2::element_text(size = legend_text_size,  face = "bold"),
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_line(colour = grid_color_major),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "grey60", linewidth = 0.3),
      strip.text       = ggplot2::element_text(size = axis_text_size, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
    )
}

# ============================================================
#  Load uniform prior fits
# ============================================================
fit_files <- c(
  fit_het_lin_unif = file.path(fit_dir, "fit_het_lin_unif.rds"),
  fit_het_log_unif = file.path(fit_dir, "fit_het_log_unif.rds"),
  fit_het_str_unif = file.path(fit_dir, "fit_het_str_unif.rds")
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
      ". Run the corresponding uniform fitting script first."
    )
  }
}

# ============================================================
# Data prep
# ============================================================
prepare_joint_data <- function(filepath = "chapter06_england_scotland/data/GB_data.csv") {
  data <- read.csv(filepath, stringsAsFactors = FALSE)
  data$Date <- lubridate::dmy(data$Date)
  data <- data[order(data$Date), ]

  analysis_start <- as.Date("2020-01-31")
  analysis_end   <- as.Date("2020-06-01")
  idx <- which(data$Date >= analysis_start & data$Date <= analysis_end)

  dates     <- data$Date[idx]
  deaths_EN <- data$Deaths_EN[idx]; deaths_EN[is.na(deaths_EN)] <- 0L
  deaths_SC <- data$Deaths_SC[idx]; deaths_SC[is.na(deaths_SC)] <- 0L

  str_EN <- data$Stringency_EN[idx]; str_EN[is.na(str_EN)] <- 0
  str_SC <- data$Stringency_SC[idx]; str_SC[is.na(str_SC)] <- 0

  norm_flip <- function(x) {
    xmin <- min(x, na.rm = TRUE); xmax <- max(x, na.rm = TRUE)
    if (isTRUE(all.equal(xmax, xmin))) return(rep(1.0, length(x)))
    pmin(pmax(1 - (x - xmin) / (xmax - xmin), 0), 1)
  }

  u_EN <- norm_flip(str_EN); u_SC <- norm_flip(str_SC)

  t1_EN <- which.min(abs(dates - as.Date("2020-03-26")))
  t2_EN <- which.min(abs(dates - as.Date("2020-05-12")))
  t1_SC <- which.min(abs(dates - as.Date("2020-03-24")))
  t2_SC <- which.min(abs(dates - as.Date("2020-05-28")))

  list(
    england = list(deaths = as.integer(deaths_EN), dates = dates,
                   n_steps = length(dates), t1 = t1_EN, t2 = t2_EN,
                   population = 56e6, u = as.numeric(u_EN)),
    scotland = list(deaths = as.integer(deaths_SC), dates = dates,
                    n_steps = length(dates), t1 = t1_SC, t2 = t2_SC,
                    population = 5.5e6, u = as.numeric(u_SC))
  )
}

joint_data <- prepare_joint_data()
dates  <- joint_data$england$dates
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
#  Helpers 
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

summarise_draw_array <- function(arr, dates, model_name,
                                 country_names = c("England","Scotland"),
                                 value_name = "Value") {
  stopifnot(length(dim(arr)) == 3, dim(arr)[2] == 2, dim(arr)[3] == length(dates))
  out <- vector("list", 2)
  for (c in 1:2) {
    mat <- arr[, c, , drop = TRUE]
    out[[c]] <- tibble::tibble(
      Country = country_names[c], Date = dates,
      Med = apply(mat, 2, stats::median),
      Lo  = apply(mat, 2, stats::quantile, probs = 0.025),
      Hi  = apply(mat, 2, stats::quantile, probs = 0.975),
      Model = model_name
    )
  }
  dplyr::bind_rows(out) %>% dplyr::rename(!!value_name := Med)
}

ct_from_fit <- function(fit, model_type, dates, u_mat = NULL, t1_idx = NULL) {
  post <- rstan::extract(fit)
  if (!is.null(post$ct) && length(dim(post$ct)) == 3) return(post$ct)

  ensure_matrix <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) x <- cbind(x, x)
    if (ncol(x) == 1)    x <- cbind(x, x)
    x
  }

  if (model_type == "stringency") {
    cstar <- ensure_matrix(post$cstar)
    n_draw <- nrow(cstar); C <- 2; T <- ncol(u_mat)
    ct <- array(NA_real_, c(n_draw, C, T))
    for (c in 1:C) for (t in 1:T)
      ct[, c, t] <- cstar[, c] + (1 - cstar[, c]) * u_mat[c, t]
    return(ct)
  }

  if (model_type == "logistic") {
    cstar <- ensure_matrix(post$cstar)
    tm    <- ensure_matrix(post$tm)
    k     <- ensure_matrix(post$k)
    n_draw <- nrow(cstar); C <- 2; T <- length(dates)
    ct <- array(NA_real_, c(n_draw, C, T))
    for (c in 1:C) for (t in 1:T)
      ct[, c, t] <- 1 - (1 - cstar[, c]) * inv_logit(k[, c] * (t - tm[, c]))
    return(ct)
  }

  if (model_type == "linear") {
    cstar <- ensure_matrix(post$cstar)
    if (!is.null(post$t0)) {
      t0 <- ensure_matrix(post$t0)
    } else {
      w  <- ensure_matrix(post$w)
      t0 <- sweep(w, 2, t1_idx, FUN = function(wc, t1c) t1c - wc)
    }
    n_draw <- nrow(cstar); C <- 2; T <- length(dates)
    ct <- array(NA_real_, c(n_draw, C, T))
    for (c in 1:C) {
      t1c <- t1_idx[c]
      for (t in 1:T)
        ct[, c, t] <- ifelse(t < t0[, c], 1.0,
          ifelse(t < t1c,
            1 - (1 - cstar[, c]) * pmin(pmax((t - t0[, c]) / pmax(t1c - t0[, c], 1e-6), 0), 1),
            cstar[, c]))
    }
    return(ct)
  }
  stop("Unknown model_type: ", model_type)
}



# ============================================================
#  Model registry
# ============================================================
models <- list(
  list(name = "Linear",     type = "linear",     fit = fit_het_lin_unif),
  list(name = "Logistic",   type = "logistic",   fit = fit_het_log_unif),
  list(name = "Stringency", type = "stringency", fit = fit_het_str_unif)
)

npi_cols <- c("Linear" = "#F8766D", "Logistic" = "#00BA38", "Stringency" = "#619CFF")

# ============================================================
#  Fitted vs observed
# ============================================================
fitted_all <- dplyr::bind_rows(lapply(models, function(m) {
  post    <- rstan::extract(m$fit)
  mu      <- post$mu
  phi_med <- phi_median_by_country(post)

  mu_EN <- apply(mu[,1,], 2, stats::median)
  mu_SC <- apply(mu[,2,], 2, stats::median)

  tibble::tibble(
    Date    = rep(dates, 2),
    Country = rep(c("England","Scotland"), each = n_days),
    Fitted  = c(mu_EN, mu_SC),
    Lower   = c(stats::qnbinom(0.025, size = phi_med[1], mu = mu_EN),
                stats::qnbinom(0.025, size = phi_med[2], mu = mu_SC)),
    Upper   = c(stats::qnbinom(0.975, size = phi_med[1], mu = mu_EN),
                stats::qnbinom(0.975, size = phi_med[2], mu = mu_SC)),
    Model   = m$name
  )
}))

p_fit <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(data = fitted_all,
    ggplot2::aes(x = Date, ymin = Lower, ymax = Upper, fill = Model), alpha = 0.15, colour = NA) +
  ggplot2::geom_line(data = fitted_all,
    ggplot2::aes(x = Date, y = Fitted, colour = Model), linewidth = 1.1) +
  ggplot2::geom_point(data = obs_df,
    ggplot2::aes(x = Date, y = Observed), size = 1, alpha = 0.7) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = Lock), linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = npi_cols) +
  ggplot2::scale_fill_manual(values = npi_cols) +
  ggplot2::labs(title = "Fitted vs observed - uniform prior models",
                x = "Date", y = "Daily deaths", colour = "Model", fill = "Model") +
  create_custom_theme(legend_position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y", ncol = 1)

print(p_fit)

plot_one_country <- function(ctry, lock_date, file_name) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = fitted_all %>% dplyr::filter(Country == ctry),
      ggplot2::aes(x = Date, ymin = Lower, ymax = Upper, fill = Model), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(data = fitted_all %>% dplyr::filter(Country == ctry),
      ggplot2::aes(x = Date, y = Fitted, colour = Model), linewidth = 1.1) +
    ggplot2::geom_point(data = obs_df %>% dplyr::filter(Country == ctry),
      ggplot2::aes(x = Date, y = Observed), size = 1, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = lock_date, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = npi_cols) +
    ggplot2::scale_fill_manual(values = npi_cols) +
    ggplot2::labs(title = paste0("Fitted vs observed - ", ctry, " (uniform priors)"),
                  x = "Date", y = "Daily deaths", colour = "Model", fill = "Model") +
    create_custom_theme(legend_position = "bottom")
  print(p)
  ggplot2::ggsave(file.path(fig_dir, file_name), p, width = 10, height = 4.5, dpi = 300)
  invisible(p)
}

plot_one_country("England",  lock_EN, "unif_compare3_fitted_vs_observed_EN.png")
plot_one_country("Scotland", lock_SC, "unif_compare3_fitted_vs_observed_SC.png")

# ============================================================
# 5) Residuals
# ============================================================
resid_all <- fitted_all %>%
  dplyr::left_join(obs_df, by = c("Date","Country")) %>%
  dplyr::mutate(Residual = Observed - Fitted)

rmse_tbl <- resid_all %>%
  dplyr::group_by(Country, Model) %>%
  dplyr::summarise(RMSE = sqrt(mean(Residual^2)), MAE = mean(abs(Residual)),
                   Bias = mean(Residual), .groups = "drop")

rmse_sub <- rmse_tbl %>%
  dplyr::mutate(txt = paste0(Model, "=", sprintf("%.2f", RMSE))) %>%
  dplyr::group_by(Country) %>%
  dplyr::summarise(line = paste(txt, collapse = ", "), .groups = "drop") %>%
  dplyr::summarise(sub = paste0(Country, ": ", line, collapse = " | "), .groups = "drop") %>%
  dplyr::pull(sub)

p_res_time <- ggplot2::ggplot(resid_all, ggplot2::aes(x = Date, y = Residual, colour = Model)) +
  ggplot2::geom_point(alpha = 0.55, size = 1) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = 0.15) +
  ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = Lock), linewidth = 0.5) +
  ggplot2::scale_colour_manual(values = npi_cols) +
  ggplot2::labs(title = "Residuals vs time - uniform prior models",
                x = "Date", y = "Observed - fitted", colour = "Model") +
  create_custom_theme(legend_position = "bottom") +
  ggplot2::facet_wrap(~Country, scales = "free_y")


print(p_res_time)

#write.csv(rmse_tbl, file.path(res_dir, "unif_compare3_rmse.csv"), row.names = FALSE)

# ============================================================
# Contact modifier c(t)
# ============================================================
ct_all <- dplyr::bind_rows(lapply(models, function(m) {
  ct_arr <- ct_from_fit(m$fit, m$type, dates=dates, u_mat=u_mat, t1_idx=t1_idx)
  summarise_draw_array(ct_arr, dates, m$name, value_name="Ct")
}))

p_ct <- ggplot2::ggplot(ct_all, ggplot2::aes(x=Date, y=Ct, colour=Model, fill=Model)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lo, ymax=Hi), alpha=0.15, colour=NA) +
  ggplot2::geom_line(linewidth=1.1) +
  ggplot2::geom_vline(data=lines_df, ggplot2::aes(xintercept=Lock), linewidth=0.5) +
  ggplot2::scale_colour_manual(values=npi_cols) + ggplot2::scale_fill_manual(values=npi_cols) +
  ggplot2::labs(title="Contact modifier c(t) - uniform prior models",
                subtitle="Line=median; ribbon=95% pointwise interval; vertical=lockdown",
                x="Date", y="c(t)", colour="Model", fill="Model") +
  create_custom_theme(legend_position="bottom") +
  ggplot2::facet_wrap(~Country, scales="fixed", ncol=1) +
  ggplot2::coord_cartesian(ylim=c(0,1.05))

print(p_ct)

# ============================================================
#  Stacked panels per country
# ============================================================

# fitted over c(t)
plot_stack_fit_ct <- function(country) {
  fit_c  <- fitted_all %>% dplyr::filter(Country == country)
  obs_c  <- obs_df     %>% dplyr::filter(Country == country)
  ct_c   <- ct_all     %>% dplyr::filter(Country == country)
  line_c <- lines_df   %>% dplyr::filter(Country == country)
  
  p_fit_c <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data=fit_c, ggplot2::aes(x=Date,ymin=Lower,ymax=Upper,fill=Model), alpha=0.15,colour=NA) +
    ggplot2::geom_line(data=fit_c, ggplot2::aes(x=Date,y=Fitted,colour=Model), linewidth=1.1) +
    ggplot2::geom_point(data=obs_c, ggplot2::aes(x=Date,y=Observed), size=1, alpha=0.7) +
    ggplot2::geom_vline(data=line_c, ggplot2::aes(xintercept=Lock), linewidth=0.5) +
    ggplot2::scale_colour_manual(values=npi_cols) + ggplot2::scale_fill_manual(values=npi_cols) +
    ggplot2::labs(title=paste0("Fitted vs observed : ",country), x=NULL, y="Daily deaths") +
    create_custom_theme(legend_position="bottom") +
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
  
  p_ct_c <- ggplot2::ggplot(ct_c, ggplot2::aes(x=Date,y=Ct,colour=Model,fill=Model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=Lo,ymax=Hi), alpha=0.15,colour=NA) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::geom_vline(data=line_c, ggplot2::aes(xintercept=Lock), linewidth=0.5) +
    ggplot2::coord_cartesian(ylim=c(0,1.05)) +
    ggplot2::scale_colour_manual(values=npi_cols) + ggplot2::scale_fill_manual(values=npi_cols) +
 
    create_custom_theme(legend_position="bottom")
  
  p_s <- (p_fit_c / p_ct_c) + patchwork::plot_layout(heights=c(1,1.05),guides="collect") &
    ggplot2::theme(legend.position="bottom")
  
  sfx <- ifelse(country=="England","EN","SC")
  ggplot2::ggsave(file.path(fig_dir,paste0("unif_stack_fit_ct_",sfx,".png")), p_s, width=10, height=9, dpi=300)
  print(p_s)
  invisible(p_s)
}

plot_stack_fit_ct("England")
plot_stack_fit_ct("Scotland")

# ============================================================
#  Posterior summaries
# ============================================================
posterior_tbl <- function(fit, model_label, pars) {
  keep <- intersect(pars, fit@sim$pars_oi)
  if (length(keep) == 0) return(tibble::tibble(Model=model_label, Parameter=character(),
    Median=numeric(), SD=numeric(), Q2.5=numeric(), Q97.5=numeric()))
  sm <- rstan::summary(fit, pars=keep, probs=c(0.025,0.5,0.975))$summary
  tibble::tibble(Model=model_label, Parameter=rownames(sm),
    Median=sm[,"50%"], SD=sm[,"sd"], Q2.5=sm[,"2.5%"], Q97.5=sm[,"97.5%"])
}

pars_common <- c("R0","cstar","nu","CV_shared","phi","log_I0","I0","w","t0","k","tm")
pars_log    <- c(pars_common, "k_pop_median","tm_pop_median")

post_summ <- dplyr::bind_rows(
  posterior_tbl(fit_het_lin_unif, "Linear",     pars_common),
  posterior_tbl(fit_het_log_unif, "Logistic",   pars_log),
  posterior_tbl(fit_het_str_unif, "Stringency", pars_common)
) %>% dplyr::arrange(Parameter, Model)

print(post_summ, n = Inf)
#write.csv(post_summ, file.path(res_dir,"unif_posterior_summaries_3NPI.csv"), row.names=FALSE)
