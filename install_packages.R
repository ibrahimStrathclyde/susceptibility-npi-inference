# =============================================================================
# install_packages.R
#
# Installs all R packages required to reproduce the analyses in this
# repository. Run this script once on a clean R installation before
# sourcing any chapter script.
#
# RStan requires a working C++ toolchain. Install it first following:
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
# =============================================================================

# --- CRAN packages -----------------------------------------------------------

cran_packages <- c(
  # ODE solving (Chapters 2-5)
  "deSolve",

  # Data manipulation and plotting (all chapters)
  "tidyverse",
  "ggplot2",
  "gridExtra",
  "grid",
  "cowplot",
  "patchwork",

  # Colour scales
  "RColorBrewer",
  "viridis",

  # Statistical utilities (Chapter 2-3 MLE)
  "MASS",

  # Correlation and pairs plots (Chapters 2-3)
  "GGally",

  # Stan interface (Chapter 6)
  "rstan",

  # LOO-CV for Stan model comparison (Chapter 6)
  "loo",

  # MCMC diagnostics (Chapter 6)
  "bayesplot",

  # Date handling (Chapter 6)
  "lubridate",

  # Tabulation helper (Chapters 2-3)
  "knitr"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)

if (length(to_install) > 0) {
  cat("Installing", length(to_install), "package(s):\n")
  cat(paste(" ", to_install, collapse = "\n"), "\n\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages are already installed.\n")
}

# --- RStan configuration -----------------------------------------------------
# Set this once after installing rstan to speed up subsequent compilations.
if (requireNamespace("rstan", quietly = TRUE)) {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  cat("RStan configured: auto_write = TRUE, mc.cores =",
      parallel::detectCores(), "\n")
}

# --- Session information ------------------------------------------------------
cat("\n=== Session information ===\n")
print(sessionInfo())
