# susceptibility-npi-inference

**Simultaneous inference of susceptibility distributions and non-pharmaceutical
interventions from epidemic trajectories**

Ibrahim Mohammed · Chris Robertson · M. Gabriela M. Gomes  
University of Strathclyde, Glasgow

---

> **Status:** Repository complete. All analysis scripts have been cleaned,
> documented, and tested for reproducibility. See the [Status](#status)
> section for the full component checklist.

---

## Overview

This repository contains the R and Stan code supporting the PhD thesis:

> Mohammed, I. (2026). *Simultaneous inference of susceptibility distributions
> and non-pharmaceutical interventions from epidemic trajectories.*
> University of Strathclyde, Glasgow. Student ID 202279132.

The thesis develops and applies a heterogeneous SEIR modelling framework in
which susceptibility varies across individuals according to a Gamma distribution.
Non-pharmaceutical interventions (NPIs) are incorporated through a
contact-reduction factor c(t). The central identifiability challenge is the
compensation ridge between the susceptibility coefficient of variation (nu) and
the NPI contact-reduction parameter (c*): the two parameters trade off along a
near-flat likelihood surface when only a single epidemic is observed. The thesis
shows that joint fitting to two epidemics with different initial conditions
resolves this ridge.

The empirical application fits the heterogeneous SEIR model under three NPI
functional forms (piecewise-linear ramp, logistic, and stringency-index driven)
to first-wave COVID-19 mortality data from England and Scotland, using both
hierarchical and uniform prior specifications.

---

## Repository structure

```
susceptibility-npi-inference/
|
|- README.md
|- LICENSE
|- .gitignore
|- install_packages.R                    # installs all required R packages
|
|- R/                                    # shared function libraries
|   |- utility_functions.R              # ODE (Reduced.m_intervene), simulation,
|   |                                   #   logit/expit  [Chs. 2-3]
|   |- mle_functions.R                  # Poisson likelihoods, MLE fitters (het and hom),
|   |                                   #   profile likelihood, generate_trajectory  [Chs. 2-3, 5]
|   `- Generalfun_mispec_distribution.R # LA-SM discretisation, K-group ODE,
|                                       #   Gamma and Lognormal fitters  [Chs. 4-5]
|
|- chapter02_single_epidemic_identifiability/
|   |- 01_baseline_cases.R              # Scenarios S3 (hom+NPI) and S4 (het+NPI):
|   |                                   #   fits both models to het data  [Ch. 2 main]
|   |- Fitting_stochastic_tau_leap_reducedm.R  # MLE under tau-leaping stochastic
|   |                                          #   simulation  [Ch. 2, Sec. 2.7]
|   `- HELPERS_TAU_SIM.R               # tau-leaping simulator and ODE helpers;
|                                      #   sourced by the fitting script above
|
|- chapter03_two_epidemic_inference/
|   |- 02_mle_single_epidemic.R        # single-epidemic ODE MLE; establishes
|   |                                  #   the identifiability problem  [Ch. 3]
|   `- 03_mle_two_epidemics.R          # two-epidemic joint MLE; resolves the
|                                      #   nu-c* compensation ridge  [Ch. 3]
|
|- appendixC_supplementary/            # Appendix C: supplementary for Chs. 2-3
|   |- 04_single_epidemic_correlation.R   # parameter correlation, single epidemic
|   |- 05_two_epidemics_correlation.R     # parameter correlation, two epidemics
|   `- 06_correlation_npi_comparison.R    # combined figure across NPI strengths;
|                                         #   run after six correlation runs complete
|
|- chapter05_distributional_misspecification/
|   |- single_epidemic_misspec.R       # Gamma vs Lognormal misspecification,
|   |                                  #   one epidemic  [Ch. 5]
|   `- two_epidemics_misspec.R         # Gamma vs Lognormal misspecification,
|                                      #   two epidemics  [Ch. 5]
|
|- chapter06_england_scotland/
|   |- data/
|   |   `- GB_data.csv                 # daily deaths and OxCGRT stringency,
|   |                                  #   England and Scotland,
|   |                                  #   31 Jan 2020 to 01 Nov 2021 (640 obs.)
|   |
|   |   --- Hierarchical priors (main analysis) ---
|   |- linear_npi_fit.R                # piecewise-linear NPI, hierarchical priors
|   |- logistic_npi_fit.R              # logistic NPI, hierarchical priors
|   |- stringency_npi_fit.R            # stringency-index NPI, hierarchical priors
|   |- compare_npi_specifications.R    # LOOIC comparison + diagnostics, hierarchical
|   |
|   |   --- Uniform priors (sensitivity analysis) ---
|   |- linear_npi_fit_uniform.R        # piecewise-linear NPI, uniform priors
|   |- logistic_npi_fit_uniform.R      # logistic NPI, uniform priors
|   |- stringency_npi_fit_uniform.R    # stringency-index NPI, uniform priors
|   `- compare_npi_specifications_uniform.R  # LOOIC comparison + diagnostics, uniform
|
`- outputs/                            # created at runtime; not tracked by git
    |- figures/
    |- results/
    `- stan_fits/                      # Stan fit objects (.rds) and compiled models
```

---

## Chapter mapping

**Chapter numbering** follows the revised (post-viva correction) thesis.

| Submitted | Revised | Content |
|-----------|---------|---------|
| 1, 2 | 1 | Background, review, and modelling framework |
| 3 | 2 | Single-epidemic identifiability |
| 4 | 3 | Joint two-epidemic inference |
| 5 | 4 | Discretisation methods (LA-SM) |
| 6 | 5 | Distributional misspecification |
| 7 | 6 | England and Scotland application |
| 8 | 7 | General conclusions |

No code folder is provided for Chapters 1, 4, or 7. Chapter 1 is theoretical.
Chapter 4 (discretisation) is fully implemented inside
`R/Generalfun_mispec_distribution.R`. Chapter 7 contains no new computation.

The tau-leaping study was moved from submitted Appendix C into the main text
as Chapter 2, Section 2.7. Its scripts are therefore in
`chapter02_single_epidemic_identifiability/`. Revised Appendix C holds
parameter-correlation supplementary material for Chapters 2 and 3.

---

## How to reproduce

### 1. Prerequisites

- R (>= 4.2.0)
- RStan (>= 2.21): follow the [RStan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
- A C++ toolchain compatible with RStan (Rtools on Windows; Xcode CLT on macOS)

Install all required CRAN packages:

```r
source("install_packages.R")
```

### 2. Working directory

Set your working directory to the repository root, or open the repository as
an RStudio project. All `source()` and file paths in the scripts are relative
to the repository root.

### 3. Source order for Chapters 2 and 3

Every Chapter 2 and 3 analysis script begins with:

```r
source("R/utility_functions.R")   # ODE, simulation, logit/expit
source("R/mle_functions.R")       # MLE fitters and profile likelihood
```

Chapter 5 scripts begin with:

```r
source("R/Generalfun_mispec_distribution.R")  # LA-SM discretisation
```

The tau-leaping scripts are self-contained and source only
`HELPERS_TAU_SIM.R` from the same directory.

### 4. Chapter-by-chapter

#### Chapter 2: Single-epidemic identifiability

```r
# Scenarios S3 (hom+NPI) and S4 (het+NPI) fitted to het data
source("chapter02_single_epidemic_identifiability/01_baseline_cases.R")

# Section 2.7: tau-leaping stochastic identifiability study
source("chapter02_single_epidemic_identifiability/Fitting_stochastic_tau_leap_reducedm.R")
```

`01_baseline_cases.R` generates the S3 and S4 fitting and prediction plots.
Scenarios S1 and S2 (no NPI) are analytically trivial and omitted.

#### Chapter 3: Two-epidemic inference

Run the single-epidemic script first; it writes a CSV that the two-epidemic
script reads for the interval-width comparison.

```r
source("chapter03_two_epidemic_inference/02_mle_single_epidemic.R")
source("chapter03_two_epidemic_inference/03_mle_two_epidemics.R")
```

#### Appendix C: Parameter correlation

Both correlation scripts must each be run **three times**, changing
`c_value2_spec` to `0.2`, `0.3`, and `0.4` in turn. Only after all six runs
are complete does the comparison script produce output.

```r
# Run each three times with c_value2_spec = 0.2, then 0.3, then 0.4
source("appendixC_supplementary/04_single_epidemic_correlation.R")
source("appendixC_supplementary/05_two_epidemics_correlation.R")

# Once all six runs are done
source("appendixC_supplementary/06_correlation_npi_comparison.R")
```

`06_correlation_npi_comparison.R` checks for the six result folders and
stops with a readable message naming any that are missing.

#### Chapter 5: Distributional misspecification

Both scripts have a `SIMULATE_WITH` switch near the top. Set it to `"gamma"`
or `"lognormal"` and run each script once per setting.

```r
source("chapter05_distributional_misspecification/single_epidemic_misspec.R")
source("chapter05_distributional_misspecification/two_epidemics_misspec.R")
```

#### Chapter 6: England and Scotland application

Stan models compile on first run and are cached automatically. Each fitting
script saves its fit object to `outputs/stan_fits/`. Run the three fitting
scripts before the comparison script.

**Hierarchical priors (main analysis):**

```r
source("chapter06_england_scotland/linear_npi_fit.R")
source("chapter06_england_scotland/logistic_npi_fit.R")
source("chapter06_england_scotland/stringency_npi_fit.R")
source("chapter06_england_scotland/compare_npi_specifications.R")
```

**Uniform priors (sensitivity analysis):**

```r
source("chapter06_england_scotland/linear_npi_fit_uniform.R")
source("chapter06_england_scotland/logistic_npi_fit_uniform.R")
source("chapter06_england_scotland/stringency_npi_fit_uniform.R")
source("chapter06_england_scotland/compare_npi_specifications_uniform.R")
```

The comparison scripts can load fits from `outputs/stan_fits/` so they run
in a fresh session without the fitting objects in memory.

Stan sampler settings across all Chapter 6 scripts: 4 chains, 2500 iterations,
1250 warmup, `adapt_delta = 0.98`, `max_treedepth = 15`. Each fit takes
roughly 1 to 4 hours on a multi-core machine.

---

## Fit object names

| Script | Object in R | RDS file |
|--------|-------------|----------|
| `linear_npi_fit.R` | `fit_het_lin_hier` | `outputs/stan_fits/fit_het_lin_hier.rds` |
| `logistic_npi_fit.R` | `fit_het_log_hier` | `outputs/stan_fits/fit_het_log_hier.rds` |
| `stringency_npi_fit.R` | `fit_het_str_hier` | `outputs/stan_fits/fit_het_str_hier.rds` |
| `linear_npi_fit_uniform.R` | `fit_het_lin_unif` | `outputs/stan_fits/fit_het_lin_unif.rds` |
| `logistic_npi_fit_uniform.R` | `fit_het_log_unif` | `outputs/stan_fits/fit_het_log_unif.rds` |
| `stringency_npi_fit_uniform.R` | `fit_het_str_unif` | `outputs/stan_fits/fit_het_str_unif.rds` |

The `_hier` and `_unif` objects can coexist in a single R session.

---

## Data

`chapter06_england_scotland/data/GB_data.csv` contains:

| Column | Description |
|--------|-------------|
| `Date` | Calendar date (DD/MM/YY or DD/MM/YYYY) |
| `Deaths_EN` | Daily COVID-19 deaths, England |
| `Deaths_SC` | Daily COVID-19 deaths, Scotland |
| `Stringency_Total` | OxCGRT composite stringency index (UK) |
| `Stringency_EN` | OxCGRT stringency index, England |
| `Stringency_SC` | OxCGRT stringency index, Scotland |

Coverage: 31 January 2020 to 01 November 2021 (640 daily observations).

Death counts are from the UK Health Security Agency and National Records of
Scotland. Stringency indices are from the Oxford COVID-19 Government Response
Tracker (Hale et al., 2021, *Nature Human Behaviour*).

---

## Key model parameters

| Symbol | Description | Typical value |
|--------|-------------|---------------|
| R0 | Basic reproduction number | 3.0 |
| nu (v) | Susceptibility coefficient of variation | estimated |
| c* (c_value2) | NPI contact-reduction factor | estimated |
| t0 | NPI onset time | estimated |
| delta | Incubation rate (1/mean exposed period) | 1/5.5 per day |
| gamma | Recovery rate (1/mean infectious period) | 1/4 per day |
| rho | Relative infectiousness of E class | 0.5 |
| K | Gamma discretisation groups | 30 (Chs. 2-3), 100 (Ch. 5), 20 (Ch. 6) |

---

## Citation

If you use this code or the associated thesis, please cite:

```
Mohammed, I. (2026). Simultaneous inference of susceptibility distributions
and non-pharmaceutical interventions from epidemic trajectories.
PhD thesis, University of Strathclyde, Glasgow. Student ID 202279132.
```

---

## Authors

- **Ibrahim Mohammed** (lead; University of Strathclyde)
- **Chris Robertson** (supervisor; University of Strathclyde)
- **M. Gabriela M. Gomes** (supervisor; University of Strathclyde)

---

## Licence

The code in this repository is released under the [MIT Licence](LICENSE).
The data file `GB_data.csv` is derived from publicly available sources
(UKHSA and NRS) and is provided here for reproducibility only; please
cite the original sources if you use the data independently.

---

## Status

| Component | State |
|-----------|-------|
| `R/utility_functions.R` | Complete |
| `R/mle_functions.R` | Complete |
| `R/Generalfun_mispec_distribution.R` | Complete |
| `chapter02/01_baseline_cases.R` (S3, S4) | Complete |
| `chapter02/Fitting_stochastic_tau_leap_reducedm.R` | Complete |
| `chapter02/HELPERS_TAU_SIM.R` | Complete |
| `chapter03/02_mle_single_epidemic.R` | Complete |
| `chapter03/03_mle_two_epidemics.R` | Complete |
| `appendixC/04_single_epidemic_correlation.R` | Complete |
| `appendixC/05_two_epidemics_correlation.R` | Complete |
| `appendixC/06_correlation_npi_comparison.R` | Complete |
| `chapter05/single_epidemic_misspec.R` | Complete |
| `chapter05/two_epidemics_misspec.R` | Complete |
| `chapter06/data/GB_data.csv` | Complete |
| `chapter06/linear_npi_fit.R` (hierarchical) | Complete |
| `chapter06/logistic_npi_fit.R` (hierarchical) | Complete |
| `chapter06/stringency_npi_fit.R` (hierarchical) | Complete |
| `chapter06/compare_npi_specifications.R` | Complete |
| `chapter06/linear_npi_fit_uniform.R` | Complete |
| `chapter06/logistic_npi_fit_uniform.R` | Complete |
| `chapter06/stringency_npi_fit_uniform.R` | Complete |
| `chapter06/compare_npi_specifications_uniform.R` | Complete |
