# susceptibility-npi-inference

**Simultaneous inference of susceptibility distributions and non-pharmaceutical
interventions from epidemic trajectories**

Ibrahim Mohammed · Chris Robertson · M. Gabriela M. Gomes  
University of Strathclyde, Glasgow

---

> **Status:** Repository under active development. Code is being cleaned and
> documented for reproducibility. All analysis scripts and shared functions are
> present; see the [Status](#status) section for what remains outstanding.

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
to first-wave COVID-19 mortality data from England and Scotland using
hierarchical Bayesian inference implemented in Stan.

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
|   |                                   #   logit/expit  [used by Chs. 2-3]
|   |- mle_functions.R                  # Poisson likelihoods, MLE fitters,
|   |                                   #   profile likelihood, generate_trajectory
|   |                                   #   [used by Chs. 2-3]
|   `- Generalfun_mispec_distribution.R # LA-SM discretisation, heterogeneous
|                                       #   ODE and likelihoods  [used by Ch. 5]
|
|- chapter02_single_epidemic_identifiability/
|   |- 01_baseline_cases.R              # Scenarios S3 (hom+NPI) and S4 (het+NPI):
|   |                                   #   fits both models to heterogeneous data
|   |                                   #   with NPIs. S1 and S2 (no NPI) are
|   |                                   #   analytically trivial and omitted.
|   |                                   #   [Ch. 2 main]
|   |- Fitting_stochastic_tau_leap_reducedm.R  # MLE under tau-leaping stochastic
|   |                                          #   simulation  [Ch. 2, Sec. 2.7]
|   `- HELPERS_TAU_SIM.R               # tau-leaping simulator and ODE helpers;
|                                      #   sourced by Fitting_stochastic_tau_leap_reducedm.R
|
|- chapter03_two_epidemic_inference/
|   |- 02_mle_single_epidemic.R        # single-epidemic ODE MLE; establishes
|   |                                  #   the identifiability problem  [Ch. 3]
|   `- 03_mle_two_epidemics.R          # two-epidemic joint MLE; resolves the
|                                      #   nu-c* compensation ridge  [Ch. 3]
|
|- appendixC_supplementary/            # Appendix C: supplementary material for Chs. 2-3
|   |- 04_single_epidemic_correlation.R  # parameter correlation across initial
|   |                                    #   conditions, single epidemic
|   `- 05_two_epidemics_correlation.R    # parameter correlation across initial
|                                        #   conditions, two epidemics
|
|- chapter05_distributional_misspecification/
|   |- single_epidemic_misspec.R       # Gamma vs Lognormal misspecification,
|   |                                  #   one epidemic  [Ch. 5]
|   `- two_epidemics_misspec.R         # Gamma vs Lognormal misspecification,
|                                      #   two epidemics  [Ch. 5]
|
|- chapter06_england_scotland/
|   |- data/
|   |   `- GB_data.csv                 # daily deaths and OxCGRT stringency index,
|   |                                  #   England and Scotland,
|   |                                  #   31 Jan 2020 to 01 Nov 2021 (640 obs.)
|   |- linear_npi_fit.R                # joint hierarchical fit, piecewise-linear NPI
|   |- logistic_npi_fit.R              # joint hierarchical fit, logistic NPI
|   |- stringency_npi_fit.R            # joint hierarchical fit, stringency-index NPI
|   `- compare_npi_specifications.R    # LOOIC comparison across the three NPI forms
|
`- outputs/                            # created at runtime; not tracked by git
    |- figures/
    |- results/
    `- stan_fits/                      # Stan fit objects (.rds)
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

**No code folder is provided for Chapters 1, 4, or 7.**
Chapter 1 is theoretical. Chapter 4 (discretisation) is implemented inside
`R/Generalfun_mispec_distribution.R` and has no standalone analysis script.
Chapter 7 contains no new computation.

**Tau-leaping study (Section 2.7).** The tau-leaping identifiability study
was moved from submitted Appendix C into the main text as Chapter 2,
Section 2.7. Its scripts (`Fitting_stochastic_tau_leap_reducedm.R` and
`HELPERS_TAU_SIM.R`) are therefore in
`chapter02_single_epidemic_identifiability/`. `HELPERS_TAU_SIM.R` is a
self-contained helper sourced by the fitting script; it carries its own ODE
and does not depend on `R/utility_functions.R`.

**Revised Appendix C** holds supplementary parameter-correlation material for
Chapters 2 and 3. Those scripts are in `appendixC_supplementary/`.

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

Set your working directory to the repository root before sourcing any script,
or open the repository as an RStudio project. All `source()` paths in the
chapter scripts are relative to the repository root.

### 3. Chapter-by-chapter reproduction

Scripts in Chapters 2 and 3 begin with:

```r
source("R/utility_functions.R")   # ODE, simulation
source("R/mle_functions.R")       # MLE, profile likelihood
```

Chapter 5 scripts begin with:

```r
source("R/Generalfun_mispec_distribution.R")  # LA-SM discretisation
```

The tau-leaping scripts are self-contained and source only `HELPERS_TAU_SIM.R`.

#### Chapter 2: Single-epidemic identifiability

```r
# Scenarios S3 (hom+NPI) and S4 (het+NPI) fitted to heterogeneous data
source("chapter02_single_epidemic_identifiability/01_baseline_cases.R")

# Section 2.7: tau-leaping stochastic identifiability study
source("chapter02_single_epidemic_identifiability/Fitting_stochastic_tau_leap_reducedm.R")
```

`01_baseline_cases.R` generates the S3 and S4 fitting and prediction plots.
S1 (hom, no NPI) and S2 (het, no NPI) are analytically trivial and are not
included as separate scripts. `Fitting_stochastic_tau_leap_reducedm.R`
reproduces the stochastic simulation study in Section 2.7; its helper file
`HELPERS_TAU_SIM.R` must be in the same directory.

#### Chapter 3: Two-epidemic inference

Run scripts in order: the single-epidemic script provides the identifiability
comparator.

```r
source("chapter03_two_epidemic_inference/02_mle_single_epidemic.R")
source("chapter03_two_epidemic_inference/03_mle_two_epidemics.R")
```

`02_mle_single_epidemic.R` demonstrates the nu-c* ridge on a single epidemic.
`03_mle_two_epidemics.R` shows how joint fitting across two epidemics with
different initial conditions resolves the ridge.

#### Appendix C: Parameter correlation supplementary material

```r
source("appendixC_supplementary/04_single_epidemic_correlation.R")
source("appendixC_supplementary/05_two_epidemics_correlation.R")
```

Both scripts loop over a grid of I0 values and are computationally intensive.
Set `n_replicates` at the top of each script to a smaller value (e.g. 20)
for a quick verification run.

#### Chapter 5: Distributional misspecification

```r
source("chapter05_distributional_misspecification/single_epidemic_misspec.R")
source("chapter05_distributional_misspecification/two_epidemics_misspec.R")
```

Both use K = 100 discretisation groups and n_replicates = 200 by default.
The random seed is fixed at 20251217.

#### Chapter 6: England and Scotland application

Stan models are compiled on first run and cached automatically.
Each fitting script saves its Stan fit object to `outputs/stan_fits/`.
Run the three fitting scripts before the comparison script.

```r
source("chapter06_england_scotland/linear_npi_fit.R")
source("chapter06_england_scotland/logistic_npi_fit.R")
source("chapter06_england_scotland/stringency_npi_fit.R")
source("chapter06_england_scotland/compare_npi_specifications.R")
```

Stan runs use `adapt_delta = 0.95` to `0.98` and `max_treedepth = 12` to `15`
depending on the NPI specification. Running all three fits at the default
settings requires several hours on a multi-core machine.
`options(mc.cores = parallel::detectCores())` is set at the top of each script.

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
| delta | Incubation rate (1/mean exposed time) | 1/5.5 per day |
| gamma | Recovery rate (1/mean infectious time) | 1/4 per day |
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
| `R/mle_functions.R` (merged, het only) | Complete |
| `R/Generalfun_mispec_distribution.R` | Complete |
| `chapter02`: `01_baseline_cases.R` (S3, S4) | Complete |
| `chapter02`: `Fitting_stochastic_tau_leap_reducedm.R` + `HELPERS_TAU_SIM.R` | Complete |
| `chapter03`: `02_mle_single_epidemic.R` | Complete |
| `chapter03`: `03_mle_two_epidemics.R` | Complete |
| `appendixC`: `04_single_epidemic_correlation.R` | Complete |
| `appendixC`: `05_two_epidemics_correlation.R` | Complete |
| `chapter05`: misspecification scripts (het only) | Complete |
| `chapter06`: Stan fitting scripts (het only) | In progress |
| `chapter06`: `compare_npi_specifications.R` | Complete |
| `chapter06/data/GB_data.csv` | Complete |
| End-to-end test on clean R installation | Pending |

