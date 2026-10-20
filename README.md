# Bornean slow loris hierarchical analysis

This repository contains the data source and analysis workflow supporting the manuscript:

**Hierarchical Modeling of Sparse Bornean Slow Loris Detections During Severe 2026 Dry-Season Conditions: Robust Inference Under Constrained Field Access**

The study concerns *Nycticebus borneanus* surveys in the Punggualas research zone, Sebangau National Park, Central Kalimantan, Indonesia, during July and September 2026.

## Scope of inference

The ecological state is **seasonal use of a 200-m segment**, not permanent occupancy. The observation model represents the probability of a **positive protocol-level observation conditional on seasonal segment use** and may combine local availability with conditional detectability.

The seven positive protocol-level records must **not** be interpreted as seven independent animals because individuals were not marked or otherwise individually identified.

## Repository files

- `Loris_data collection Sheet.xlsx` — original field workbook used as the observation-data source.
- `loris_analysis_v2.R` — complete analysis workflow, including data audit, model fitting, posterior summaries, directional posterior probabilities, posterior predictive checks, and figures.
- `loris_nonspatial_v2.stan` — non-spatial sensitivity model.
- `loris_transect_RE_v2.stan` — primary transect random-intercept model.
- `loris_AR1_v2.stan` — stationary AR1 spatial sensitivity model.

Running the R workflow also exports transparent intermediate audit/data tables and model summaries.

## Expected field-data audit

Before fitting any model, the R script checks the workbook structure and stops if these expected totals are not reproduced:

- 5 transects
- 5 contiguous 200-m segments per transect
- 2 primary visits (July and September)
- 2 sequential protocols per segment/visit (moving VES, then stationary Point Count)
- 100 planned protocol-level opportunities
- 73 completed protocol-level observations
- 27 unavailable or terminated opportunities
- 7 positive protocol-level records
- 5 segments with at least one positive record

Missing or terminated surveys remain missing and are never recoded as zero.

## Model formulations

All three formulations use the same observation model:

`logit(p_r) = gamma0 + gamma_ves * VES_r + gamma_sept * Sept_r`

where VES = 1 for moving VES and 0 for Point Count, and Sept = 1 for September and 0 for July.

### 1. Non-spatial sensitivity model

`logit(psi_i) = beta0 + beta_dist * Dist_i`

### 2. Primary transect random-intercept model

`logit(psi_i) = beta0 + beta_dist * Dist_i + b_k[i]`

with `b_k ~ Normal(0, sigma_transect)`.

### 3. Stationary AR1 spatial sensitivity model

`logit(psi_k,s) = beta0 + beta_dist * Dist_k,s + omega_k,s`

with a stationary AR1 residual process within each ordered transect:

- `omega_k,1 ~ Normal(0, sigma_omega)`
- `omega_k,s ~ Normal(rho * omega_k,s-1, sigma_omega * sqrt(1-rho^2))`

The AR1 model is a **sensitivity analysis**, not a claim that the spatial correlation scale is precisely identifiable from these sparse data.

## Priors

- Intercepts: Normal(0, 1.5)
- Shared slopes: Normal(0, 1.0)
- Random-effect standard deviations: half-Normal(0, 0.75)
- AR1 correlation: Normal(0, 0.5), constrained to [-0.95, 0.95]

Distance from the dry-riverbed reference is assigned to the midpoint of each 200-m segment and standardized before model fitting.

## Software requirements

The workflow uses R and RStan. Required R packages are:

- readxl
- dplyr
- purrr
- stringr
- tidyr
- tibble
- rstan
- readr
- ggplot2

RStan requires a working C++ toolchain compatible with the local R installation.

## Reproducing the analysis

1. Clone or download this repository.
2. Keep the Excel workbook, R script, and all three Stan files in the same working directory.
3. Open R in that directory.
4. Install any missing packages listed above.
5. Run:

```r
source("loris_analysis_v2.R")
```

The script uses fixed random seeds for the three model formulations.

### MCMC settings

- 4 chains
- 4,000 iterations per chain
- 2,000 warm-up iterations per chain
- Non-spatial: `adapt_delta = 0.98`, `max_treedepth = 12`
- Transect random-intercept: `adapt_delta = 0.99`, `max_treedepth = 13`
- AR1: `adapt_delta = 0.99`, `max_treedepth = 13`

## Main generated outputs

The workflow writes, among others:

- `loris_100_opportunities_audit.csv`
- `loris_73_completed_events.csv`
- `loris_25_segment_covariates.csv`
- `loris_model_parameter_summary_v2.csv`
- `loris_directional_posterior_probabilities_v2.csv`
- `loris_posterior_predictive_summary_v2.csv`
- `loris_segment_posterior_use_v2.csv`
- `Figure3_parallel_coefficients.png`
- `FigureS1_directional_probabilities.png`
- `FigureS2_posterior_predictive_checks.png`
- `FigureS3_variance_autocorrelation.png`
- `sessionInfo.txt`

The fitted model objects are saved as RDS files and are not required for reproducing the analysis because they can be regenerated from the public data and model specifications.

## Model-checking targets reported in the manuscript

The posterior predictive checks compare replicated data with:

- observed positive protocol-level records = 7
- observed segments with at least one positive record = 5

The analysis emphasizes that satisfactory R-hat and effective sample size diagnostics demonstrate computational sampling stability, not ecological identifiability.

## Missing September observations

TR-01 and TR-02 could not be revisited in September because river access failed. These observation opportunities are treated as unavailable, not as nondetections. Therefore, the September coefficient is interpreted as a visit-level observation/availability contrast for the available survey opportunities and not as a causal drought, peat-dryness, lunar, or whole-study-area temporal effect.

## Contact

Adventus Panda  
Department of Biology, Faculty of Mathematics and Natural Sciences, Universitas Palangka Raya  
Email: apanda@mipa.upr.ac.id
