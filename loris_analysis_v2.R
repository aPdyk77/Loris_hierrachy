# ==============================================================================
# BORNEAN SLOW LORIS — REBUILT HIERARCHICAL ANALYSIS (v2)
# ==============================================================================
# FIELD-DATA SOURCE
#   Only the five transect sheets are read as observation data:
#     TR-01 PA, TR-02 PA, TR-03 PA, TR-04 PA, TR-05 PA
#   Other workbook sheets are NOT ingested as survey observations.
#
# DESIGN
#   5 transects x 5 segments x 2 visits x 2 protocols = 100 opportunities
#   Current workbook audit expectation:
#     73 completed protocol-level observations
#     27 missing/unconducted observations
#      7 positive protocol-level detection records
#      5 segments with >=1 detection (naive segment use = 20%)
#
# ECOLOGICAL STATE
#   psi_i = probability that 200-m segment i belongs to the set of locations
#   used by slow lorises during the July–September study period.
#   This is seasonal segment use, NOT literal continuous presence.
#
# OBSERVATION PROCESS
#   p_r = probability of a positive protocol-level observation conditional on
#   seasonal segment use. Because VES and Point Count are sequential within the
#   same nightly visit, p should not be interpreted as pure detectability alone;
#   it may also absorb within-visit local availability.
#
# PRIMARY MODEL
#   State:       distance from dry riverbed + transect random intercept
#   Observation: protocol (VES vs Point Count) + visit (September vs July)
#
# SENSITIVITY MODELS
#   1) Non-spatial state model
#   2) Stationary AR(1) residual process along ordered 200-m segments
#
# IMPORTANT
#   Canopy is intentionally excluded from the quantitative model because the
#   existing logbook does not contain a systematic, comparable segment-level
#   structural canopy-closure measurement across all 25 segments.
# ============================================================================== 

# ---- 0. PACKAGES --------------------------------------------------------------
required <- c("readxl", "dplyr", "purrr", "stringr", "tidyr", "tibble", "rstan", "readr", "ggplot2")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install required packages first: ", paste(missing_pkgs, collapse = ", "))
}

library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)
library(rstan)

options(mc.cores = min(4L, parallel::detectCores()))
rstan_options(auto_write = TRUE)

# ---- 1. FILE PATHS -------------------------------------------------------------
# Put this script, the Excel workbook, and the three .stan files in the same
# working directory, or edit the paths below.
excel_file <- "Loris_data collection Sheet.xlsx"
stan_nonspatial_file <- "loris_nonspatial_v2.stan"
stan_re_file         <- "loris_transect_RE_v2.stan"
stan_ar1_file        <- "loris_AR1_v2.stan"

stopifnot(file.exists(excel_file))
stopifnot(file.exists(stan_nonspatial_file))
stopifnot(file.exists(stan_re_file))
stopifnot(file.exists(stan_ar1_file))

# ---- 2. TRANSECT METADATA USED ONLY FOR DISTANCE CALIBRATION ------------------
# The observation rows are read ONLY from TR-01 PA ... TR-05 PA.
# These starting distances are fixed metadata used to calibrate the ecological
# distance gradient. They correspond to the recorded distance from each transect
# origin to the dry-riverbed reference.
transect_meta <- tribble(
  ~transect, ~start_dist_river_m,
  "TR-01",   5.2,
  "TR-02",   4.3,
  "TR-03",   7.1,
  "TR-04",   1.7,
  "TR-05",   2.4
)

stopifnot(nrow(transect_meta) == 5L)
stopifnot(all(is.finite(transect_meta$start_dist_river_m)))

# ---- 3. READ THE FIVE TRANSECT SHEETS -----------------------------------------
sheet_names <- sprintf("TR-%02d PA", 1:5)

read_transect_sheet <- function(sheet_name) {
  # In these sheets, rows 1–7 are metadata/blank and row 8 is the table header.
  x <- read_excel(
    excel_file,
    sheet = sheet_name,
    skip = 7,
    na = c("", "NA")
  ) %>%
    select(1:10)

  names(x) <- c(
    "transect_raw", "rep", "date", "station", "time_block",
    "peat", "cloud", "loris", "stopwatch", "notes"
  )

  x %>%
    mutate(
      source_sheet = sheet_name,
      transect = str_extract(sheet_name, "TR-\\d{2}"),
      station = str_trim(as.character(station)),
      rep = suppressWarnings(as.integer(rep)),
      peat = suppressWarnings(as.numeric(peat)),
      cloud = suppressWarnings(as.numeric(cloud)),
      loris = suppressWarnings(as.integer(loris))
    ) %>%
    filter(
      rep %in% c(1L, 2L),
      str_detect(station, "^(VES|ST)-0?[1-5]$")
    ) %>%
    mutate(
      segment = as.integer(str_extract(station, "\\d+$")),
      protocol = if_else(str_starts(station, "VES"), "VES", "PointCount"),
      protocol_ves = if_else(protocol == "VES", 1L, 0L),
      visit = if_else(rep == 1L, "July", "September"),
      visit_sept = if_else(rep == 2L, 1L, 0L),
      completed = !is.na(loris)
    )
}

survey_all <- map_dfr(sheet_names, read_transect_sheet) %>%
  arrange(transect, rep, segment, desc(protocol_ves))

# ---- 4. RAW-DATA STRUCTURE AUDIT ----------------------------------------------
# Each transect sheet should contribute 20 intended opportunities:
# 2 visits x 5 segments x 2 protocols.
per_sheet_n <- survey_all %>% count(source_sheet, name = "n")
print(per_sheet_n)
stopifnot(nrow(per_sheet_n) == 5L)
stopifnot(all(per_sheet_n$n == 20L))

# Every transect x visit x segment should contain exactly VES + Point Count rows,
# even when the actual survey result is missing.
cell_audit <- survey_all %>%
  count(transect, rep, segment, name = "n_protocol_rows")
stopifnot(nrow(cell_audit) == 50L)
stopifnot(all(cell_audit$n_protocol_rows == 2L))

stopifnot(nrow(survey_all) == 100L)
stopifnot(all(survey_all$segment %in% 1:5))
stopifnot(all(is.na(survey_all$loris) | survey_all$loris %in% c(0L, 1L)))

n_completed <- sum(!is.na(survey_all$loris))
n_missing   <- sum(is.na(survey_all$loris))
n_positive  <- sum(survey_all$loris == 1L, na.rm = TRUE)

cat("\nRAW DATA AUDIT\n")
cat("Total intended opportunities :", nrow(survey_all), "\n")
cat("Completed observations       :", n_completed, "\n")
cat("Missing/unconducted          :", n_missing, "\n")
cat("Positive protocol records    :", n_positive, "\n\n")

# These deliberately fail if the workbook changes.
stopifnot(n_completed == 73L)
stopifnot(n_missing == 27L)
stopifnot(n_positive == 7L)

cat("Completed/missing by visit and protocol:\n")
print(
  survey_all %>%
    mutate(status = if_else(completed, "completed", "missing")) %>%
    count(visit, protocol, status, name = "n") %>%
    arrange(visit, protocol, status)
)

cat("\nPositive protocol-level records:\n")
print(
  survey_all %>%
    filter(loris == 1L) %>%
    select(transect, visit, station, peat, cloud, notes)
)

# ---- 5. DEFINE THE 25 ECOLOGICAL SEGMENTS -------------------------------------
# Transects move away from the dry riverbed. Because psi is defined for the full
# 200-m segment, use the midpoint as its ecological distance support.
# For origin distance d0, segment midpoints are d0 + 100, 300, 500, 700, 900 m.
site_tbl <- expand_grid(
  transect = sprintf("TR-%02d", 1:5),
  segment = 1:5
) %>%
  left_join(transect_meta, by = "transect") %>%
  arrange(transect, segment) %>%
  mutate(
    site_id = row_number(),
    transect_id = match(transect, sprintf("TR-%02d", 1:5)),
    distance_mid_m = start_dist_river_m + 100 + 200 * (segment - 1)
  )

# Standardize with an explicit mean and SD so the transformation is reproducible.
dist_mean <- mean(site_tbl$distance_mid_m)
dist_sd   <- sd(site_tbl$distance_mid_m)
site_tbl <- site_tbl %>%
  mutate(dist_z = (distance_mid_m - dist_mean) / dist_sd)

stopifnot(nrow(site_tbl) == 25L)
stopifnot(!anyNA(site_tbl$distance_mid_m))
stopifnot(isTRUE(all.equal(mean(site_tbl$dist_z), 0, tolerance = 1e-12)))

cat("\nDistance standardization:\n")
cat("Mean midpoint distance =", round(dist_mean, 3), "m\n")
cat("SD midpoint distance   =", round(dist_sd, 3), "m\n\n")

# ---- 6. JOIN OBSERVATIONS TO SEGMENTS -----------------------------------------
survey_all <- survey_all %>%
  left_join(
    site_tbl %>% select(transect, segment, site_id, transect_id,
                        distance_mid_m, dist_z),
    by = c("transect", "segment")
  )

stopifnot(!anyNA(survey_all$site_id))

# Missing/unconducted rows remain in the 100-opportunity audit table, but are
# excluded from the observation likelihood. They are NEVER recoded as zeros.
events <- survey_all %>%
  filter(completed) %>%
  arrange(site_id, rep, segment, desc(protocol_ves))

stopifnot(nrow(events) == 73L)
stopifnot(sum(events$loris == 1L) == 7L)

n_detected_segments <- events %>%
  group_by(site_id) %>%
  summarise(ever_detected = any(loris == 1L), .groups = "drop") %>%
  summarise(n = sum(ever_detected)) %>%
  pull(n)

stopifnot(n_detected_segments == 5L)

cat("Detected segments          :", n_detected_segments, "of 25\n")
cat("Naive segment use          :", sprintf("%.1f%%", 100 * n_detected_segments / 25), "\n\n")

# Export transparent audit tables.
write.csv(survey_all, "loris_100_opportunities_audit.csv", row.names = FALSE)
write.csv(events,     "loris_73_completed_events.csv",   row.names = FALSE)
write.csv(site_tbl,   "loris_25_segment_covariates.csv", row.names = FALSE)

# ---- 7. DESCRIPTIVE OBSERVATION CHECKS ----------------------------------------
# Descriptive only; do not treat these as adjusted detection probabilities.
cat("Raw protocol detection proportions (descriptive):\n")
print(
  events %>%
    group_by(protocol) %>%
    summarise(n = n(), positives = sum(loris), raw_positive_rate = mean(loris),
              .groups = "drop")
)

cat("\nRaw visit detection proportions (descriptive):\n")
print(
  events %>%
    group_by(visit) %>%
    summarise(n = n(), positives = sum(loris), raw_positive_rate = mean(loris),
              .groups = "drop")
)

cat("\nPeat dryness by visit (completed observations):\n")
print(events %>% count(visit, peat, name = "n"))

cat("\nCloud index by visit (completed observations):\n")
print(events %>% count(visit, cloud, name = "n"))

# ---- 8. STAN DATA --------------------------------------------------------------
# Stronger regularization than the legacy Normal(0, 2.5) is used because there
# are only seven positive protocol-level records. All slope priors are on the
# log-odds scale; distance is standardized.
prior_intercept_sd <- 1.5
prior_slope_sd     <- 1.0
prior_re_sd        <- 0.75
prior_rho_sd       <- 0.50

common <- list(
  N = nrow(site_tbl),
  R = nrow(events),
  site_id = as.integer(events$site_id),
  y = as.integer(events$loris),
  dist_z = as.vector(site_tbl$dist_z),
  protocol_ves = as.numeric(events$protocol_ves),
  visit_sept = as.numeric(events$visit_sept),
  prior_intercept_sd = prior_intercept_sd,
  prior_slope_sd = prior_slope_sd
)

stan_data_nonspatial <- common

stan_data_re <- c(
  common,
  list(
    K = 5L,
    transect_id = as.integer(site_tbl$transect_id),
    prior_re_sd = prior_re_sd
  )
)

stan_data_ar1 <- c(
  common,
  list(
    K = 5L,
    S = 5L,
    transect_id = as.integer(site_tbl$transect_id),
    segment_id = as.integer(site_tbl$segment),
    prior_re_sd = prior_re_sd,
    prior_rho_sd = prior_rho_sd
  )
)

# ---- 9. FIT THE THREE SPATIAL FORMULATIONS ------------------------------------
# Primary candidate = transect random intercept.
# Non-spatial and AR1 are sensitivity formulations, not a winner-take-all contest.
fit_nonspatial <- stan(
  file = stan_nonspatial_file,
  data = stan_data_nonspatial,
  chains = 4, iter = 4000, warmup = 2000,
  seed = 20260921,
  control = list(adapt_delta = 0.98, max_treedepth = 12),
  refresh = 250
)

fit_re <- stan(
  file = stan_re_file,
  data = stan_data_re,
  chains = 4, iter = 4000, warmup = 2000,
  seed = 20260922,
  control = list(adapt_delta = 0.99, max_treedepth = 13),
  refresh = 250
)

fit_ar1 <- stan(
  file = stan_ar1_file,
  data = stan_data_ar1,
  chains = 4, iter = 4000, warmup = 2000,
  seed = 20260923,
  control = list(adapt_delta = 0.99, max_treedepth = 13),
  refresh = 250
)

saveRDS(fit_nonspatial, "fit_loris_nonspatial_v2.rds")
saveRDS(fit_re,         "fit_loris_transect_RE_v2.rds")
saveRDS(fit_ar1,        "fit_loris_AR1_v2.rds")

# ---- 10. PARAMETER SUMMARIES ---------------------------------------------------
base_pars <- c("beta0", "beta_dist", "gamma0", "gamma_ves", "gamma_sept")

summarise_fit <- function(fit, model_name, extra_pars = character()) {
  pars <- c(base_pars, extra_pars)
  s <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  tibble(
    model = model_name,
    parameter = rownames(s),
    mean = s[, "mean"],
    sd = s[, "sd"],
    q2.5 = s[, "2.5%"],
    median = s[, "50%"],
    q97.5 = s[, "97.5%"],
    n_eff = s[, "n_eff"],
    Rhat = s[, "Rhat"]
  )
}

sum_nonspatial <- summarise_fit(fit_nonspatial, "Non-spatial")
sum_re         <- summarise_fit(fit_re, "Transect RE", "sigma_transect")
sum_ar1        <- summarise_fit(fit_ar1, "AR1", c("sigma_omega", "rho"))

parameter_summary <- bind_rows(sum_nonspatial, sum_re, sum_ar1)
print(parameter_summary, n = Inf)
write.csv(parameter_summary, "loris_model_parameter_summary_v2.csv", row.names = FALSE)

# ---- 11. DIRECTIONAL POSTERIOR PROBABILITIES ----------------------------------
# Report posterior mass on both sides of zero rather than forcing a binary
# significant/non-significant interpretation.
directional_probs <- function(fit, model_name) {
  d <- as.matrix(fit, pars = c("beta_dist", "gamma_ves", "gamma_sept"))
  tibble(
    model = model_name,
    Pr_beta_distance_gt_0 = mean(d[, "beta_dist"] > 0),
    Pr_beta_distance_lt_0 = mean(d[, "beta_dist"] < 0),
    Pr_gamma_VES_gt_0 = mean(d[, "gamma_ves"] > 0),
    Pr_gamma_VES_lt_0 = mean(d[, "gamma_ves"] < 0),
    Pr_gamma_September_gt_0 = mean(d[, "gamma_sept"] > 0),
    Pr_gamma_September_lt_0 = mean(d[, "gamma_sept"] < 0)
  )
}

posterior_direction <- bind_rows(
  directional_probs(fit_nonspatial, "Non-spatial"),
  directional_probs(fit_re,         "Transect RE"),
  directional_probs(fit_ar1,        "AR1")
)

print(posterior_direction)
write.csv(posterior_direction,
          "loris_directional_posterior_probabilities_v2.csv",
          row.names = FALSE)

# ---- 12. POSTERIOR PREDICTIVE CHECKS ------------------------------------------
ppc_summary <- function(fit, model_name) {
  d <- as.matrix(fit, pars = c("n_positive_rep", "n_detected_segments_rep"))
  tibble(
    model = model_name,
    obs_positive_records = 7L,
    pred_positive_median = median(d[, "n_positive_rep"]),
    pred_positive_l95 = unname(quantile(d[, "n_positive_rep"], 0.025)),
    pred_positive_u95 = unname(quantile(d[, "n_positive_rep"], 0.975)),
    obs_detected_segments = 5L,
    pred_segments_median = median(d[, "n_detected_segments_rep"]),
    pred_segments_l95 = unname(quantile(d[, "n_detected_segments_rep"], 0.025)),
    pred_segments_u95 = unname(quantile(d[, "n_detected_segments_rep"], 0.975))
  )
}

ppc <- bind_rows(
  ppc_summary(fit_nonspatial, "Non-spatial"),
  ppc_summary(fit_re,         "Transect RE"),
  ppc_summary(fit_ar1,        "AR1")
)
print(ppc)
write.csv(ppc, "loris_posterior_predictive_summary_v2.csv", row.names = FALSE)

# ---- 13. SEGMENT-LEVEL POSTERIOR USE ------------------------------------------
site_use_summary <- function(fit, model_name) {
  x <- rstan::extract(fit, pars = "post_prob_use")$post_prob_use
  tibble(
    site_id = 1:ncol(x),
    model = model_name,
    post_use_mean = colMeans(x),
    post_use_median = apply(x, 2, median),
    post_use_l95 = apply(x, 2, quantile, 0.025),
    post_use_u95 = apply(x, 2, quantile, 0.975)
  ) %>%
    left_join(
      site_tbl %>% select(site_id, transect, segment, distance_mid_m),
      by = "site_id"
    )
}

site_use <- bind_rows(
  site_use_summary(fit_nonspatial, "Non-spatial"),
  site_use_summary(fit_re,         "Transect RE"),
  site_use_summary(fit_ar1,        "AR1")
)
write.csv(site_use, "loris_segment_posterior_use_v2.csv", row.names = FALSE)

# ---- 14. CONVERGENCE / SAMPLER CHECK ------------------------------------------
check_fit <- function(fit, label, max_treedepth) {
  s <- summary(fit)$summary
  sampler <- get_sampler_params(fit, inc_warmup = FALSE)

  div_all <- sum(vapply(
    sampler,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))

  treedepth_hits <- sum(vapply(
    sampler,
    function(x) sum(x[, "treedepth__"] >= max_treedepth),
    numeric(1)
  ))

  cat("\n", label, "\n", sep = "")
  cat("Max Rhat                :", max(s[, "Rhat"], na.rm = TRUE), "\n")
  cat("Min n_eff               :", min(s[, "n_eff"], na.rm = TRUE), "\n")
  cat("Total divergences       :", div_all, "\n")
  cat("Max-treedepth hits      :", treedepth_hits, "\n")
}

check_fit(fit_nonspatial, "NON-SPATIAL", 12L)
check_fit(fit_re,         "TRANSECT RANDOM INTERCEPT — PRIMARY", 13L)
check_fit(fit_ar1,        "AR1 — SPATIAL SENSITIVITY", 13L)

cat("\nIMPORTANT INTERPRETATION NOTES\n")
cat("1. Good Rhat/ESS indicates computational convergence, not ecological identifiability.\n")
cat("2. Compare beta_dist/gamma_ves/gamma_sept across all three state structures.\n")
cat("3. Treat gamma_ves as an overall protocol contrast, not a causal movement/noise effect.\n")
cat("4. Treat gamma_sept as a visit-level observation/availability contrast, not a causal\n")
cat("   moonlight, peat-dryness, or month effect.\n")
cat("5. Distance is also ordered along nightly transect progression, so beta_dist remains\n")
cat("   an association rather than a causal river-distance effect.\n")
cat("6. Peat dryness is strongly aligned with visit in these data. Do not include peat and\n")
cat("   visit simultaneously in the primary sparse-data model.\n")

library(ggplot2)
library(dplyr)
library(readr)

param <- read_csv("loris_model_parameter_summary_v2.csv", show_col_types = FALSE)

# data koefisien utama untuk Figure 3
fig_main_df <- param %>%
  filter(parameter %in% c("beta_dist", "gamma_ves", "gamma_sept")) %>%
  mutate(
    model = factor(model, levels = c("Non-spatial", "Transect RE", "AR1")),
    parameter = factor(
      parameter,
      levels = c("beta_dist", "gamma_ves", "gamma_sept"),
      labels = c(
        "Distance from riverbed",
        "VES vs Point Count",
        "September vs July"
      )
    )
  )

# mean antar-formulasi model pada tiap parameter
fig_main_mean <- fig_main_df %>%
  group_by(parameter) %>%
  summarise(overall_mean = mean(mean), .groups = "drop")

fig_main <- ggplot(fig_main_df, aes(x = model, y = mean)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey35"
  ) +
  geom_hline(
    data = fig_main_mean,
    aes(yintercept = overall_mean),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "#d62728"
  ) +
  geom_errorbar(
    aes(ymin = q2.5, ymax = q97.5),
    width = 0.08,
    linewidth = 0.8,
    colour = "#1f77b4"
  ) +
  geom_point(
    size = 2.6,
    colour = "#1f77b4"
  ) +
  facet_wrap(~ parameter, nrow = 1, scales = "free_y") +
  labs(
    x = "Model formulation",
    y = "Posterior estimate (log-odds)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(
  "Figure3_parallel_coefficients_journal_ready.png",
  fig_main,
  width = 10.2,
  height = 3.7,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure3_parallel_coefficients_journal_ready.pdf",
  fig_main,
  width = 10.2,
  height = 3.7,
  bg = "white"
)

fig_main

dirprob <- read_csv("loris_directional_posterior_probabilities_v2.csv")

dir_long <- bind_rows(
  dirprob %>%
    transmute(model, parameter = "Distance from riverbed",
              prob = Pr_beta_distance_gt_0),
  dirprob %>%
    transmute(model, parameter = "VES vs Point Count",
              prob = Pr_gamma_VES_gt_0),
  dirprob %>%
    transmute(model, parameter = "September vs July",
              prob = Pr_gamma_September_lt_0)
) %>%
  mutate(
    model = factor(model, levels = c("Non-spatial", "Transect RE", "AR1")),
    parameter = factor(
      parameter,
      levels = c("Distance from riverbed",
                 "VES vs Point Count",
                 "September vs July")
    )
  )

fig_prob <- ggplot(dir_long, aes(x = model, y = prob)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ parameter, nrow = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Model formulation",
    y = "Directional posterior probability",
    title = "Directional posterior probabilities across model formulations"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("FigureS1_directional_probabilities.png", fig_prob, width = 10.5, height = 4.2, dpi = 300)
fig_prob

ppc <- read_csv("loris_posterior_predictive_summary_v2.csv")

ppc_long <- bind_rows(
  ppc %>%
    transmute(
      model,
      metric = "Positive protocol-level records",
      obs = obs_positive_records,
      med = pred_positive_median,
      l95 = pred_positive_l95,
      u95 = pred_positive_u95
    ),
  ppc %>%
    transmute(
      model,
      metric = "Detected segments",
      obs = obs_detected_segments,
      med = pred_segments_median,
      l95 = pred_segments_l95,
      u95 = pred_segments_u95
    )
) %>%
  mutate(
    model = factor(model, levels = c("Non-spatial", "Transect RE", "AR1")),
    metric = factor(metric, levels = c("Positive protocol-level records",
                                       "Detected segments"))
  )

fig_ppc <- ggplot(ppc_long, aes(x = model, y = med)) +
  geom_errorbar(aes(ymin = l95, ymax = u95), width = 0.12, linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_point(aes(y = obs), shape = 4, stroke = 1.2, size = 3) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  labs(
    x = "Model formulation",
    y = "Posterior predictive summary",
    title = "Posterior predictive checks across model formulations"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("FigureS2_posterior_predictive_checks.png", fig_ppc, width = 9.5, height = 4.2, dpi = 300)
fig_ppc

fig_var <- param %>%
  filter(parameter %in% c("sigma_transect", "sigma_omega", "rho")) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("sigma_transect", "sigma_omega", "rho"),
      labels = c("Transect SD", "Spatial SD", "AR1 autocorrelation")
    )
  ) %>%
  ggplot(aes(x = parameter, y = median)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.12, linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  labs(
    x = NULL,
    y = "Posterior estimate",
    title = "Variance and autocorrelation parameters"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("FigureS3_variance_autocorrelation.png", fig_var, width = 10.5, height = 4.2, dpi = 300)
fig_var

library(rstan)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# Helper: summarize posterior quantity
# ------------------------------------------------------------

summ_prob <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.5, 0.975))
  
  sprintf(
    "%.2f [%.2f, %.2f]",
    q[2], q[1], q[3]
  )
}
capture.output(sessionInfo(), file = "sessionInfo.txt")
cat("\nSession information written to sessionInfo.txt\n")

