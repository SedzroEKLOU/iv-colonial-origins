# =============================================================================
# 00_master.R
# =============================================================================
# Purpose : Master script for the R replication and extension of
#           Acemoglu, D., Johnson, S. & Robinson, J. A. (2001).
#           "The Colonial Origins of Comparative Development:
#            An Empirical Investigation."
#           American Economic Review, 91(5), 1369–1401.
#
# Design  : Instrumental Variables (IV / 2SLS)
#           Instrument Z : log settler mortality (logmort0)
#           Endogenous D : expropriation risk = institutional quality (risk)
#           Outcome    Y : log GDP per capita (loggdp)
#           N           : 64 former colonies
#
# Data    : ajrcomment.dta — AJR dataset augmented with Albouy (2012) variables
#           Source: AER replication files (ICPSR 112563)
#           Place in: data/raw/ajrcomment.dta
#
# Pipeline (4 phases):
#   Phase 1 — Data & descriptive    : scripts 01–02
#   Phase 2 — Core replication      : scripts 03–05
#   Phase 3 — Extensions            : scripts 06–09
#   Phase 4 — README
#
# Key packages : fixest, ivreg, sandwich, modelsummary, gt, tidyverse
#
# References:
#   AJR (2001)   — The Colonial Origins (main paper)
#   Albouy (2012)— AER comment (critique of settler mortality data)
#   AJR (2012)   — AER reply
# =============================================================================

library(here)
setwd(here())

run_script <- function(path, label) {
  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  %s\n", label))
  cat(sprintf("%s\n", strrep("=", 70)))
  t0 <- proc.time()
  source(path)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("\n[DONE] %s (%.1f sec)\n", label, elapsed))
}

# =============================================================================
# PHASE 1 — Data preparation and descriptive analysis
# =============================================================================
run_script(here("scripts", "R", "01_data_prep.R"),
           "[01/09] Data preparation → ajr_clean.rds")
run_script(here("scripts", "R", "02_descriptive.R"),
           "[02/09] Descriptive statistics and correlation analysis")

# =============================================================================
# PHASE 2 — Core replication
# =============================================================================
run_script(here("scripts", "R", "03_ols_naive.R"),
           "[03/09] OLS (naive, biased) — endogeneity illustration")
run_script(here("scripts", "R", "04_first_stage.R"),
           "[04/09] First stage: settler mortality → institutions (7 specs)")
run_script(here("scripts", "R", "05_iv_2sls.R"),
           "[05/09] IV/2SLS main results: institutions → GDP (7 specs)")

# =============================================================================
# PHASE 3 — Extensions
# =============================================================================
run_script(here("scripts", "R", "06_weak_iv_tests.R"),
           "[06/09] Weak IV tests + Anderson-Rubin CIs")
run_script(here("scripts", "R", "07_albouy_critique.R"),
           "[07/09] Albouy (2012) critique: campaign & slave mortality")
run_script(here("scripts", "R", "08_robustness.R"),
           "[08/09] LIML, jackknife IV, leave-one-out influence")
run_script(here("scripts", "R", "09_heterogeneity.R"),
           "[09/09] Regional heterogeneity + alternative institutions")

cat(sprintf("\n%s\n", strrep("=", 70)))
cat("  REPLICATION COMPLETE\n")
cat(sprintf("%s\n\n", strrep("=", 70)))
