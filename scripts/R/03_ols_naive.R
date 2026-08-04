# =============================================================================
# 03_ols_naive.R
# =============================================================================
# Purpose : OLS estimates — illustrate the endogeneity problem.
#
# Why OLS is biased here:
#   The structural equation is: log(GDP) = α + β·institutions + ε
#   But institutions are ENDOGENOUS: richer countries can afford better
#   property rights and rule of law. So Cov(institutions, ε) ≠ 0.
#   OLS estimate of β is biased UPWARD (positive omitted variable bias).
#
# In addition, institutions are measured with error (subjective risk scores),
#   which creates attenuation bias (towards zero). The net direction of OLS
#   bias in the AJR context is ambiguous in theory, but OLS is far larger
#   than the IV estimate in practice — suggesting reverse causality dominates.
#
# Specifications (matching AJR Table 2 / table1.do logic):
#   1. Baseline: loggdp ~ risk
#   2. + latitude
#   3. + continent dummies
#   4. + continent dummies + latitude
#   5. + % European descent (edes1975)
#   6. + malaria index
#
# Output:
#   output/tables/tab_ols.tex
# =============================================================================

library(tidyverse)
library(fixest)
library(modelsummary)
library(gt)
library(here)

cat("Loading data...\n")
df <- readRDS(here("data", "clean", "ajr_clean.rds"))

# =============================================================================
# OLS specifications
# =============================================================================

cat("Estimating OLS...\n")

ols1 <- feols(loggdp ~ risk,
              data = df, vcov = "HC1")

ols2 <- feols(loggdp ~ risk + lat,
              data = df, vcov = "HC1")

ols3 <- feols(loggdp ~ risk + asia + africa + other,
              data = df, vcov = "HC1")

ols4 <- feols(loggdp ~ risk + asia + africa + other + lat,
              data = df, vcov = "HC1")

ols5 <- feols(loggdp ~ risk + edes1975,
              data = df %>% filter(!is.na(edes1975)), vcov = "HC1")

ols6 <- feols(loggdp ~ risk + malaria,
              data = df %>% filter(!is.na(malaria)), vcov = "HC1")

ols_models <- list(
  "(1) Baseline"        = ols1,
  "(2) + Latitude"      = ols2,
  "(3) + Continents"    = ols3,
  "(4) + Both"          = ols4,
  "(5) + Euro descent"  = ols5,
  "(6) + Malaria"       = ols6
)

cat("\n=== OLS results ===\n")
modelsummary(
  ols_models,
  coef_map = c(
    "risk"     = "Institutions (risk)",
    "lat"      = "Latitude",
    "asia"     = "Asia",
    "africa"   = "Africa",
    "other"    = "Other",
    "edes1975" = "% European descent",
    "malaria"  = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title    = "OLS Estimates: Institutions → Log GDP per Capita",
  notes    = paste0(
    "Notes: OLS with HC1 heteroskedasticity-robust SEs. ",
    "Outcome: log GDP per capita (PPP, 1995). ",
    "Institutions: average protection against expropriation risk, 1985-1995 (0-10). ",
    "OLS is upward-biased due to reverse causality (rich countries invest in institutions). ",
    "IV estimates in Table 3 correct for this endogeneity."
  ),
  output = here("output", "tables", "tab_ols.tex")
)
modelsummary(
  ols_models,
  coef_map = c(
    "risk"     = "Institutions (risk)",
    "lat"      = "Latitude",
    "asia"     = "Asia",
    "africa"   = "Africa",
    "other"    = "Other",
    "edes1975" = "% European descent",
    "malaria"  = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  output   = here("output", "tables", "tab_ols.html")
)

cat("OLS table saved.\n")
cat(sprintf("\nBaseline OLS: β(risk) = %.3f (SE = %.3f)\n",
            coef(ols1)["risk"], se(ols1)["risk"]))
cat("Interpretation: 1-unit increase in institutional quality → ",
    sprintf("%.1f%% higher GDP per capita\n", 100 * coef(ols1)["risk"]))
cat("(Compare with IV estimate in script 05 — OLS should be larger if reverse causality dominates)\n")

saveRDS(ols_models, here("output", "ols_models.rds"))
