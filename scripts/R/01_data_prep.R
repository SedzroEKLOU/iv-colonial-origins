# =============================================================================
# 01_data_prep.R
# =============================================================================
# Purpose : Load ajrcomment.dta and prepare the analysis dataset.
#
# Variable inventory (from glimpse output):
#   longname  — country name
#   shortnam  — ISO-3 code
#   logmort0  — log settler mortality per 1000 per annum (instrument)
#   risk      — average protection against expropriation risk, 1985-1995
#               (0 = worst, 10 = best; endogenous institutions variable)
#   loggdp    — log GDP per capita (PPP) in 1995 (outcome)
#   latitude  — absolute latitude / 90 (geography control)
#   asia      — continent dummy
#   africa    — continent dummy
#   other     — Australia/NZ/Canada/US not counted as neo-Europe in some specs
#   neoeuro   — 1 = neo-European settler country (AUS, CAN, NZL, USA)
#   edes1975  — % European descent in 1975 (alternative control)
#   malaria   — malaria index (alternative control)
#   campaign  — Albouy: 1 = mortality data from military campaigns (not settlers)
#   slave     — Albouy: 1 = mortality partly from African slave labour records
#   source0   — 1 = AJR's own primary source data (vs. conjectured by AJR)
#   cons90    — constraints on executive in 1900 (alternative institutions var)
#   lado1995  — labour rights index 1995 (alternative institutions var)
#   mort      — settler mortality (raw, not log)
#   step      — AJR's step variable for data quality
#
# Output: data/clean/ajr_clean.rds
# =============================================================================

library(tidyverse)
library(haven)
library(here)

cat("Loading ajrcomment.dta...\n")

df <- read_dta(here("data", "raw", "ajrcomment.dta")) %>%
  as_tibble() %>%
  # Drop Zap attributes from Stata labels (keep values only)
  mutate(across(where(is.labelled), as_factor)) %>%
  mutate(across(where(is.factor), as.character))

cat(sprintf("Loaded: %d observations, %d variables\n", nrow(df), ncol(df)))

# =============================================================================
# 1.  Derived variables
# =============================================================================

df <- df %>%
  mutate(
    # Region identifier (mutually exclusive)
    region = case_when(
      africa  == 1 ~ "Africa",
      asia    == 1 ~ "Asia",
      neoeuro == 1 ~ "Neo-Europe",
      other   == 1 ~ "Other",
      TRUE         ~ "Americas/Other"
    ),

    # Albouy critique flag: data quality concerns
    albouy_concern = as.integer(campaign == 1 | slave == 1),

    # Primary source only (AJR's own data, not conjectured)
    primary_source = as.integer(source0 == 1),

    # Log of alternative institutions measures (where possible)
    logcons90 = log(cons90 + 1),

    # Normalised latitude (already 0-1 in the data as latitude/90)
    lat = latitude
  )

# =============================================================================
# 2.  Sample summary
# =============================================================================

cat("\n=== Sample summary ===\n")
cat(sprintf("Total countries      : %d\n", nrow(df)))
cat(sprintf("By region:\n"))
df %>% count(region) %>% print()
cat(sprintf("\nWith settler mortality data (source0=1): %d\n",
            sum(df$primary_source)))
cat(sprintf("Albouy concern (campaign or slave)    : %d\n",
            sum(df$albouy_concern)))
cat(sprintf("Neo-Europe countries                   : %d\n",
            sum(df$neoeuro)))

cat("\n=== Key variables ===\n")
df %>%
  select(logmort0, risk, loggdp, lat, edes1975, malaria) %>%
  summary() %>%
  print()

# =============================================================================
# 3.  Save
# =============================================================================

saveRDS(df, here("data", "clean", "ajr_clean.rds"))
cat("\nSaved: data/clean/ajr_clean.rds\n")
