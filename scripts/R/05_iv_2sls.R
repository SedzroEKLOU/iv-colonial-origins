# =============================================================================
# 05_iv_2sls.R
# =============================================================================
# Purpose : Main IV/2SLS results — institutions → log GDP per capita.
#           Replicates the core of AJR (2001) Tables 4 and 5.
#
# 2SLS estimator:
#   First stage : risk_hat = α + γ·logmort0 + controls
#   Second stage: loggdp   = α + β·risk_hat + controls
#
#   β_IV = β_RF / γ_FS = (Cov(logmort0, loggdp) / Var(logmort0))
#                       / (Cov(logmort0, risk)   / Var(logmort0))
#
# Interpretation (LATE / Local Average Treatment Effect):
#   β_IV estimates the causal effect of institutions on GDP for the
#   "compliers" — countries whose institutional quality changed because
#   of differences in settler mortality. With a single continuous instrument,
#   LATE ≈ ATE under linearity.
#
# 7 specifications (same as first stage):
#   1. No controls
#   2. + Latitude
#   3. No Neo-Europes
#   4. + Continent dummies
#   5. + Continent dummies + latitude
#   6. + % European descent
#   7. + Malaria index
#
# Outputs:
#   output/tables/tab_iv_main.tex
#   output/figures/fig_iv_coefs.pdf
#   output/figures/fig_ols_vs_iv.pdf
# =============================================================================

library(tidyverse)
library(fixest)
library(ivreg)          # for Sargan-Hansen, Hausman test
library(modelsummary)
library(ggthemes)
library(patchwork)
library(here)

cat("Loading data...\n")
df <- readRDS(here("data", "clean", "ajr_clean.rds"))

theme_set(
  theme_clean() +
    theme(
      plot.background = element_blank(),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10)
    )
)

# =============================================================================
# 1.  IV/2SLS: 7 specifications (fixest)
# =============================================================================

cat("Estimating IV/2SLS (7 specifications)...\n")

iv1 <- feols(loggdp ~ 1 | risk ~ logmort0,
             data = df, vcov = "HC1")

iv2 <- feols(loggdp ~ lat | risk ~ logmort0,
             data = df, vcov = "HC1")

iv3 <- feols(loggdp ~ 1 | risk ~ logmort0,
             data = df %>% filter(neoeuro == 0), vcov = "HC1")

iv4 <- feols(loggdp ~ asia + africa + other | risk ~ logmort0,
             data = df, vcov = "HC1")

iv5 <- feols(loggdp ~ asia + africa + other + lat | risk ~ logmort0,
             data = df, vcov = "HC1")

iv6 <- feols(loggdp ~ edes1975 | risk ~ logmort0,
             data = df %>% filter(!is.na(edes1975)), vcov = "HC1")

iv7 <- feols(loggdp ~ malaria | risk ~ logmort0,
             data = df %>% filter(!is.na(malaria)), vcov = "HC1")

iv_models <- list(
  "(1) No controls"    = iv1,
  "(2) Latitude"       = iv2,
  "(3) No neo-Europes" = iv3,
  "(4) Continents"     = iv4,
  "(5) Cont. + Lat."   = iv5,
  "(6) Euro descent"   = iv6,
  "(7) Malaria"        = iv7
)

cat("\n=== IV/2SLS results ===\n")
lapply(seq_along(iv_models), function(i) {
  mod <- iv_models[[i]]
  coefs <- coef(mod)
  ses   <- se(mod)
  nm    <- names(coefs)[str_detect(names(coefs), "fit_risk|risk")]
  cat(sprintf("Spec %d: β(institutions) = %.3f (SE = %.3f)\n",
              i, coefs[nm[1]], ses[nm[1]]))
})

# =============================================================================
# 2.  Export main results table
# =============================================================================

modelsummary(
  iv_models,
  coef_map = c(
    "fit_risk"  = "Institutions (risk) [IV]",
    "lat"       = "Latitude",
    "asia"      = "Asia",
    "africa"    = "Africa",
    "other"     = "Other",
    "edes1975"  = "% European descent",
    "malaria"   = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title    = "IV/2SLS: Effect of Institutions on Log GDP per Capita (Tables 4–5)",
  notes    = paste0(
    "Notes: 2SLS. Endogenous: expropriation risk (institutions). ",
    "Instrument: log settler mortality. HC1 robust SEs. ",
    "Spec (3) excludes neo-European countries. ",
    "IV estimates correct for endogeneity of institutional quality. ",
    "Interpretation: 1-unit increase in expropriation risk → β × 100% higher GDP per capita."
  ),
  output = here("output", "tables", "tab_iv_main.tex")
)
modelsummary(
  iv_models,
  coef_map = c(
    "fit_risk"  = "Institutions (risk) [IV]",
    "lat"       = "Latitude",
    "asia"      = "Asia",
    "africa"    = "Africa",
    "other"     = "Other",
    "edes1975"  = "% European descent",
    "malaria"   = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  output   = here("output", "tables", "tab_iv_main.html")
)
cat("IV table saved.\n")

# =============================================================================
# 3.  OLS vs. IV comparison figure
# =============================================================================

ols_models <- tryCatch(readRDS(here("output", "ols_models.rds")),
                        error = function(e) NULL)

extract_beta <- function(mod, nm_pattern, method) {
  coefs <- coef(mod)
  ses   <- se(mod)
  nm    <- names(coefs)[str_detect(names(coefs), nm_pattern)]
  if (length(nm) == 0) return(NULL)
  tibble(
    method   = method,
    estimate = coefs[nm[1]],
    se       = ses[nm[1]]
  )
}

iv_betas <- map2_dfr(
  iv_models, seq_along(iv_models),
  ~extract_beta(.x, "fit_risk|risk", sprintf("IV spec %d", .y))
)

ols_betas <- if (!is.null(ols_models)) {
  map2_dfr(
    ols_models, seq_along(ols_models),
    ~extract_beta(.x, "^risk$", sprintf("OLS spec %d", .y))
  )
} else NULL

comparison <- bind_rows(ols_betas, iv_betas) %>%
  mutate(
    ci_lo  = estimate - 1.96 * se,
    ci_hi  = estimate + 1.96 * se,
    type   = if_else(str_detect(method, "^IV"), "IV (2SLS)", "OLS (biased)"),
    spec_n = str_extract(method, "\\d+") %>% as.integer()
  )

p_compare <- comparison %>%
  ggplot(aes(x = estimate, y = fct_rev(factor(method)),
             color = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("IV (2SLS)" = "#2F3D70", "OLS (biased)" = "#D04E59"),
    name   = NULL
  ) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  labs(
    title    = "OLS vs. IV: Effect of Institutions on Log GDP per Capita",
    subtitle = paste0(
      "OLS (red): endogenous — biased upward by reverse causality\n",
      "IV (blue): corrects endogeneity using settler mortality as instrument"
    ),
    x        = "Coefficient on institutional quality (risk)",
    y        = NULL,
    caption  = "Specifications 1–7 as in table2.do / table3.do. HC1 robust SEs."
  ) +
  theme(
    legend.position = "none",
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
  )

ggsave(p_compare,
       filename = here("output", "figures", "fig_ols_vs_iv.pdf"),
       width = 10, height = 7, dpi = 300)
cat("OLS vs IV comparison figure saved.\n")

# =============================================================================
# 4.  IV coefficients stability figure
# =============================================================================

p_iv <- iv_betas %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    method = fct_rev(factor(method))
  ) %>%
  ggplot(aes(x = estimate, y = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, color = "#2F3D70", linewidth = 0.8) +
  geom_point(color = "#2F3D70", size = 3) +
  labs(
    title    = "IV/2SLS Stability: Effect of Institutions on Log GDP per Capita",
    subtitle = "7 specifications with different control variables | Instrument: log settler mortality",
    x        = "IV coefficient on expropriation risk",
    y        = NULL,
    caption  = "HC1 robust SEs. All specs point to large positive effect of institutions on income."
  ) +
  theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))

ggsave(p_iv,
       filename = here("output", "figures", "fig_iv_coefs.pdf"),
       width = 9, height = 5, dpi = 300)
cat("IV coefficient stability figure saved.\n")

saveRDS(iv_models, here("output", "iv_models.rds"))
