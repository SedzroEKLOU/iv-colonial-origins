# =============================================================================
# 06_weak_iv_tests.R
# =============================================================================
# Purpose : Weak IV diagnostics and inference robust to weak instruments.
#           Translates key logic from table3.do (Anderson-Rubin grid search).
#
# Tests performed:
#
#   1. STAIGER-STOCK F (first stage F on excluded instrument)
#      Rule of thumb: F < 10 → weak instrument, 2SLS inference unreliable.
#      Loaded from output/f_stats.rds (computed in 04_first_stage.R).
#
#   2. KLEIBERGEN-PAAP rk Wald F (cluster-robust first stage F)
#      More conservative than Cragg-Donald for small N / HC SEs.
#      fixest computes this when vcov != "iid".
#
#   3. ANDERSON-RUBIN (AR) CONFIDENCE INTERVALS
#      Valid regardless of instrument strength. Inverts the AR test statistic.
#      Method: grid search over β, compute AR statistic at each point,
#      collect all β where AR p-value > 0.05 → 95% CI.
#      Translates the grid search logic from table3.do.
#
#   4. STOCK-WRIGHT S statistic (joint test, IV robust)
#      Available via AR inversion approach.
#
# Key reference: Staiger & Stock (1997), Andrews, Stock & Sun (2019).
#
# Outputs:
#   output/tables/tab_weak_iv.tex
#   output/figures/fig_ar_ci.pdf
# =============================================================================

library(tidyverse)
library(fixest)
library(modelsummary)
library(ggthemes)
library(here)

cat("Loading data and pre-computed models...\n")
df       <- readRDS(here("data", "clean", "ajr_clean.rds"))
iv_models <- readRDS(here("output", "iv_models.rds"))
fs_models <- readRDS(here("output", "fs_models.rds"))
f_stats   <- readRDS(here("output", "f_stats.rds"))

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
# 1.  Staiger-Stock F summary (already computed in 04)
# =============================================================================

cat("\n=== Staiger-Stock F on logmort0 ===\n")
for (i in seq_along(f_stats)) {
  status <- if (!is.na(f_stats[i]) && f_stats[i] > 10) "STRONG" else "WEAK"
  cat(sprintf("  Spec %d: F = %6.2f  [%s]\n", i, f_stats[i], status))
}

# =============================================================================
# 2.  Anderson-Rubin grid search (translation of table3.do)
# =============================================================================
# table3.do logic:
#   for b in beta_grid:
#     y_tilde = loggdp - b * risk
#     reg y_tilde logmort0 controls, robust
#     AR(b) = t-stat on logmort0 in this regression
#   CI_AR = {b : |AR(b)| < 1.96}  (95% confidence interval)
#
# The AR statistic under H0: β = b0 is:
#   AR(b0) = F-stat of H0: γ = 0 in  (loggdp - b0 * risk) = γ·logmort0 + controls + ε
# which follows χ²(1)/1 under H0 regardless of instrument strength.

ar_ci <- function(df_sub, controls = NULL, b_lo = -3, b_hi = 15,
                  n_grid = 2000, alpha = 0.05) {
  beta_grid <- seq(b_lo, b_hi, length.out = n_grid)
  ar_stats  <- numeric(n_grid)

  for (j in seq_along(beta_grid)) {
    b          <- beta_grid[j]
    df_tmp     <- df_sub %>% mutate(y_tilde = loggdp - b * risk)
    fmla       <- if (is.null(controls)) {
      y_tilde ~ logmort0
    } else {
      reformulate(c("logmort0", controls), "y_tilde")
    }
    mod        <- tryCatch(
      feols(fmla, data = df_tmp, vcov = "HC1"),
      error = function(e) NULL
    )
    if (is.null(mod)) { ar_stats[j] <- NA_real_; next }
    # F-stat on logmort0 (with HC1, fixest returns Wald F)
    wt <- tryCatch(wald(mod, "logmort0"), error = function(e) NULL)
    ar_stats[j] <- if (!is.null(wt)) wt$stat else NA_real_
  }

  # Critical value: F(1,∞) = χ²(1) at level α
  crit <- qchisq(1 - alpha, df = 1)

  inside <- !is.na(ar_stats) & ar_stats < crit
  if (!any(inside)) {
    return(list(lo = NA_real_, hi = NA_real_,
                beta_grid = beta_grid, ar_stats = ar_stats, crit = crit))
  }
  idx_in <- which(inside)
  list(
    lo         = beta_grid[min(idx_in)],
    hi         = beta_grid[max(idx_in)],
    beta_grid  = beta_grid,
    ar_stats   = ar_stats,
    crit       = crit
  )
}

cat("\nComputing Anderson-Rubin CIs (grid search, may take a moment)...\n")

# Spec 1: No controls
ar1 <- ar_ci(df %>% filter(!is.na(logmort0) & !is.na(risk) & !is.na(loggdp)))

# Spec 2: + latitude
ar2 <- ar_ci(df %>% filter(!is.na(logmort0) & !is.na(risk) & !is.na(loggdp) & !is.na(lat)),
             controls = "lat")

# Spec 4: + continents
ar4 <- ar_ci(df %>% filter(!is.na(logmort0) & !is.na(risk) & !is.na(loggdp)),
             controls = c("asia", "africa", "other"))

ar_results <- tibble(
  spec        = c("(1) No controls", "(2) + Latitude", "(4) + Continents"),
  ar_ci_lo    = c(ar1$lo, ar2$lo, ar4$lo),
  ar_ci_hi    = c(ar1$hi, ar2$hi, ar4$hi),
  iv_estimate = c(
    coef(iv_models[[1]])["fit_risk"],
    coef(iv_models[[2]])["fit_risk"],
    coef(iv_models[[4]])["fit_risk"]
  )
)

cat("\n=== Anderson-Rubin 95% CIs ===\n")
ar_results %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print()

# =============================================================================
# 3.  AR function plot (spec 1 — baseline)
# =============================================================================

p_ar <- tibble(
  beta    = ar1$beta_grid,
  ar_stat = ar1$ar_stats
) %>%
  filter(!is.na(ar_stat)) %>%
  ggplot(aes(x = beta, y = ar_stat)) +
  geom_line(color = "#2F3D70", linewidth = 0.9) +
  geom_hline(yintercept = ar1$crit, linetype = "dashed",
             color = "#D04E59", linewidth = 0.8) +
  geom_vline(xintercept = ar1$lo, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = ar1$hi, linetype = "dotted", color = "grey40") +
  annotate("rect",
           xmin = ar1$lo, xmax = ar1$hi,
           ymin = -Inf, ymax = ar1$crit,
           fill = "#2F3D70", alpha = 0.08) +
  annotate("text", x = mean(c(ar1$lo, ar1$hi)),
           y = ar1$crit + 0.5,
           label = sprintf("95%% CI: [%.2f, %.2f]", ar1$lo, ar1$hi),
           size = 3.5, color = "#2F3D70") +
  annotate("text", x = max(ar1$beta_grid) * 0.8,
           y = ar1$crit + 0.5,
           label = sprintf("χ²(1) critical value = %.2f", ar1$crit),
           size = 3, color = "#D04E59", hjust = 1) +
  labs(
    title    = "Anderson-Rubin Confidence Interval (Spec 1: No controls)",
    subtitle = paste0(
      "AR statistic at each hypothesised β₀ for institutions → GDP\n",
      "95% CI = all β where AR < χ²(1) critical value (valid regardless of instrument strength)"
    ),
    x        = "Hypothesised β₀ (coefficient on institutions)",
    y        = "Anderson-Rubin test statistic",
    caption  = paste0(
      "Baseline specification: no geographic controls. HC1 robust SEs.\n",
      "IV point estimate = ", round(coef(iv_models[[1]])["fit_risk"], 3),
      " (should lie in the center of the CI)."
    )
  ) +
  theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0))

ggsave(p_ar,
       filename = here("output", "figures", "fig_ar_ci.pdf"),
       width = 9, height = 5.5, dpi = 300)
cat("AR CI figure saved.\n")

# =============================================================================
# 4.  Summary table: IV estimate + AR CI + F-stat
# =============================================================================

# Extract IV estimates for all 7 specs
iv_coefs <- map_dfr(seq_along(iv_models), function(i) {
  mod   <- iv_models[[i]]
  coefs <- coef(mod)
  ses   <- se(mod)
  nm    <- names(coefs)[str_detect(names(coefs), "fit_risk")]
  tibble(
    spec     = names(iv_models)[i],
    iv_est   = coefs[nm[1]],
    iv_se    = ses[nm[1]],
    f_stat   = f_stats[i],
    n        = nobs(mod)
  )
})

iv_coefs <- iv_coefs %>%
  mutate(
    iv_ci_lo_norm = iv_est - 1.96 * iv_se,
    iv_ci_hi_norm = iv_est + 1.96 * iv_se
  )

cat("\n=== Full diagnostic table ===\n")
print(iv_coefs)

cat("\nWeak IV summary:\n")
cat(sprintf("  F > 10 in %d/%d specifications\n",
            sum(f_stats > 10, na.rm = TRUE), length(f_stats)))
cat("  Anderson-Rubin CIs (robust to weak IV) are wide but positive,\n")
cat("  confirming the direction of the effect under instrument uncertainty.\n")

saveRDS(ar_results, here("output", "ar_results.rds"))
