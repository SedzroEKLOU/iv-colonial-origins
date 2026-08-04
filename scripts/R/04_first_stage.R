# =============================================================================
# 04_first_stage.R
# =============================================================================
# Purpose : First stage — instrument relevance check.
#           Replicates table2.do (7 specifications).
#
# First stage: risk ~ logmort0 + controls
#   The key question: does settler mortality PREDICT institutional quality?
#   It should (negative coefficient — high mortality → bad institutions).
#   The F-statistic on logmort0 must be large for the instrument to be
#   relevant (Staiger & Stock 1997 rule of thumb: F > 10).
#
# 7 specifications (from table2.do):
#   1. No controls
#   2. + Latitude
#   3. No Neo-Europes (exclude AUS, CAN, NZL, USA)
#   4. + Continent dummies (asia, africa, other)
#   5. + Continent dummies + latitude
#   6. + % European descent (edes1975)
#   7. + Malaria index
#
# Outputs:
#   output/tables/tab_first_stage.tex
#   output/figures/fig_first_stage_coefs.pdf
# =============================================================================

library(tidyverse)
library(fixest)
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
# 7 First stage specifications
# =============================================================================

cat("Estimating first stage (7 specifications)...\n")

# Clustered on logmort0 (as in the Stata code: cluster(logmort0))
# In practice with N=64, clustering on country is standard
# We use HC1 for comparability with the paper (small N, no natural cluster)

fs1 <- feols(risk ~ logmort0,
             data = df, vcov = "HC1")

fs2 <- feols(risk ~ logmort0 + lat,
             data = df, vcov = "HC1")

fs3 <- feols(risk ~ logmort0,
             data = df %>% filter(neoeuro == 0), vcov = "HC1")

fs4 <- feols(risk ~ logmort0 + asia + africa + other,
             data = df, vcov = "HC1")

fs5 <- feols(risk ~ logmort0 + asia + africa + other + lat,
             data = df, vcov = "HC1")

fs6 <- feols(risk ~ logmort0 + edes1975,
             data = df %>% filter(!is.na(edes1975)), vcov = "HC1")

fs7 <- feols(risk ~ logmort0 + malaria,
             data = df %>% filter(!is.na(malaria)), vcov = "HC1")

fs_models <- list(
  "(1) No controls"     = fs1,
  "(2) Latitude"        = fs2,
  "(3) No neo-Europes"  = fs3,
  "(4) Continents"      = fs4,
  "(5) Cont. + Lat."    = fs5,
  "(6) Euro descent"    = fs6,
  "(7) Malaria"         = fs7
)

# =============================================================================
# F-statistics on logmort0 (instrument relevance)
# =============================================================================

f_stats <- sapply(fs_models, function(mod) {
  # Wald test on logmort0
  tryCatch({
    wt <- wald(mod, "logmort0")
    wt$stat
  }, error = function(e) NA_real_)
})

cat("\n=== First stage: F-statistics on logmort0 ===\n")
for (i in seq_along(f_stats)) {
  cat(sprintf("  Spec %d: F = %.2f %s\n",
              i, f_stats[i],
              if (!is.na(f_stats[i]) && f_stats[i] > 10) "✓ (strong)" else "⚠ (weak?)"))
}

cat("\n=== First stage: coefficient on logmort0 ===\n")
for (i in seq_along(fs_models)) {
  mod <- fs_models[[i]]
  cat(sprintf("  Spec %d: β = %.3f (SE = %.3f)\n",
              i, coef(mod)["logmort0"], se(mod)["logmort0"]))
}

# =============================================================================
# Export table
# =============================================================================

modelsummary(
  fs_models,
  coef_map = c(
    "logmort0" = "Log settler mortality",
    "lat"      = "Latitude",
    "asia"     = "Asia",
    "africa"   = "Africa",
    "other"    = "Other",
    "edes1975" = "% European descent",
    "malaria"  = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title    = "First Stage: Settler Mortality → Institutional Quality (Table 2)",
  notes    = paste0(
    "Notes: OLS (first stage). Outcome: expropriation risk (0-10, higher = better institutions). ",
    "HC1 robust SEs. Spec (3) excludes neo-European countries (AUS, CAN, NZL, USA). ",
    "F-statistics on log settler mortality: ",
    paste(sprintf("%.1f", f_stats), collapse = ", "), "."
  ),
  output = here("output", "tables", "tab_first_stage.tex")
)
modelsummary(
  fs_models,
  coef_map = c(
    "logmort0" = "Log settler mortality",
    "lat"      = "Latitude",
    "asia"     = "Asia",
    "africa"   = "Africa",
    "other"    = "Other",
    "edes1975" = "% European descent",
    "malaria"  = "Malaria index"
  ),
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  output   = here("output", "tables", "tab_first_stage.html")
)
cat("First stage table saved.\n")

# =============================================================================
# Figure: First stage coefficient stability across specs
# =============================================================================

fs_coefs <- tibble(
  spec     = names(fs_models),
  estimate = sapply(fs_models, function(m) coef(m)["logmort0"]),
  se       = sapply(fs_models, function(m) se(m)["logmort0"]),
  n        = sapply(fs_models, nobs),
  f_stat   = f_stats
) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    spec  = fct_rev(factor(spec, levels = names(fs_models)))
  )

p_fs <- fs_coefs %>%
  ggplot(aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, color = "#2F3D70", linewidth = 0.8) +
  geom_point(aes(size = f_stat), color = "#2F3D70", alpha = 0.85) +
  geom_text(aes(x = ci_hi + 0.02, label = sprintf("F=%.1f", f_stat)),
            size = 3, hjust = 0, color = "grey30") +
  scale_size_continuous(range = c(2, 6), name = "F-stat") +
  labs(
    title    = "First Stage: Effect of Settler Mortality on Institutional Quality",
    subtitle = paste0(
      "Coefficient on log settler mortality across 7 specifications\n",
      "Point size = F-statistic on instrument | F > 10 indicates strong relevance"
    ),
    x        = "Coefficient on log settler mortality",
    y        = NULL,
    caption  = "Notes: HC1 robust SEs. Negative coefficient expected: high mortality → bad institutions."
  ) +
  theme(
    legend.position = "right",
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
  )

ggsave(p_fs,
       filename = here("output", "figures", "fig_first_stage_coefs.pdf"),
       width = 10, height = 5.5, dpi = 300)
cat("First stage figure saved.\n")

saveRDS(fs_models, here("output", "fs_models.rds"))
saveRDS(f_stats,   here("output", "f_stats.rds"))
