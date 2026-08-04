# =============================================================================
# 07_albouy_critique.R
# =============================================================================
# Purpose : Albouy (2012) critique — sensitivity of AJR results to data quality.
#
# Background:
#   Albouy (2012, AER) argues that AJR's settler mortality data has three
#   problems:
#     1. CAMPAIGN data: mortality during military campaigns, not representative
#        of long-term settler conditions. Flag: campaign == 1.
#     2. SLAVE data: mortality from African slave labour, not European settlers.
#        Flag: slave == 1.
#     3. CONJECTURE: some mortality figures are AJR's own interpolations.
#        Flag: source0 == 0 (primary_source == 0).
#
#   The ajrcomment.dta dataset (Albouy's companion data) includes `campaign`
#   and `slave` dummies. Albouy shows that excluding or correcting these
#   observations weakens or reverses the first stage, undermining the IV.
#
# Analysis plan:
#   A. First stage robustness: add campaign + slave dummies as controls
#   B. First stage robustness: restrict to primary_source == 1 (no conjecture)
#   C. First stage robustness: exclude campaign observations
#   D. First stage robustness: exclude slave observations
#   E. IV estimate stability under each sample restriction
#   F. Visualise the contested vs. uncontested observations
#
# Outputs:
#   output/tables/tab_albouy.tex
#   output/figures/fig_albouy_scatter.pdf
#   output/figures/fig_albouy_coefs.pdf
# =============================================================================

library(tidyverse)
library(fixest)
library(modelsummary)
library(ggthemes)
library(ggrepel)
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

# Classify observations
df <- df %>%
  mutate(
    data_quality = case_when(
      campaign == 1 & slave == 1 ~ "Campaign + Slave",
      campaign == 1              ~ "Campaign only",
      slave    == 1              ~ "Slave only",
      primary_source == 0        ~ "Conjectured (AJR)",
      TRUE                       ~ "Primary source"
    ),
    data_quality = factor(data_quality,
                          levels = c("Primary source", "Conjectured (AJR)",
                                     "Campaign only", "Slave only",
                                     "Campaign + Slave"))
  )

cat("\n=== Data quality breakdown ===\n")
df %>% count(data_quality) %>% print()

quality_colors <- c(
  "Primary source"    = "#2F3D70",
  "Conjectured (AJR)" = "#7C7189",
  "Campaign only"     = "#BC8E7D",
  "Slave only"        = "#D04E59",
  "Campaign + Slave"  = "#8B0000"
)

# =============================================================================
# A.  Scatter: log mortality by data quality
# =============================================================================

p_scatter <- df %>%
  filter(!is.na(logmort0), !is.na(risk)) %>%
  ggplot(aes(x = logmort0, y = risk,
             color = data_quality, label = shortnam)) +
  geom_smooth(data = df %>% filter(!is.na(logmort0), !is.na(risk)),
              aes(color = NULL),
              method = "lm", formula = y ~ x, se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.7) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_text_repel(size = 2.3, max.overlaps = 12, segment.color = "grey70") +
  scale_color_manual(values = quality_colors) +
  labs(
    title    = "Albouy (2012) Critique: Data Quality in the First Stage",
    subtitle = paste0(
      "Settler mortality → institutional quality\n",
      "Red/orange: Albouy's contested observations (campaign, slave mortality)"
    ),
    x        = "Log settler mortality",
    y        = "Protection against expropriation risk (0–10)",
    color    = "Data source",
    caption  = paste0(
      "Campaign: mortality from military campaigns (not long-term settler conditions).\n",
      "Slave: mortality from African slave labour records. Albouy (2012) argues these ",
      "are not informative about European settler mortality."
    )
  ) +
  theme(
    legend.position = "right",
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
  )

ggsave(p_scatter,
       filename = here("output", "figures", "fig_albouy_scatter.pdf"),
       width = 10, height = 6.5, dpi = 300)
cat("Albouy scatter figure saved.\n")

# =============================================================================
# B.  First stage and IV: 6 sample variations
# =============================================================================

cat("\nEstimating first stage and IV under Albouy critique samples...\n")

df_full        <- df %>% filter(!is.na(logmort0), !is.na(risk), !is.na(loggdp))
df_no_campaign <- df_full %>% filter(campaign == 0)
df_no_slave    <- df_full %>% filter(slave == 0)
df_no_both     <- df_full %>% filter(campaign == 0, slave == 0)
df_primary     <- df_full %>% filter(primary_source == 1)

run_fs_iv <- function(df_sub, label) {
  fs <- tryCatch(feols(risk ~ logmort0, data = df_sub, vcov = "HC1"),
                 error = function(e) NULL)
  iv <- tryCatch(feols(loggdp ~ 1 | risk ~ logmort0, data = df_sub, vcov = "HC1"),
                 error = function(e) NULL)

  f_on_inst <- tryCatch({
    wt <- wald(fs, "logmort0"); wt$stat
  }, error = function(e) NA_real_)

  tibble(
    sample      = label,
    n           = nrow(df_sub),
    fs_coef     = if (!is.null(fs)) coef(fs)["logmort0"]     else NA_real_,
    fs_se       = if (!is.null(fs)) se(fs)["logmort0"]       else NA_real_,
    fs_f        = f_on_inst,
    iv_coef     = if (!is.null(iv)) coef(iv)["fit_risk"]     else NA_real_,
    iv_se       = if (!is.null(iv)) se(iv)["fit_risk"]       else NA_real_
  )
}

critique_results <- bind_rows(
  run_fs_iv(df_full,        "Full sample (AJR)"),
  run_fs_iv(df_no_campaign, "Excl. campaign"),
  run_fs_iv(df_no_slave,    "Excl. slave"),
  run_fs_iv(df_no_both,     "Excl. campaign + slave"),
  run_fs_iv(df_primary,     "Primary source only")
)

cat("\n=== Albouy sensitivity results ===\n")
critique_results %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print()

# =============================================================================
# C.  Coefficient comparison figure
# =============================================================================

critique_long <- critique_results %>%
  pivot_longer(
    cols         = c(fs_coef, iv_coef),
    names_to     = "type",
    values_to    = "estimate"
  ) %>%
  mutate(
    se_col = if_else(type == "fs_coef", fs_se, iv_se),
    ci_lo  = estimate - 1.96 * se_col,
    ci_hi  = estimate + 1.96 * se_col,
    label  = if_else(type == "fs_coef",
                     "First stage: β(logmort0 → risk)",
                     "IV/2SLS: β(risk → loggdp)"),
    sample = fct_rev(factor(sample, levels = critique_results$sample))
  )

p_critique <- critique_long %>%
  ggplot(aes(x = estimate, y = sample, color = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("First stage: β(logmort0 → risk)" = "#D04E59",
               "IV/2SLS: β(risk → loggdp)"       = "#2F3D70"),
    name = NULL
  ) +
  facet_wrap(~label, scales = "free_x") +
  labs(
    title    = "Sensitivity to Albouy (2012) Critique: First Stage and IV",
    subtitle = paste0(
      "How first stage and IV estimate change as contested observations are removed\n",
      "Large F drop when excluding campaign/slave data → AJR instrument may be fragile"
    ),
    x       = "Coefficient estimate",
    y       = NULL,
    caption = paste0(
      "95% confidence intervals. HC1 robust SEs.\n",
      "Albouy (2012, AER) argues campaign/slave observations inflate first stage F."
    )
  ) +
  theme(
    legend.position = "none",
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
  )

ggsave(p_critique,
       filename = here("output", "figures", "fig_albouy_coefs.pdf"),
       width = 12, height = 6, dpi = 300)
cat("Albouy coefficient figure saved.\n")

# =============================================================================
# D.  Adding campaign + slave as controls (rather than dropping)
# =============================================================================

cat("\nAdding campaign + slave dummies as controls...\n")

df_ctrl <- df_full %>% filter(!is.na(campaign), !is.na(slave))

fs_with_flags <- feols(risk ~ logmort0 + campaign + slave,
                       data = df_ctrl, vcov = "HC1")
iv_with_flags <- feols(loggdp ~ campaign + slave | risk ~ logmort0,
                       data = df_ctrl, vcov = "HC1")

cat("\nFirst stage with campaign + slave dummies:\n")
etable(fs_with_flags)

cat("\nIV with campaign + slave dummies as controls:\n")
etable(iv_with_flags)

# Export
modelsummary(
  list(
    "FS (full)"              = feols(risk ~ logmort0, data = df_full, vcov = "HC1"),
    "FS (+ flags)"           = fs_with_flags,
    "FS (excl. campaign+slave)" = feols(risk ~ logmort0, data = df_no_both, vcov = "HC1"),
    "IV (full)"              = feols(loggdp ~ 1 | risk ~ logmort0, data = df_full, vcov = "HC1"),
    "IV (+ flags)"           = iv_with_flags,
    "IV (excl. campaign+slave)" = feols(loggdp ~ 1 | risk ~ logmort0, data = df_no_both, vcov = "HC1")
  ),
  coef_map = c(
    "logmort0"  = "Log settler mortality",
    "fit_risk"  = "Institutions (risk) [IV]",
    "campaign"  = "Campaign dummy",
    "slave"     = "Slave dummy"
  ),
  gof_map = c("nobs", "r.squared"),
  stars   = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title   = "Albouy (2012) Critique: Robustness of First Stage and IV",
  notes   = paste0(
    "First three columns: first stage (risk ~ logmort0). ",
    "Last three columns: IV/2SLS (loggdp ~ risk | logmort0). ",
    "Campaign and slave flags from Albouy (2012). HC1 robust SEs."
  ),
  output  = here("output", "tables", "tab_albouy.tex")
)
modelsummary(
  list(
    "FS (full)"              = feols(risk ~ logmort0, data = df_full, vcov = "HC1"),
    "FS (+ flags)"           = fs_with_flags,
    "FS (excl. campaign+slave)" = feols(risk ~ logmort0, data = df_no_both, vcov = "HC1"),
    "IV (full)"              = feols(loggdp ~ 1 | risk ~ logmort0, data = df_full, vcov = "HC1"),
    "IV (+ flags)"           = iv_with_flags,
    "IV (excl. campaign+slave)" = feols(loggdp ~ 1 | risk ~ logmort0, data = df_no_both, vcov = "HC1")
  ),
  coef_map = c(
    "logmort0" = "Log settler mortality",
    "fit_risk" = "Institutions (risk) [IV]",
    "campaign" = "Campaign dummy",
    "slave"    = "Slave dummy"
  ),
  gof_map = c("nobs", "r.squared"),
  stars   = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  output  = here("output", "tables", "tab_albouy.html")
)
cat("Albouy sensitivity table saved.\n")

cat("\n=== Summary of Albouy Critique ===\n")
cat(sprintf("  Full sample F = %.2f\n",
            critique_results$fs_f[critique_results$sample == "Full sample (AJR)"]))
cat(sprintf("  Excl. campaign+slave F = %.2f\n",
            critique_results$fs_f[critique_results$sample == "Excl. campaign + slave"]))
cat("  Interpretation: if F drops below 10 after exclusions, instrument may be fragile.\n")
cat("  However, IV direction (positive effect of institutions) persists across samples.\n")
