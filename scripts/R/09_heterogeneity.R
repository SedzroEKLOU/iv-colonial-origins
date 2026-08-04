# =============================================================================
# 09_heterogeneity.R
# =============================================================================
# Purpose : Regional heterogeneity and alternative institution measures.
#
# Two extensions:
#
#   A. REGIONAL SUB-SAMPLES
#      AJR's result may be driven by a specific region. We estimate the
#      IV separately for Africa, Asia, and Americas/Other. Neo-Europes
#      (AUS, CAN, NZL, USA) are excluded from sub-group analysis since
#      their settler mortality experience differs fundamentally.
#
#   B. ALTERNATIVE INSTITUTION MEASURES
#      AJR use average expropriation risk (1985-1995). But institutions
#      have multiple dimensions. The ajrcomment.dta includes:
#        - cons90  : constraints on the executive in 1900 (political institutions)
#        - lado1995: labour rights index 1995 (labour market institutions)
#      These are instrumentd the same way (logmort0) to check if the
#      institutional channel is specific to property rights or general.
#
# Outputs:
#   output/tables/tab_heterogeneity.tex
#   output/figures/fig_regional_iv.pdf
#   output/figures/fig_alt_institutions.pdf
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

region_colors <- c(
  "Africa"        = "#D04E59",
  "Asia"          = "#2F3D70",
  "Americas/Other"= "#BC8E7D",
  "Full sample"   = "#5B7C99"
)

# =============================================================================
# A.  Regional sub-sample IV
# =============================================================================

cat("Regional sub-sample IV estimates...\n")

# Remove neo-Europes from sub-group analysis
df_reg <- df %>% filter(neoeuro == 0, !is.na(logmort0), !is.na(risk), !is.na(loggdp))

run_regional_iv <- function(df_sub, region_label) {
  if (nrow(df_sub) < 5) {
    cat(sprintf("  %s: too few observations (%d), skipping\n", region_label, nrow(df_sub)))
    return(NULL)
  }
  mod <- tryCatch(
    feols(loggdp ~ 1 | risk ~ logmort0, data = df_sub, vcov = "HC1"),
    error = function(e) NULL
  )
  fs  <- tryCatch(
    feols(risk ~ logmort0, data = df_sub, vcov = "HC1"),
    error = function(e) NULL
  )
  f_on_inst <- tryCatch({
    wt <- wald(fs, "logmort0"); wt$stat
  }, error = function(e) NA_real_)

  if (is.null(mod)) return(NULL)
  coefs <- coef(mod)
  nm    <- names(coefs)[str_detect(names(coefs), "fit_risk")]
  tibble(
    region    = region_label,
    n         = nrow(df_sub),
    beta_iv   = coefs[nm[1]],
    se_iv     = se(mod)[nm[1]],
    fs_f      = f_on_inst
  )
}

regional_results <- bind_rows(
  run_regional_iv(df_reg, "Full sample (excl. neo-Eur.)"),
  run_regional_iv(df_reg %>% filter(africa == 1), "Africa"),
  run_regional_iv(df_reg %>% filter(asia   == 1), "Asia"),
  run_regional_iv(df_reg %>% filter(africa == 0, asia == 0), "Americas/Other")
) %>%
  mutate(
    ci_lo  = beta_iv - 1.96 * se_iv,
    ci_hi  = beta_iv + 1.96 * se_iv,
    region = factor(region, levels = c(
      "Full sample (excl. neo-Eur.)", "Africa", "Asia", "Americas/Other"
    ))
  )

cat("\n=== Regional IV results ===\n")
print(regional_results)

p_regional <- regional_results %>%
  ggplot(aes(x = beta_iv, y = fct_rev(region),
             color = region, size = n)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, linewidth = 0.8) +
  geom_point() +
  geom_text(aes(x = ci_hi + 0.2,
                label = sprintf("N=%d, F=%.1f", n, fs_f)),
            size = 3, hjust = 0, color = "grey30") +
  scale_color_manual(values = region_colors) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(
    title    = "Regional Heterogeneity: IV Estimate by Sub-Sample",
    subtitle = paste0(
      "2SLS estimate of institutions → log GDP per capita\n",
      "Neo-European countries excluded from all sub-groups"
    ),
    x        = "IV coefficient on institutions (risk → loggdp)",
    y        = NULL,
    color    = NULL,
    caption  = paste0(
      "95% CIs. HC1 robust SEs. F-statistic on logmort0 shown per region.\n",
      "Small sub-samples → wide CIs; interpret with caution."
    )
  ) +
  theme(
    legend.position = "none",
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
  )

ggsave(p_regional,
       filename = here("output", "figures", "fig_regional_iv.pdf"),
       width = 10, height = 5, dpi = 300)
cat("Regional figure saved.\n")

# =============================================================================
# B.  Alternative institution measures
# =============================================================================

cat("\nAlternative institution measures...\n")

df_alt <- df %>% filter(!is.na(logmort0), !is.na(loggdp))

# Baseline: expropriation risk (risk, 1985-1995)
iv_risk <- feols(loggdp ~ 1 | risk ~ logmort0,
                 data = df_alt %>% filter(!is.na(risk)), vcov = "HC1")

# Constraints on executive in 1900 (cons90)
iv_cons90 <- tryCatch(
  feols(loggdp ~ 1 | cons90 ~ logmort0,
        data = df_alt %>% filter(!is.na(cons90)), vcov = "HC1"),
  error = function(e) {
    cat("  Note: cons90 IV failed (possible collinearity). Trying with lat control.\n")
    feols(loggdp ~ lat | cons90 ~ logmort0,
          data = df_alt %>% filter(!is.na(cons90), !is.na(lat)), vcov = "HC1")
  }
)

# Labour rights (lado1995) — if sufficient variation
iv_lado <- tryCatch(
  feols(loggdp ~ 1 | lado1995 ~ logmort0,
        data = df_alt %>% filter(!is.na(lado1995)), vcov = "HC1"),
  error = function(e) NULL
)

# First stages for alternative instruments
fs_cons90 <- tryCatch(
  feols(cons90 ~ logmort0, data = df_alt %>% filter(!is.na(cons90)), vcov = "HC1"),
  error = function(e) NULL
)
fs_lado <- tryCatch(
  feols(lado1995 ~ logmort0, data = df_alt %>% filter(!is.na(lado1995)), vcov = "HC1"),
  error = function(e) NULL
)

alt_list <- list(
  "Expropriation risk (risk)"         = iv_risk,
  "Exec. constraints, 1900 (cons90)"  = iv_cons90
)
if (!is.null(iv_lado)) {
  alt_list[["Labour rights, 1995 (lado1995)"]] <- iv_lado
}

# Display
cat("\n=== Alternative institution IV results ===\n")
lapply(names(alt_list), function(nm) {
  mod   <- alt_list[[nm]]
  coefs <- coef(mod)
  ses   <- se(mod)
  endog <- names(coefs)[str_detect(names(coefs), "fit_")]
  if (length(endog) == 0) endog <- names(coefs)[2]
  cat(sprintf("  %s: β = %.3f (SE = %.3f), N = %d\n",
              nm, coefs[endog[1]], ses[endog[1]], nobs(mod)))
})

modelsummary(
  alt_list,
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title    = "Alternative Institution Measures: IV/2SLS Results",
  notes    = paste0(
    "Instrument: log settler mortality (logmort0) in all specifications. HC1 robust SEs. ",
    "Expropriation risk: average 1985-1995 (AJR's main measure). ",
    "Exec. constraints in 1900: political dimension of institutions. ",
    "Labour rights 1995: labour market dimension."
  ),
  output = here("output", "tables", "tab_heterogeneity.tex")
)
modelsummary(
  alt_list,
  gof_map  = c("nobs", "r.squared"),
  stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  output   = here("output", "tables", "tab_heterogeneity.html")
)
cat("Heterogeneity tables saved.\n")

# =============================================================================
# C.  Scatter: logmort0 → cons90 (alternative first stage visual)
# =============================================================================

if (!is.null(fs_cons90)) {
  p_cons90 <- df %>%
    filter(!is.na(logmort0), !is.na(cons90)) %>%
    ggplot(aes(x = logmort0, y = cons90,
               color = region, label = shortnam)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "grey40", fill = "grey85", linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_text_repel(size = 2.3, max.overlaps = 12, segment.color = "grey60") +
    scale_color_manual(values = c(
      "Africa"         = "#D04E59",
      "Asia"           = "#2F3D70",
      "Neo-Europe"     = "#5B7C99",
      "Americas/Other" = "#BC8E7D",
      "Other"          = "#7C7189"
    )) +
    labs(
      title    = "Alternative First Stage: Settler Mortality → Exec. Constraints (1900)",
      subtitle = "Does settler mortality predict political institutions in 1900?",
      x        = "Log settler mortality",
      y        = "Constraints on executive, 1900 (cons90)",
      color    = NULL,
      caption  = "Source: Albouy (2012). Political institutions measured before modern economic development."
    ) +
    theme(
      legend.position = "bottom",
      plot.caption    = element_text(size = 8, color = "grey40", hjust = 0)
    )

  ggsave(p_cons90,
         filename = here("output", "figures", "fig_alt_institutions.pdf"),
         width = 9, height = 6, dpi = 300)
  cat("Alternative institutions figure saved.\n")
}

cat("\n=== Pipeline complete ===\n")
cat("All 9 scripts run successfully. Check output/ for tables and figures.\n")
