# =============================================================================
# 02_descriptive.R
# =============================================================================
# Purpose : Descriptive statistics and the core intuition of AJR (2001).
#
# Three key facts that motivate the IV strategy:
#
#   1. CORRELATION: settler mortality correlates negatively with institutions
#      (high mortality → bad institutions). This is the relevance condition.
#
#   2. CORRELATION: institutions correlate positively with GDP per capita.
#      But this is OLS — endogenous (rich countries invest in institutions).
#
#   3. EXCLUSION: settler mortality should affect GDP ONLY through institutions,
#      not directly. This is untestable but defended via the history:
#      settlers who survived built lasting institutions; those who didn't
#      set up extractive ones (plantation colonies, mining enclaves).
#
# The three scatter plots below tell the whole story visually.
#
# Outputs:
#   output/figures/fig_mortality_institutions.pdf
#   output/figures/fig_institutions_gdp.pdf
#   output/figures/fig_mortality_gdp_reduced.pdf
#   output/figures/fig_three_panels.pdf
#   output/tables/tab_descriptive.tex
# =============================================================================

library(tidyverse)
library(gt)
library(ggthemes)
library(patchwork)
library(ggrepel)
library(modelsummary)
library(here)

cat("Loading data...\n")
df <- readRDS(here("data", "clean", "ajr_clean.rds"))

theme_set(
  theme_clean() +
    theme(
      plot.background   = element_blank(),
      legend.background = element_rect(color = "white"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title        = element_text(size = 12),
      axis.text         = element_text(size = 10)
    )
)

region_colors <- c(
  "Africa"          = "#D04E59",
  "Asia"            = "#2F3D70",
  "Neo-Europe"      = "#5B7C99",
  "Americas/Other"  = "#BC8E7D",
  "Other"           = "#7C7189"
)

# =============================================================================
# 1.  Summary statistics table
# =============================================================================

sum_vars <- c("loggdp", "risk", "logmort0", "lat",
              "africa", "asia", "neoeuro", "edes1975")

datasummary(
  All(df %>% select(all_of(sum_vars))) ~ Mean + SD + Min + Max + N,
  data   = df,
  title  = "Descriptive Statistics: AJR (2001) Sample",
  notes  = "N = 64 former colonies. logmort0 = log settler mortality. risk = expropriation risk (0–10). loggdp = log GDP per capita (PPP, 1995).",
  output = here("output", "tables", "tab_descriptive.tex")
)
datasummary(
  All(df %>% select(all_of(sum_vars))) ~ Mean + SD + Min + Max + N,
  data   = df,
  output = here("output", "tables", "tab_descriptive.html")
)
cat("Descriptive table saved.\n")

# =============================================================================
# 2.  Figure A — Settler mortality → Institutions (first stage intuition)
# =============================================================================

p_mort_inst <- df %>%
  filter(!is.na(logmort0), !is.na(risk)) %>%
  ggplot(aes(x = logmort0, y = risk, color = region, label = shortnam)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 15, segment.color = "grey60") +
  scale_color_manual(values = region_colors) +
  labs(
    title    = "First Stage: Settler Mortality → Institutional Quality",
    subtitle = paste0(
      "Colonies with high settler mortality developed extractive institutions (low risk score)\n",
      "Colonies where Europeans survived built inclusive institutions (high risk score)"
    ),
    x        = "Log settler mortality (per 1,000 per annum)",
    y        = "Protection against expropriation risk (0–10)",
    color    = NULL,
    caption  = "Source: Acemoglu, Johnson & Robinson (2001). Each point is a former colony."
  ) +
  theme(legend.position = "bottom",
        plot.caption    = element_text(size = 8, color = "grey40", hjust = 0))

ggsave(p_mort_inst,
       filename = here("output", "figures", "fig_mortality_institutions.pdf"),
       width = 9, height = 6, dpi = 300)
cat("Figure 1 (mortality → institutions) saved.\n")

# =============================================================================
# 3.  Figure B — Institutions → GDP (OLS, endogenous)
# =============================================================================

p_inst_gdp <- df %>%
  filter(!is.na(risk), !is.na(loggdp)) %>%
  ggplot(aes(x = risk, y = loggdp, color = region, label = shortnam)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 15, segment.color = "grey60") +
  scale_color_manual(values = region_colors) +
  labs(
    title    = "OLS (Endogenous): Institutions → Log GDP per Capita",
    subtitle = paste0(
      "Positive correlation — but OLS is biased: rich countries can afford better institutions\n",
      "Reverse causality + omitted variables → IV needed"
    ),
    x        = "Protection against expropriation risk (0–10)",
    y        = "Log GDP per capita (PPP, 1995)",
    color    = NULL,
    caption  = "OLS estimate is upward-biased. AJR (2001) instrument: settler mortality."
  ) +
  theme(legend.position = "bottom",
        plot.caption    = element_text(size = 8, color = "grey40", hjust = 0))

ggsave(p_inst_gdp,
       filename = here("output", "figures", "fig_institutions_gdp.pdf"),
       width = 9, height = 6, dpi = 300)
cat("Figure 2 (institutions → GDP) saved.\n")

# =============================================================================
# 4.  Figure C — Settler mortality → GDP (reduced form)
# =============================================================================

p_mort_gdp <- df %>%
  filter(!is.na(logmort0), !is.na(loggdp)) %>%
  ggplot(aes(x = logmort0, y = loggdp, color = region, label = shortnam)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey40", fill = "grey85", linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text_repel(size = 2.5, max.overlaps = 15, segment.color = "grey60") +
  scale_color_manual(values = region_colors) +
  labs(
    title    = "Reduced Form: Settler Mortality → Log GDP per Capita",
    subtitle = paste0(
      "ITT: colonies with higher settler mortality are poorer today\n",
      "IV estimate = Reduced form / First stage"
    ),
    x        = "Log settler mortality (per 1,000 per annum)",
    y        = "Log GDP per capita (PPP, 1995)",
    color    = NULL,
    caption  = "Reduced form: direct effect of the instrument on the outcome (valid under exclusion restriction)."
  ) +
  theme(legend.position = "bottom",
        plot.caption    = element_text(size = 8, color = "grey40", hjust = 0))

ggsave(p_mort_gdp,
       filename = here("output", "figures", "fig_mortality_gdp_reduced.pdf"),
       width = 9, height = 6, dpi = 300)
cat("Figure 3 (mortality → GDP, reduced form) saved.\n")

# =============================================================================
# 5.  Combined three-panel figure (publication-ready)
# =============================================================================

p_combined <- (p_mort_inst | p_mort_gdp) / p_inst_gdp +
  plot_annotation(
    title   = "The AJR (2001) IV Strategy: Three Key Relationships",
    subtitle = paste0(
      "Left: Instrument relevance (Z → D) | ",
      "Right: Reduced form (Z → Y) | ",
      "Bottom: OLS/IV second stage (D → Y)"
    ),
    caption = paste0(
      "Data: Acemoglu, Johnson & Robinson (2001), augmented by Albouy (2012). N = 64 former colonies.\n",
      "Z = log settler mortality | D = expropriation risk | Y = log GDP per capita (PPP, 1995)."
    ),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey30"),
      plot.caption  = element_text(size = 8,  color = "grey40", hjust = 0)
    )
  )

ggsave(p_combined,
       filename = here("output", "figures", "fig_three_panels.pdf"),
       width = 14, height = 12, dpi = 300)
cat("Combined three-panel figure saved.\n")

# =============================================================================
# 6.  Correlation matrix of key variables
# =============================================================================

cat("\n=== Correlation matrix (key variables) ===\n")
df %>%
  select(logmort0, risk, loggdp, lat, africa, asia) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(3) %>%
  print()
