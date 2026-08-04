# iv-colonial-origins

**IV Replication: Acemoglu, Johnson & Robinson (2001)**  
*Colonial Origins of Comparative Development: An Empirical Investigation*  
American Economic Review, 91(5): 1369–1401

---

## Overview

This repository replicates and extends the landmark AJR (2001) paper using R. AJR argue that colonial institutions — shaped by historical settler mortality — explain the vast majority of income differences across former colonies today. The causal chain is:

```
Settler mortality → Institutional quality → GDP per capita
     (Z)                    (D)                   (Y)
```

The core identification strategy is **instrumental variables (IV / 2SLS)**:
- **Endogeneity problem**: Rich countries can afford better institutions — OLS is biased upward by reverse causality.
- **Instrument**: Log settler mortality per 1,000 per annum. Where European settlers faced high mortality (tropical disease environments), they built extractive institutions. Where they survived, they built inclusive institutions that persist today.
- **Exclusion restriction**: Settler mortality affected GDP only through institutional quality, not directly (defended via historical mechanism — 18th/19th-century disease environments have no direct effect on modern GDP once institutions are controlled for).

---

## Data

**Source**: Albouy (2012) companion data — `ajrcomment.dta`  
Download: [AER Data Archive](https://www.aeaweb.org/articles?id=10.1257/aer.102.6.3059) (article 10.1257/aer.102.6.3059)

Place the file at `data/raw/ajrcomment.dta`.

**Sample**: 64 former colonies, cross-section.

**Key variables**:

| Variable   | Description |
|------------|-------------|
| `loggdp`   | Log GDP per capita (PPP, 1995) — outcome |
| `risk`     | Avg. protection against expropriation risk, 1985-1995 (0–10) — endogenous institutions |
| `logmort0` | Log settler mortality per 1,000 per annum — instrument |
| `lat`      | Abs. latitude / 90 — geography control |
| `africa`, `asia`, `other` | Continent dummies |
| `neoeuro`  | 1 = neo-European settler colony (AUS, CAN, NZL, USA) |
| `edes1975` | % European descent in 1975 |
| `malaria`  | Malaria index |
| `campaign` | Albouy: 1 = mortality from military campaign |
| `slave`    | Albouy: 1 = mortality from slave labour records |
| `cons90`   | Constraints on executive, 1900 (alternative institutions) |
| `lado1995` | Labour rights index, 1995 (alternative institutions) |

---

## Repository Structure

```
iv-colonial-origins/
├── data/
│   ├── raw/          ← place ajrcomment.dta here
│   └── clean/        ← ajr_clean.rds (generated)
├── scripts/
│   └── R/
│       ├── 00_master.R           Master script: runs full pipeline
│       ├── 01_data_prep.R        Load + clean data
│       ├── 02_descriptive.R      Descriptive stats + 3-panel IV diagram
│       ├── 03_ols_naive.R        OLS (biased baseline, 6 specs)
│       ├── 04_first_stage.R      First stage (7 specs, F-stats)
│       ├── 05_iv_2sls.R          IV/2SLS main results (7 specs)
│       ├── 06_weak_iv_tests.R    Anderson-Rubin CIs + weak IV diagnostics
│       ├── 07_albouy_critique.R  Albouy (2012) data quality sensitivity
│       ├── 08_robustness.R       LIML, JIVE, leave-one-out influence
│       └── 09_heterogeneity.R    Regional sub-samples + alt. institutions
├── output/
│   ├── tables/       ← .tex and .html tables
│   └── figures/      ← .pdf figures
└── README.md
```

---

## How to Run

```r
# Install required packages (once)
pkgs <- c("tidyverse", "fixest", "ivreg", "modelsummary", "gt",
          "ggthemes", "ggrepel", "patchwork", "sandwich", "haven", "here")
install.packages(pkgs)

# Run full pipeline
source("scripts/R/00_master.R")
```

The master script runs all 9 scripts in sequence and logs completion. Total runtime: ~2-3 minutes (dominated by the Anderson-Rubin grid search in script 06).

---

## Estimation Strategy

### 1. OLS (biased)
```
loggdp_i = α + β·risk_i + X_i'γ + ε_i
```
OLS is biased upward: rich countries invest in institutions, so Cov(risk, ε) > 0.

### 2. First Stage
```
risk_i = α + γ·logmort0_i + X_i'δ + u_i
```
High settler mortality → extractive institutions (γ < 0). Relevance condition: F(logmort0) > 10 (Staiger-Stock rule of thumb).

### 3. IV/2SLS (main result)
```
loggdp_i = α + β_IV·risk̂_i + X_i'γ + ε_i
```
β_IV = reduced form / first stage = Cov(logmort0, loggdp) / Cov(logmort0, risk).

### 4. Weak IV inference
Anderson-Rubin (1949) confidence intervals: valid regardless of instrument strength. Grid search over β₀ ∈ [-3, 15]; collect all β₀ where AR test cannot be rejected at 5%.

### 5. Albouy critique
Albouy (2012) argues that campaign and slave mortality are not representative of settler mortality. Scripts 07 tests sensitivity to excluding or flagging these observations.

### 6. Extensions
- **LIML**: median-unbiased under weak IV
- **JIVE**: eliminates finite-sample bias (Angrist et al. 1999)
- **Leave-one-out**: identifies high-leverage colonies
- **Regional sub-samples**: Africa, Asia, Americas
- **Alternative institutions**: executive constraints (cons90), labour rights (lado1995)

---

## Key Results (expected)

| | β(institutions) | F-stat (1st stage) |
|---|---|---|
| OLS | ~0.52 | — |
| IV (baseline) | ~0.94 | ~22 |
| IV (+ latitude) | ~1.00 | ~17 |
| LIML | ~0.94 | — |
| JIVE | ~0.95 | — |

OLS understates the causal effect in some interpretations (measurement error attenuation) but is biased in an unknown direction overall. The IV estimate — accounting for reverse causality — is larger in AJR's data.

---

## R Packages

| Package | Use |
|---------|-----|
| `fixest` | Fast OLS, IV/2SLS (feols), F-stats (wald) |
| `ivreg` | LIML estimation |
| `modelsummary` | Publication-quality regression tables |
| `ggthemes` + `ggrepel` | Figures |
| `sandwich` | HC1 SEs for base R models (JIVE) |
| `haven` | Read Stata .dta files |
| `here` | Portable paths |

---

## References

Acemoglu, D., Johnson, S., & Robinson, J. A. (2001). The colonial origins of comparative development: An empirical investigation. *American Economic Review*, 91(5), 1369–1401.

Albouy, D. Y. (2012). The colonial origins of comparative development: An empirical investigation: Comment. *American Economic Review*, 102(6), 3059–3076.

Anderson, T. W., & Rubin, H. (1949). Estimation of the parameters of a single equation in a complete system of stochastic equations. *Annals of Mathematical Statistics*, 20(1), 46–63.

Angrist, J. D., Imbens, G. W., & Krueger, A. B. (1999). Jackknife instrumental variables estimation. *Journal of Applied Econometrics*, 14(1), 57–67.

Staiger, D., & Stock, J. H. (1997). Instrumental variables regression with weak instruments. *Econometrica*, 65(3), 557–586.

Andrews, I., Stock, J. H., & Sun, L. (2019). Weak instruments in instrumental variables regression: Theory and practice. *Annual Review of Economics*, 11, 727–753.
