# Core analysis code: Sensory impairment and cardiovascular–kidney–metabolic syndrome stages

Code accompanying the manuscript *"Sensory impairment and the risk of death and stroke across cardiovascular–kidney–metabolic syndrome stages 0–3: a two-stage individual-participant-data meta-analysis of seven international cohorts"* (submitted to BMC Medicine).

## Scope of this repository

This repository contains the **core statistical pipeline** — cohort-level models, two-stage pooling — and the code that produces the **main-text figures and tables** (Figure 1, Figure 2, Table 1). It does **not** contain participant data, aggregate result files, or the data-preparation and sensitivity/validation code, which are available from the corresponding author upon reasonable request.

## Overview

Two-stage individual-participant-data (IPD) meta-analysis. In stage 1, cohort-specific Cox proportional-hazards models (all-cause death) and Fine–Gray competing-risk models (stroke, with death as competing event) were fitted in seven cohorts (CHARLS, ELSA, HRS, KLoSA, MHAS, SHARE, HAALSI). In stage 2, cohort estimates were pooled with random-effects meta-analysis (restricted maximum-likelihood between-study variance, Knapp–Hartung adjustment).

## Software

- R 4.6.1 — packages: `survival`, `cmprsk`, `survey`, `metafor`
- Scripts read and write local files with absolute paths; adjust the path constants at the top of each script to your own directory layout.

## Pipeline order

Leading numbers reflect the chronological build order of the project.

| Folder | Script | Role |
|---|---|---|
| `01_cohort_models/` | `20_step04_cohort_models.R` | Stage-1 Cox models, all-cause death endpoint (M1–M3) |
| | `29_step04b_cvd_endpoint.R` | Stage-1 stroke endpoint, Fine–Gray with death as competing event |
| | `23_step06_interaction_subgroups.R` | Stage-specific models, subgroups and pooled multiplicative interactions (Table 3, Figure 2) |
| `02_two_stage_pooling/` | `22_step05_pool.R` | Random-effects pooling (REML + Knapp–Hartung), main results (Table 2, Figure 1) |
| | `31_step07_repool.R` | Re-pooling for sensitivity and interaction sets |
| | `66_step05b_viklosaout.R` | Pooled estimates excluding KLoSA (cohort staged without laboratory components) |
| `03_main_figures_tables/` | `35_table1_baseline.R` | Table 1 (baseline characteristics) |
| | `34_step08_figures.R` | Main forest plots, Figures 1 and 2 |

Main-text Tables 2–5 report the pooled estimates produced by the scripts above (Tables 2–3 from `22`/`31` with `20`/`29`/`23`; Table 4 from the mediation analysis and Table 5 from the sensitivity analyses, whose code is available on request).

## Data access

No participant data are included. All cohort data must be obtained from their custodians:

- CHARLS: https://charls.charlsdata.com/pages/data/111/en.html
- ELSA: https://www.elsa-project.ac.uk/accessing-elsa-data
- SHARE: https://share-eric.eu/data/data-access
- KLoSA: https://survey.keis.or.kr/eng/myinfo/login.jsp
- MHAS: https://mhasweb.org/DataProducts/Home.aspx
- HAALSI: https://haalsi.org
- HRS: https://hrsdata.isr.umich.edu/data-products
- NHANES (public): https://wwwn.cdc.gov/nchs/nhanes/
- UKHLS (public): https://www.understandingsociety.ac.uk

Harmonized cross-cohort variables were built with the Gateway to Global Aging Data harmonized datasets (https://g2aging.org).

## License

Released under the MIT License (see `LICENSE`).
