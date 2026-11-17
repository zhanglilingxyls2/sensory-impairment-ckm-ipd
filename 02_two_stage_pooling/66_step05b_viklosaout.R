# =====================================================================
# 66_step05b_viklosaout.R  (v1.0, 2026-11-14)
# Topic 6 - VI-only sensitivity EXCLUDING KLoSA (no blood assays; staged
# without laboratory components). Owner ruling 2026-11-14 ("做 VI-only 剔
# KLoSA 敏感性; 暂缓 HRS 2024 死亡侧").
#
# Question answered: does the cohort staged without laboratory components
# drive the VI stage pattern / interaction? All stage cells remain above
# the k>=5 pooling gate (overall 6; stage 0-1 5, stage 2 6, stage 3 6 for
# death; stroke stage 0-1 was k=4 descriptive and stays descriptive).
#
# Pooling machinery mirrors 22_step05_pool.R / 23_step06_interaction_
# subgroups.R exactly: MI caliber, REML tau2 + Knapp-Hartung test,
# t-based CI, prediction interval; HR blocks on log scale, RERI on the
# linear (delta-method) scale; k<5 cells descriptive ranges.
#
# Input : 04_结果/STEP04_库内模型_v1.0.csv
#         04_结果/STEP06_逐库_交互与亚组_v1.0.csv
# Output: 04_结果/STEP05b_VI_KLoSAout_v1.0.csv
# ASCII-only source (GBK console).
# =====================================================================
suppressWarnings(suppressMessages({ library(data.table); library(metafor) }))

ROOT <- normalizePath(".", winslash = "/")
D_OUT <- file.path(ROOT, "04_结果"); D_QC <- file.path(ROOT, "05_质控")
LOGF <- file.path(D_QC, "STEP05b_VIKLoSAout_log_20261114.txt")
logf <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "
", file = LOGF, append = TRUE), silent = TRUE)
unlink(LOGF); cat("STEP05b log ->", LOGF, "\n")

res <- fread(file.path(D_OUT, "STEP04_库内模型_v1.0.csv"))
res[, se := (log(ucl) - log(lcl)) / (2 * 1.96)]
mi <- res[caliber == "mi" & !is.na(HR)]
d6 <- fread(file.path(D_OUT, "STEP06_逐库_交互与亚组_v1.0.csv"))
# Fine-Gray interaction rows carry compound exposure names; normalize as in 23
d6[block == "interaction" & grepl("_stg", exposure),
   `:=`(level = sub(".*_stg", "stg", exposure), exposure = sub("_stg.*", "", exposure))]
d6 <- d6[!(cohort == "HAALSI" & exposure %in% c("e_hi", "e_ds"))]

pool_hr <- function(d, label) {
  d <- d[!is.na(HR) & !is.na(se) & se > 0]
  k <- nrow(d)
  out <- data.table(cell_id = label, k = k,
                    cohorts = paste(sort(d$cohort), collapse = "+"))
  if (k < 5) {
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
               desc_min = min(d$HR), desc_max = max(d$HR), note = "k_lt5_descriptive")]
    return(out)
  }
  r <- tryCatch(rma(yi = log(HR), sei = se, data = d, method = "REML", test = "knha"),
                error = function(e) e)
  if (inherits(r, "error")) {
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
               desc_min = min(HR), desc_max = max(HR), note = paste0("rma_failed"))]
    return(out)
  }
  pi <- predict(r)
  out[, `:=`(pooled_HR = exp(r$b[1]), pooled_lcl = exp(r$ci.lb[1]), pooled_ucl = exp(r$ci.ub[1]),
             p = r$pval[1], I2 = r$I2, pi_lo = exp(pi$pi.lb), pi_hi = exp(pi$pi.ub),
             desc_min = NA_real_, desc_max = NA_real_, note = NA_character_)]
  out
}

pool_stage2b <- function(d, scale) {
  d <- copy(d)[!is.na(est) & !is.na(lcl) & !is.na(ucl)]
  d[, se := if (scale == "log") (log(ucl) - log(lcl)) / (2 * 1.96) else (ucl - lcl) / (2 * 1.96)]
  d <- d[is.finite(se) & se > 0 & is.finite(est)]
  d[, yi := if (scale == "log") log(est) else est]
  k <- nrow(d)
  out <- data.table(cell_id = paste("s2", d$block[1], d$outcome[1], d$method[1],
                                    d$exposure[1], d$level[1], sep = "|"),
                    k = k, cohorts = paste(sort(d$cohort), collapse = "+"))
  if (k < 3) {
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
               desc_min = mean(d$est), desc_max = NA_real_, note = "k_lt3_descriptive")]
    return(out)
  }
  r <- tryCatch(rma(yi = yi, sei = se, data = as.data.frame(d), method = "REML", test = "knha"),
                error = function(e) e)
  if (inherits(r, "error")) {
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
               desc_min = NA_real_, desc_max = NA_real_, note = "rma_failed")]
    return(out)
  }
  pi <- predict(r)
  tr <- if (scale == "log") exp else identity
  out[, `:=`(pooled_HR = tr(r$b[1]), pooled_lcl = tr(r$ci.lb[1]), pooled_ucl = tr(r$ci.ub[1]),
             p = r$pval[1], I2 = r$I2, pi_lo = tr(pi$pi.lb), pi_hi = tr(pi$pi.ub),
             desc_min = NA_real_, desc_max = NA_real_, note = NA_character_)]
  out
}

out <- list()
# ---- 1. overall VI-only (M1-M3), death Cox / stroke FG ----
for (oc in c("death", "stroke")) {
  me <- if (oc == "death") "Cox" else "FineGray"
  for (mo in c("M1", "M2", "M3")) {
    sel <- mi[outcome == oc & method == me & exposure == "e_vi" & model == mo &
              substratum == "overall" & cohort != "KLoSA"]
    out[[length(out) + 1]] <- pool_hr(sel, paste("overall", oc, me, "VI-only", mo, sep = "|"))
  }
}
# ---- 2. stage-stratified M3 VI-only ----
for (stg in c("ckm0_1", "ckm2", "ckm3")) {
  for (oc in c("death", "stroke")) {
    me <- if (oc == "death") "Cox" else "FineGray"
    sel <- mi[outcome == oc & method == me & exposure == "e_vi" & model == "M3" &
              substratum == stg & cohort != "KLoSA"]
    out[[length(out) + 1]] <- pool_hr(sel, paste(stg, oc, me, "VI-only", "M3", sep = "|"))
  }
}
# ---- 3. interaction (VI x stage2 / stage3) and RERI, VI-only ----
ix <- d6[block == "interaction" & exposure == "e_vi" & cohort != "KLoSA"]
for (oc in c("death", "stroke")) {
  me <- if (oc == "death") "Cox" else "FineGray"
  for (lv in c("stg2", "stg3")) {
    sel <- ix[outcome == oc & method == me & level == lv]
    if (nrow(sel)) out[[length(out) + 1]] <- pool_stage2b(sel, "log")
  }
}
rer <- d6[block == "RERI" & exposure == "e_vi" & cohort != "KLoSA" & outcome == "death"]
if (nrow(rer)) out[[length(out) + 1]] <- pool_stage2b(rer, "linear")

fin <- rbindlist(out, use.names = TRUE, fill = TRUE)
fwrite(fin, file.path(D_OUT, "STEP05b_VI_KLoSAout_v1.0.csv"))
logf("[DONE] cells = ", nrow(fin),
     "; pooled = ", sum(!is.na(fin$pooled_HR)),
     "; descriptive = ", sum(!is.na(fin$note) & grepl("descriptive", fin$note)))
cat("STEP05b DONE\n")
