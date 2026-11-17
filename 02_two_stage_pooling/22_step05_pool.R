# =====================================================================
# 22_step05_pool.R  (v1.0, 2026-11-06)
# Topic 6 - STEP05 two-stage pooling (Stage 2)
#
# Input : 04_结果/STEP04_库内模型_v1.0.csv  (in-cohort MI/CC estimates)
# Rules (frozen v0.3 8.2 / v1.0 10.1; user rulings 2026-11-06):
#   - Pool MI caliber only (main); CC kept for S6 later
#   - REML tau2 + Knapp-Hartung test + prediction interval (always reported)
#   - Construct exclusion: MHAS has no HI/DSI (absent rows);
#     HAALSI HI/DSI are binary-construct -> directional only, NOT in HI/DSI main pools (v1.0 6.2)
#   - k<5 cells: descriptive only (no pooled estimate; per-cohort range reported)
#   - LOO: for pooled cells, re-pool dropping each cohort (HR range)
#   - Sparse/degenerate rows (NA HR in STEP04) are excluded; k counts estimable cohorts
# Output: 04_结果/STEP05_合并_v1.0.csv
#         04_结果/STEP05_森林图_合并汇总.pdf/.png  (pooled diamonds, exposure x model x outcome)
#         04_结果/STEP05_森林图_逐库M3.pdf/.png    (per-cohort + pooled, M3)
#         04_结果/STEP05_森林图_CKM分层M3.pdf/.png (stratified pooled, M3)
# ASCII-only source (GBK console).
# =====================================================================

suppressWarnings(suppressMessages({ library(data.table); library(metafor) }))

SEED <- 20260904
ROOT <- normalizePath(".", winslash = "/")
D_OUT <- file.path(ROOT, "04_结果"); D_QC <- file.path(ROOT, "05_质控")
LOGF <- file.path(D_QC, "STEP05_run_log_20261106.txt")
logf <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "
", file = LOGF, append = TRUE), silent = TRUE)
unlink(LOGF); cat("STEP05 log ->", LOGF, "\n")

logf("metafor", as.character(packageVersion("metafor")))

res <- fread(file.path(D_OUT, "STEP04_库内模型_v1.0.csv"))
res[, se := (log(ucl) - log(lcl)) / (2 * 1.96)]
mi <- res[caliber == "mi" & !is.na(HR)]
logf("[1] STEP04 rows = ", nrow(res), "; MI estimable rows = ", nrow(mi),
     "; sparse/degenerate rows excluded = ", nrow(res) - nrow(mi) - sum(res$caliber == "cc"))

# construct exclusion: HAALSI HI/DSI binary (directional only, frozen v1.0 6.2)
mi[, in_main_pool := !(cohort == "HAALSI" & exposure %in% c("e_hi", "e_ds"))]

EXPLAB <- c(e_vi = "VI-only", e_hi = "HI-only", e_ds = "DSI")
pool_cell <- function(d, label) {
  d <- d[!is.na(HR) & !is.na(se) & se > 0]
  k <- nrow(d)
  if (k == 0) {
    return(data.table(cell = label, k = 0L, cohorts = NA_character_,
      pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
      p = NA_real_, I2 = NA_real_, tau2 = NA_real_, Q = NA_real_, Q_p = NA_real_,
      pi_lo = NA_real_, pi_hi = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
      desc_min_HR = NA_real_, desc_max_HR = NA_real_, note = "no_estimable_rows"))
  }
  out <- data.table(cell = label, k = k,
                    cohorts = paste(sort(d$cohort), collapse = "+"))
  if (k < 5) {
    dmin <- min(d$HR); dmax <- max(d$HR)
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, tau2 = NA_real_, Q = NA_real_, Q_p = NA_real_,
               pi_lo = NA_real_, pi_hi = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
               desc_min_HR = dmin, desc_max_HR = dmax,
               note = "k_lt5_descriptive")]
    return(out)
  }
  r <- tryCatch(rma(yi = log(HR), sei = se, data = d, method = "REML", test = "knha"),
                error = function(e) e)
  if (inherits(r, "error")) {
    out[, `:=`(pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
               p = NA_real_, I2 = NA_real_, tau2 = NA_real_, Q = NA_real_, Q_p = NA_real_,
               pi_lo = NA_real_, pi_hi = NA_real_, loo_min = NA_real_, loo_max = NA_real_,
               desc_min_HR = min(HR), desc_max_HR = max(HR), note = paste0("rma_failed: ", conditionMessage(r)))]
    return(out)
  }
  pi <- predict(r)
  out[, `:=`(pooled_HR = exp(r$b[1]),
             pooled_lcl = exp(r$ci.lb[1]), pooled_ucl = exp(r$ci.ub[1]),
             p = r$pval[1], I2 = r$I2, tau2 = r$tau2, Q = r$QE, Q_p = r$QEp,
             pi_lo = exp(pi$pi.lb), pi_hi = exp(pi$pi.ub),
             desc_min_HR = NA_real_, desc_max_HR = NA_real_, note = NA_character_)]
  # leave-one-out
  if (k >= 5) {
    loo <- vapply(seq_len(k), function(i) {
      tryCatch({
        ri <- rma(yi = log(HR), sei = se, data = d[-i], method = "REML", test = "knha")
        exp(ri$b[1])
      }, error = function(e) NA_real_)
    }, numeric(1))
    out[, `:=`(loo_min = min(loo, na.rm = TRUE), loo_max = max(loo, na.rm = TRUE))]
  }
  out
}

pool_all <- function(d, stratum_name) {
  out <- list()
  for (oc in c("death", "stroke")) {
    for (me in c("Cox", "FineGray", "CauseSpecificCox")) {
      for (ex in c("e_vi", "e_hi", "e_ds")) {
        for (mo in c("M1", "M2", "M3")) {
          sel <- d[outcome == oc & method == me & exposure == ex & model == mo & in_main_pool == TRUE]
          if (nrow(sel) == 0) next
          out[[length(out) + 1]] <- pool_cell(
            sel, paste(stratum_name, oc, me, EXPLAB[ex], mo, sep = "|"))
        }
      }
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

# ---- 1. overall pools (main: Cox death / FG stroke) ----
overall <- pool_all(mi[substratum == "overall"], "overall")
# ---- 2. CKM strata pools (M3, per stage) ----
strata <- rbindlist(lapply(c("ckm0_1", "ckm2", "ckm3"), function(stg) {
  pool_all(mi[substratum == stg & model == "M3"], stg)
}), use.names = TRUE, fill = TRUE)

allp <- rbindlist(list(overall[, method := NULL], strata), use.names = TRUE, fill = TRUE)
setnames(allp, "cell", "cell_id")
# split cell_id into columns
parts <- tstrsplit(allp$cell_id, "|", fixed = TRUE)
allp[, `:=`(substratum = parts[[1]], outcome = parts[[2]], method = parts[[3]],
            exposure = parts[[4]], model = parts[[5]])]
allp <- allp[order(substratum, outcome, method, exposure, model)]
fwrite(allp, file.path(D_OUT, "STEP05_合并_v1.0.csv"))
logf("[2] pooled rows = ", nrow(allp),
     " (k>=5 pooled: ", sum(!is.na(allp$pooled_HR)),
     "; k<5 descriptive: ", sum(allp$note == "k_lt5_descriptive", na.rm = TRUE), ")")

# =====================================================================
# Forest plots (working figures; final paper figures come in STEP08)
# =====================================================================
png_dev <- function(path, w, h) png(path, width = w, height = h, units = "px", res = 300)
pdf_dev <- function(path, w, h) pdf(path, width = w, height = h, onefile = TRUE)

# ---- Fig A: pooled summary (overall, Cox death + FG stroke, M1-M3 x 3 exposures) ----
summ <- allp[substratum == "overall" & ((outcome == "death" & method == "Cox") |
                                        (outcome == "stroke" & method == "FineGray"))]
summ[, panel := fifelse(outcome == "death", "All-cause death (Cox)", "Incident stroke (Fine-Gray)")]
summ[, rowlab := paste0(EXPLAB[exposure], " | ", model)]
make_forest_summary <- function(dev_open, path_pdf, path_png) {
  panels <- c("All-cause death (Cox)", "Incident stroke (Fine-Gray)")
  models <- c("M1", "M2", "M3")
  labs <- unlist(lapply(panels, function(pn)
    c(pn, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
      paste0("  ", ex, " | ", models))))))
  est <- unlist(lapply(panels, function(pn) c(NA, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pooled_HR[1]) && r$pooled_HR[1] < 5 && r$pooled_HR[1] > 0.2) r$pooled_HR[1] else NA
    }))))))
  lo <- unlist(lapply(panels, function(pn) c(NA, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pooled_lcl[1]) && r$pooled_HR[1] < 5 && r$pooled_HR[1] > 0.2) r$pooled_lcl[1] else NA
    }))))))
  hi <- unlist(lapply(panels, function(pn) c(NA, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pooled_ucl[1]) && r$pooled_HR[1] < 5 && r$pooled_HR[1] > 0.2) r$pooled_ucl[1] else NA
    }))))))
  pi_l <- unlist(lapply(panels, function(pn) c(NA, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pi_lo[1]) && is.finite(r$pi_lo[1]) && r$pi_lo[1] > 0.05) r$pi_lo[1] else NA
    }))))))
  pi_h <- unlist(lapply(panels, function(pn) c(NA, unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pi_hi[1]) && is.finite(r$pi_hi[1]) && r$pi_hi[1] < 20) r$pi_hi[1] else NA
    }))))))
  txt <- unlist(lapply(panels, function(pn) c("", unlist(lapply(c("VI-only", "HI-only", "DSI"), function(ex)
    lapply(models, function(mo) {
      r <- summ[panel == pn & rowlab == paste0(ex, " | ", mo)]
      if (nrow(r) && !is.na(r$pooled_HR[1])) {
        if (!is.na(r$note[1]) && r$note[1] == "k_lt5_descriptive") "k<5" else
          sprintf("%.2f (%.2f-%.2f)", r$pooled_HR[1], r$pooled_lcl[1], r$pooled_ucl[1])
      } else if (nrow(r)) "-" else ""
    }))))))
  n <- length(labs)
  dev_open; par(mar = c(4.5, 4, 2, 6))
  plot(1, 1, type = "n", log = "x", xlim = c(0.3, 4), ylim = c(0.5, n + 1),
       xlab = "Pooled HR (REML+KH, 95% CI)", yaxt = "n", ylab = "")
  abline(v = 1, lty = 3, col = "grey60")
  for (i in seq_len(n)) {
    y <- n + 1 - i
    if (labs[i] %in% panels) {
      text(0.3, y, labs[i], font = 2, cex = 0.8, adj = 0)
    } else {
      text(0.3, y, labs[i], cex = 0.7, adj = 0)
      if (!is.na(est[i])) {
        arrows(pi_l[i], y, pi_h[i], y, angle = 90, length = 0.02, code = 3, col = "grey70")
        segments(lo[i], y, hi[i], y, lwd = 2)
        points(est[i], y, pch = 18, cex = 1.1)
      }
      text(4.6, y, txt[i], cex = 0.7, adj = 1)
    }
  }
}
make_forest_summary(pdf_dev(file.path(D_OUT, "STEP05_森林图_合并汇总.pdf"), 8, 10),
                    NULL, NULL); dev.off()
png_dev(file.path(D_OUT, "STEP05_森林图_合并汇总.png"), 2400, 3000)
make_forest_summary(NULL, NULL, NULL); dev.off()
logf("[3] fig summary written")

# ---- Fig B: per-cohort + pooled, M3 ----
plot_percohort <- function(df_pool, df_cohort, panels, fname) {
  pdf(file.path(D_OUT, paste0(fname, ".pdf")), width = 8, height = 9)
  draw_percohort(df_pool, df_cohort, panels); dev.off()
  png_dev(file.path(D_OUT, paste0(fname, ".png")), 2400, 2700)
  draw_percohort(df_pool, df_cohort, panels); dev.off()
}
draw_percohort <- function(df_pool, df_cohort, panels) {
  rows <- list()
  for (pn in panels) {
    oc <- if (grepl("death", pn)) "death" else "stroke"
    me <- if (oc == "death") "Cox" else "FineGray"
    rows[[length(rows) + 1]] <- data.table(lab = pn, panel = pn, is_hdr = TRUE,
                                           HR = NA, lcl = NA, ucl = NA)
    for (ex in c("e_vi", "e_hi", "e_ds")) {
      cd <- df_cohort[outcome == oc & method == me & exposure == ex & model == "M3"]
      cd <- cd[order(cohort)]
      for (i in seq_len(nrow(cd))) {
        rows[[length(rows) + 1]] <- data.table(lab = paste0("    ", cd$cohort[i]),
          panel = pn, is_hdr = FALSE, HR = cd$HR[i], lcl = cd$lcl[i], ucl = cd$ucl[i])
      }
      pd <- df_pool[substratum == "overall" & outcome == oc & method == me &
                      exposure == ex & model == "M3"]
      if (nrow(pd) && !is.na(pd$pooled_HR[1])) {
        rows[[length(rows) + 1]] <- data.table(lab = paste0("    Pooled (k=", pd$k[1], ")"),
          panel = pn, is_hdr = FALSE, HR = pd$pooled_HR[1],
          lcl = pd$pooled_lcl[1], ucl = pd$pooled_ucl[1], pooled = TRUE)
      } else {
        rows[[length(rows) + 1]] <- data.table(lab = "    Pooled (k<5)", panel = pn,
          is_hdr = FALSE, HR = NA, lcl = NA, ucl = NA, pooled = TRUE)
      }
    }
  }
  dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  n <- nrow(dt)
  rng <- range(c(dt$lcl, dt$ucl), na.rm = TRUE); rng <- pmin(pmax(rng, 0.2), 5)
  par(mar = c(4.5, 10, 2, 4))
  plot(1, 1, type = "n", log = "x", xlim = rng, ylim = c(0.5, n + 1),
       xlab = "HR (M3, MI caliber)", yaxt = "n", ylab = "")
  abline(v = 1, lty = 3, col = "grey60")
  for (i in seq_len(n)) {
    y <- n + 1 - i
    if (dt$is_hdr[i]) { text(0.2, y, dt$lab[i], font = 2, cex = 0.8, adj = 0) }
    else if (is.na(dt$HR[i])) { text(0.2, y, dt$lab[i], cex = 0.7, adj = 0); text(3, y, "k<5", cex = 0.7) }
    else {
      text(0.2, y, dt$lab[i], cex = if (grepl("Pooled", dt$lab[i])) 0.72 else 0.7,
           font = if (grepl("Pooled", dt$lab[i])) 2 else 1, adj = 0)
      segments(dt$lcl[i], y, dt$ucl[i], y, lwd = if (grepl("Pooled", dt$lab[i])) 3 else 1.6)
      points(dt$HR[i], y, pch = if (grepl("Pooled", dt$lab[i])) 18 else 16,
             cex = if (grepl("Pooled", dt$lab[i])) 1.3 else 0.8)
    }
  }
}
mi_m3 <- mi[model == "M3" & substratum == "overall"]
pool_m3 <- allp[substratum == "overall" & model == "M3"]
plot_percohort(pool_m3, mi_m3,
               c("All-cause death (Cox, M3)", "Incident stroke (Fine-Gray, M3)"),
               "STEP05_森林图_逐库M3")
logf("[4] fig per-cohort written")

# ---- Fig C: CKM strata (M3, pooled where k>=5) ----
draw_strata <- function() {
  panels <- c("All-cause death (Cox, by CKM stage)", "Incident stroke (Fine-Gray, by CKM stage)")
  rows <- list()
  for (pn in panels) {
    oc <- if (grepl("death", pn)) "death" else "stroke"
    me <- if (oc == "death") "Cox" else "FineGray"
    rows[[length(rows) + 1]] <- data.table(lab = pn, is_hdr = TRUE, HR = NA, lcl = NA, ucl = NA, k = NA)
    for (stg in c("ckm0_1", "ckm2", "ckm3")) {
      for (ex in c("e_vi", "e_hi", "e_ds")) {
        pd <- allp[substratum == stg & outcome == oc & method == me & exposure == ex & model == "M3"]
        if (nrow(pd) == 0) next
        lab <- sprintf("    %s | %s", sub("ckm", "Stage ", stg), EXPLAB[ex])
        if (!is.na(pd$pooled_HR[1])) {
          rows[[length(rows) + 1]] <- data.table(lab = lab, is_hdr = FALSE,
            HR = pd$pooled_HR[1], lcl = pd$pooled_lcl[1], ucl = pd$pooled_ucl[1], k = pd$k[1])
        } else {
          rows[[length(rows) + 1]] <- data.table(lab = paste0(lab, " (k=", pd$k[1], ")"),
            is_hdr = FALSE, HR = NA, lcl = NA, ucl = NA, k = pd$k[1])
        }
      }
    }
  }
  dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  n <- nrow(dt)
  rng <- range(c(dt$lcl, dt$ucl), na.rm = TRUE); rng <- pmin(pmax(rng, 0.2), 5)
  par(mar = c(4.5, 12, 2, 4))
  plot(1, 1, type = "n", log = "x", xlim = rng, ylim = c(0.5, n + 1),
       xlab = "Pooled HR (M3, by CKM stage)", yaxt = "n", ylab = "")
  abline(v = 1, lty = 3, col = "grey60")
  for (i in seq_len(n)) {
    y <- n + 1 - i
    if (dt$is_hdr[i]) { text(0.2, y, dt$lab[i], font = 2, cex = 0.8, adj = 0) }
    else if (is.na(dt$HR[i])) { text(0.2, y, dt$lab[i], cex = 0.7, adj = 0); text(3, y, "k<5", cex = 0.7) }
    else {
      text(0.2, y, dt$lab[i], cex = 0.7, adj = 0)
      segments(dt$lcl[i], y, dt$ucl[i], y, lwd = 2)
      points(dt$HR[i], y, pch = 18, cex = 1.1)
    }
  }
}
pdf(file.path(D_OUT, "STEP05_森林图_CKM分层M3.pdf"), width = 8, height = 9); draw_strata(); dev.off()
png_dev(file.path(D_OUT, "STEP05_森林图_CKM分层M3.png"), 2400, 2700); draw_strata(); dev.off()
logf("[5] fig strata written")
logf("[DONE] STEP05 complete")
cat("STEP05 DONE\n")
