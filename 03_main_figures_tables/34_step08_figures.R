# 34_step08_figures.R - STEP08 论文版 Figure 1 / Figure 2（160 mm 印刷宽，Times，无内嵌标题，完整图内图例）
# Figure 1：暴露（VI only / HI only / DSI）× 终点（死亡 Cox / 卒中 Fine-Gray）× M1–M3 合并森林图（STEP05 overall）
# Figure 2：按 CKM 期（0–1 / 2 / 3）分层的暴露 HR（核心图；M3）；k<5 格画描述性范围（虚线，不解读）
#   右列附合并交互 P（STEP06：死亡 Cox / 卒中 Fine-Gray 的 stage3 交互项）
# 输出：04_结果/Figure1_论文版_森林图.pdf / Figure2_论文版_分期森林图.pdf（PNG 由 PDF 300dpi 栅格化，另步）
suppressWarnings(suppressMessages({ library(data.table) }))
ROOT <- normalizePath(".", winslash = "/")
D_OUT <- file.path(ROOT, "04_结果")
logf <- try(file(file.path(ROOT, "05_质控", "STEP08_fig_run_log_20261108.txt"), open = "wt", encoding = "UTF-8"), silent = TRUE)
lg <- function(...) { s <- paste0(...); cat(s, "\n"); if (!inherits(logf, "try-error")) try(cat(s, "\n", file = logf), silent = TRUE) }

p5 <- fread(file.path(D_OUT, "STEP05_合并_v1.0.csv"))
p6 <- fread(file.path(D_OUT, "STEP06_合并_交互与亚组_v1.0.csv"))

EXLBL <- c("VI-only" = "Vision impairment only", "HI-only" = "Hearing impairment only", "DSI" = "Dual sensory impairment")
EXORD <- c("VI-only", "HI-only", "DSI")
MODLBL <- c(M1 = "M1", M2 = "M2", M3 = "M3")
STGLBL <- c("ckm0_1" = "Stage 0\u20131", "ckm2" = "Stage 2", "ckm3" = "Stage 3")

fmt_p <- function(p) { if (is.na(p)) return("\u2014"); if (p < 0.001) return("<.001"); sub("^0", "", sprintf("%.3f", p)) }
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f\u2013%.2f)", h, l, u)
fmt_i2 <- function(x) { if (is.na(x)) return("\u2014"); sprintf("%.0f%%", x) }

# ---------------- 公共绘制函数 ----------------
draw_forest <- function(rows, xlim_hr, ticks, colspec, legend, path, w_in, axis_lab, unit_pt = 16) {
  H_rows <- 0
  for (r in rows) H_rows <- H_rows + if (!is.null(r$panel_hdr)) 1.15 else
    (1.0 + if (!is.null(r$grp) && !is.na(r$grp)) 0.85 else 0)
  total <- 2.2 + H_rows + 1.6 + length(legend) + 0.6
  # cairo_pdf: the WinAnsi pdf() device transliterates >= (U+2265) to "=" in the legend
  cairo_pdf(path, width = w_in, height = round(total * unit_pt / 72, 2), family = "Times")
  par(mar = c(0.1, 0.1, 0.1, 0.1), xpd = NA)
  plot(NA, xlim = c(0, 10), ylim = c(-(total - 2.2), 2.2), axes = FALSE, xlab = "", ylab = "")
  # 顶线 + 表头线
  segments(0.05, 1.9, 9.7, 1.9, lwd = 1.6)
  segments(0.05, 1.15, 9.7, 1.15, lwd = 1.0)
  for (i in seq_along(colspec$heads)) text(colspec$x[i], 1.55, colspec$heads[i], cex = 1.0, adj = c(0, 0.5), font = 3)
  text(colspec$x_axis_mid, 1.55, axis_lab, cex = 1.0, adj = c(0.5, 0.5))
  # 参考线 1 的 x 位置（行循环前算好）
  x1 <- colspec$x_axis_left + (log(1) - log(xlim_hr[1])) / (log(xlim_hr[2]) - log(xlim_hr[1])) * (colspec$x_axis_right - colspec$x_axis_left)
  # 行
  y <- 0.3
  for (r in rows) {
    if (!is.null(r$panel_hdr)) {
      text(0.05, y, r$panel_hdr, cex = 1.08, adj = c(0, 0.5), font = 2)
      y <- y - 1.15; next
    }
    if (!is.null(r$grp) && !is.na(r$grp)) {
      text(0.05, y, r$grp, cex = 1.0, adj = c(0, 0.5), font = 2)
      y <- y - 0.85
    }
    if (!is.null(r$sub) && !is.na(r$sub)) text(0.55, y, r$sub, cex = 0.98, adj = c(0, 0.5))
    if (!is.na(r$lcl) && !is.na(r$ucl) && !is.null(r$lcl)) {
      segments(x1, y - 0.38, x1, y + 0.38, lty = 3, col = "grey55")   # 参考线 1：仅数据行高度，标题行断开
      mfun <- function(v) colspec$x_axis_left + (log(v) - log(xlim_hr[1])) / (log(xlim_hr[2]) - log(xlim_hr[1])) * (colspec$x_axis_right - colspec$x_axis_left)
      xm <- if (!is.na(r$hr)) mfun(r$hr) else NA_real_; xl <- mfun(r$lcl); xr <- mfun(r$ucl)
      if (r$pooled) {
        segments(xl, y, xr, y, lwd = 2.0)
        rect(xm - 0.045, y - 0.16, xm + 0.045, y + 0.16, col = "black", border = "black")
      } else {
        segments(xl, y, xr, y, lwd = 1.6, lty = 2)   # 描述格：仅虚线范围，无点
      }
    }
    if (!is.na(r$s1)) text(colspec$x[1], y, r$s1, cex = 0.95, adj = c(0, 0.5))
    if (!is.na(r$s2)) text(colspec$x[2], y, r$s2, cex = 0.95, adj = c(0, 0.5))
    if (!is.na(r$s3)) text(colspec$x[3], y, r$s3, cex = 0.95, adj = c(0, 0.5))
    if (!is.na(r$s4)) text(colspec$x[4], y, r$s4, cex = 0.95, adj = c(0, 0.5))
    y <- y - 1.0
  }
  # 刻度（x1 已在行循环前算好）
  mfun1 <- function(v) colspec$x_axis_left + (log(v) - log(xlim_hr[1])) / (log(xlim_hr[2]) - log(xlim_hr[1])) * (colspec$x_axis_right - colspec$x_axis_left)
  tick_y <- y + 0.45
  for (tv in ticks) {
    segments(mfun1(tv), tick_y, mfun1(tv), tick_y - 0.16, lwd = 1)
    text(mfun1(tv), tick_y - 0.42, format(tv, drop0trailing = TRUE), cex = 0.85, adj = c(0.5, 0.5))
  }
  # 图例
  y_leg_top <- tick_y - 1.15
  segments(0.05, y_leg_top + 0.45, 9.7, y_leg_top + 0.45, lwd = 1.0)
  yy <- y_leg_top
  for (ln in legend) { text(0.05, yy, ln, cex = 0.75, adj = c(0, 0.5)); yy <- yy - 1.0 }
  dev.off()
  lg("written:", path)
}
# ---------------- Figure 1 ----------------
ov <- p5[substratum == "overall"]
f1_rows <- list()
i <- 1
for (oc in c("death", "stroke")) {
  f1_rows[[length(f1_rows) + 1]] <- list(panel_hdr = if (oc == "death") "A   All-cause death (Cox)" else "B   Incident stroke (Fine-Gray)")
  for (ex in EXORD) {
    f1_rows[[length(f1_rows) + 1]] <- list(grp = EXLBL[[ex]], sub = NA_character_,
      hr = NA_real_, lcl = NA_real_, ucl = NA_real_, s1 = NA_character_, s2 = NA_character_,
      s3 = NA_character_, s4 = NA_character_, pooled = TRUE)
    for (md in c("M1", "M2", "M3")) {
      mth <- if (oc == "death") "Cox" else "FineGray"
      r <- ov[outcome == oc & exposure == ex & model == md & method == mth]
      f1_rows[[length(f1_rows) + 1]] <- list(
        grp = NA_character_,
        sub = MODLBL[[md]],
        hr = if (nrow(r)) r$pooled_HR else NA_real_, lcl = if (nrow(r)) r$pooled_lcl else NA_real_, ucl = if (nrow(r)) r$pooled_ucl else NA_real_,
        s1 = if (nrow(r)) fmt_hr(r$pooled_HR, r$pooled_lcl, r$pooled_ucl) else "\u2014",
        s2 = if (nrow(r)) fmt_p(r$p) else "\u2014", s3 = if (nrow(r)) fmt_i2(r$I2) else "\u2014",
        s4 = if (nrow(r)) as.character(r$k) else "\u2014", pooled = TRUE)
    }
  }
}

leg1 <- c(
  "Pooled two-stage IPD meta-analysis (REML, Knapp-Hartung). Squares = pooled HRs; lines = 95% CIs.",
  "Exposures (reference: no sensory impairment): VI only = vision impairment without hearing impairment;",
  "HI only = hearing impairment without vision impairment; DSI = dual sensory impairment.",
  "M1: age, sex.  M2: additionally education, income below cohort median, wealth quintiles.",
  "M3: additionally smoking, alcohol, physical inactivity, cancer, chronic lung disease, arthritis,",
  "and antihypertensive medication (reduced covariate set in HAALSI).",
  "Death: Cox models on calendar time.  Stroke: Fine-Gray subdistribution hazards; death competing.",
  "k = cohorts contributing (VI: 7; HI/DSI: 5; MHAS no general hearing item; HAALSI hearing binary).",
  "I² = heterogeneity.  Reference line at HR = 1; axis on the log scale.")
cs1 <- list(x = c(5.9, 7.95, 8.72, 9.42), heads = c("HR (95% CI)", "P", "I²", "k"),
            x_axis_left = 2.95, x_axis_right = 5.55, x_axis_mid = 4.25)
draw_forest(f1_rows, c(0.75, 2.0), c(0.8, 1, 1.5, 2), cs1, leg1,
            file.path(D_OUT, "Figure1_论文版_森林图.pdf"), 6.3, "Hazard ratio (log scale)")

# ---------------- Figure 2 ----------------
stg <- p5[substratum %in% c("ckm0_1", "ckm2", "ckm3") & model == "M3"]
ix <- p6[block == "interaction"]
ix_death <- ix[outcome == "death" & method == "Cox" & level == "stg3"]
ix_stroke <- ix[outcome == "stroke" & method == "FineGray" & level == "stg3"]
exp_key <- c("VI-only" = "e_vi", "HI-only" = "e_hi", "DSI" = "e_ds")

f2_rows <- list()
for (oc in c("death", "stroke")) {
  f2_rows[[length(f2_rows) + 1]] <- list(panel_hdr = if (oc == "death") "A   All-cause death (Cox, M3)" else "B   Incident stroke (Fine-Gray, M3)")
  for (ex in EXORD) {
    f2_rows[[length(f2_rows) + 1]] <- list(grp = EXLBL[[ex]], sub = NA_character_,
      hr = NA_real_, lcl = NA_real_, ucl = NA_real_, s1 = NA_character_, s2 = NA_character_,
      s3 = NA_character_, s4 = NA_character_, pooled = TRUE)
    for (sg in c("ckm0_1", "ckm2", "ckm3")) {
      mth <- if (oc == "death") "Cox" else "FineGray"
      r <- stg[outcome == oc & exposure == ex & substratum == sg & method == mth]
      pooled <- nrow(r) > 0 && !is.na(r$pooled_HR[1])
      irow <- if (oc == "death") ix_death[exposure == exp_key[[ex]]] else ix_stroke[exposure == exp_key[[ex]]]
      pint <- if (sg == "ckm3" && nrow(irow) && !is.na(irow$p_pool[1])) fmt_p(irow$p_pool[1]) else "\u2014"
      if (pooled) {
        f2_rows[[length(f2_rows) + 1]] <- list(
          grp = NA_character_, sub = STGLBL[[sg]],
          hr = r$pooled_HR[1], lcl = r$pooled_lcl[1], ucl = r$pooled_ucl[1],
          s1 = fmt_hr(r$pooled_HR[1], r$pooled_lcl[1], r$pooled_ucl[1]), s2 = fmt_p(r$p[1]),
          s3 = paste0("k=", r$k[1]), s4 = pint, pooled = TRUE)
      } else {
        f2_rows[[length(f2_rows) + 1]] <- list(
          grp = NA_character_, sub = STGLBL[[sg]],
          hr = NA_real_, lcl = r$desc_min_HR[1], ucl = r$desc_max_HR[1],
          s1 = paste0(format(round(r$desc_min_HR[1], 2), nsmall = 2), "\u2013", format(round(r$desc_max_HR[1], 2), nsmall = 2)),
          s2 = "\u2014", s3 = paste0("k=", r$k[1]), s4 = pint, pooled = FALSE)
      }
    }
  }
}

leg2 <- c(
  "HRs for sensory impairment within CKM stages (ASN 2026 lower-bound staging, as in the main analysis).",
  "Solid squares = pooled estimates (REML, Knapp-Hartung; pooled when k ≥ 5, stratum events ≥ 10,",
  "exposed events ≥ 5).  Dashed lines = descriptive range of cohort-specific HRs",
  "when k < 5 (not pooled; not interpreted).  M3 covariates as in Figure 1.",
  "Stroke: Fine-Gray subdistribution hazard models, death as the competing event.",
  "P int = pooled multiplicative interaction of exposure with stage 3 (vs stages 0–1).",
  "Reference line at HR = 1; axis on the log scale.")
cs2 <- list(x = c(5.9, 7.95, 8.62, 9.42), heads = c("HR (95% CI)", "P", "k", "P int"),
            x_axis_left = 2.95, x_axis_right = 5.55, x_axis_mid = 4.25)
draw_forest(f2_rows, c(0.65, 2.6), c(0.75, 1, 1.5, 2, 2.5), cs2, leg2,
            file.path(D_OUT, "Figure2_论文版_分期森林图.pdf"), 6.3, "Hazard ratio (log scale)")

if (!inherits(logf, "try-error")) close(logf)
