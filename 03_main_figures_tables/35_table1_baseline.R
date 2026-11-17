# 35_table1_baseline.R - Table 1: baseline characteristics by four exposure classes (7 cohorts)
# Replicates STEP04 data assembly (20_step04_cohort_models.R section 1) to obtain the SAME main pool:
#   中介_界标表_v1.0 + 死亡时间_v2.3 + 卒中时间_v1.0, HRS half-sample, N = 45,949
# Exposure groups follow the model dummies exactly (e_vi/e_hi/e_ds; MHAS has no hearing item -> VI-only vs neither)
# Output: 04_结果/Table1_基线特征_四类暴露_v1.0.csv (long format) + Table1_逐库_四类暴露_v1.0.csv
suppressWarnings(suppressMessages({ library(data.table) }))
ROOT <- normalizePath(".", winslash = "/")
D_IN <- file.path(ROOT, "03_数据"); D_OUT <- file.path(ROOT, "04_结果")
logf <- file.path(ROOT, "05_质控", "Table1_run_log_20261110.txt")
lg <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n", file = logf, append = TRUE), silent = TRUE)
unlink(logf)

lm <- fread(file.path(D_IN, "中介_界标表_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(nrow(lm) == 50183)
dt <- fread(file.path(D_IN, "死亡时间_v2.3.csv"), colClasses = list(character = "id"))
dt <- unique(dt, by = c("cohort", "id"))
lm <- merge(lm, dt[, .(cohort, id, fu_d = fu_months, ev_d = event)], by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
sk <- fread(file.path(D_IN, "卒中时间_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(!any(duplicated(sk[, .(cohort, id)])))
lm <- merge(lm, sk[, .(cohort, id, fu_s = fu_s_months, ev_s = stroke_event, stroke_prev = stroke_prevalent)],
            by = c("cohort", "id"), all.x = TRUE, sort = FALSE)

lm[, main := (cohort != "HRS") | (hrs_half == 1)]
stopifnot(lm[main == TRUE, .N] == 45949)
stopifnot(lm[main == TRUE, sum(ev_d == 1, na.rm = TRUE)] == 9438)  # v2.3 CHARLS 2026-11-12 + HAALSI W4 2026-11-13
lg("main pool =", lm[main == TRUE, .N], "| deaths =", lm[main == TRUE, sum(ev_d == 1, na.rm = TRUE)])

# ---- exposure groups (identical to model dummies) ----
lm[, e_vi := fifelse(is.na(vi), NA_real_, as.numeric(vi == 1 & (is.na(hi) | hi == 0)))]
lm[, e_hi := fifelse(is.na(hi), NA_real_, as.numeric(hi == 1 & (is.na(vi) | vi == 0)))]
lm[, e_ds := fifelse(is.na(vi) | is.na(hi), NA_real_, as.numeric(vi == 1 & hi == 1))]
lm[cohort == "MHAS", `:=`(e_hi = NA_real_, e_ds = NA_real_)]
lm[, g4 := fcase(e_ds == 1, "DSI", e_hi == 1, "HI only", e_vi == 1, "VI only",
                 default = "Neither")]
M <- lm[main == TRUE]
lg("main pool exposure-missing (vi NA):", M[is.na(e_vi), .N],
   "| hi-missing with vi available (non-MHAS):", M[!is.na(e_vi) & is.na(e_hi) & cohort != "MHAS", .N],
   "| MHAS rows:", M[cohort == "MHAS", .N])
# model-complete exposure frame, identical to STEP04's exposure-class rule at model fitting:
# MHAS needs only e_vi defined (VI-only design); all other cohorts need all three dummies defined.
M <- M[!is.na(e_vi) & (cohort == "MHAS" | (!is.na(e_hi) & !is.na(e_ds)))]
lg("Table1 frame N =", nrow(M))
stopifnot(nrow(M) == 45565)                                   # 45,949 - 371 (vi NA) - 13 (hi missing, non-MHAS)
stopifnot(M[, sum(ev_d == 1, na.rm = TRUE)] == 9332)          # 9,438 - 105 - 1
stopifnot(M[stroke_prev == 0, sum(ev_s == 1, na.rm = TRUE)] == 2613)
lg("group counts:\n", paste(capture.output(print(M[, .N, by = g4])), collapse = "\n"))

M[, female := as.numeric(ragender == 2)]
M[, low_wealth := as.numeric(fifelse(cohort == "HAALSI", ses_wealth_q == 1, ses_wealth_q == 0))]
M[, stage := factor(ckm_stage_asn, levels = 0:3)]

GRP <- c("Neither", "VI only", "HI only", "DSI")
rows <- list()
add <- function(grp, char, level, stat, value, p = NA_real_, n_avail = NA_integer_) {
  rows[[length(rows) + 1]] <<- data.table(part = "A_pooled", group = grp, characteristic = char,
    level = level, stat = stat, value = value, p = p, n_avail = n_avail)
}
fmt_np <- function(k, n) sprintf("%d (%.1f)", k, 100 * k / n)
p_cat <- function(var, d) {
  tb <- table(d$g4, d[[var]]); if (any(dim(tb) < 2)) return(NA_real_)
  tryCatch(suppressWarnings(chisq.test(tb)$p.value), error = function(e) NA_real_)
}
p_num <- function(var, d) {
  tryCatch(summary(aov(d[[var]] ~ factor(d$g4)))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
}
add_num <- function(char, var, dec = 1) {
  d <- M[!is.na(get(var))]
  pv <- p_num(var, d)
  for (g in c(GRP, "Total")) {
    dd <- if (g == "Total") d else d[g4 == g]
    add(g, char, "", "mean (SD)", sprintf(paste0("%.", dec, "f (%.", dec, "f)"), mean(dd[[var]]), sd(dd[[var]])), pv, nrow(dd))
  }
}
add_bin <- function(char, var, label = "n (%)") {
  d <- M[!is.na(get(var))]
  pv <- p_cat(var, d)
  for (g in c(GRP, "Total")) {
    dd <- if (g == "Total") d else d[g4 == g]
    add(g, char, "", label, fmt_np(sum(dd[[var]] == 1, na.rm = TRUE), nrow(dd)), pv, nrow(dd))
  }
}
add_cat <- function(char, var, levels, labels = as.character(levels)) {
  d <- M[!is.na(get(var))]
  pv <- p_cat(var, d)
  for (i in seq_along(levels)) {
    for (g in c(GRP, "Total")) {
      dd <- if (g == "Total") d else d[g4 == g]
      add(g, char, labels[i], "n (%)", fmt_np(sum(dd[[var]] == levels[i], na.rm = TRUE), nrow(dd)), pv, nrow(dd))
    }
  }
}

# ---- N and demographics ----
for (g in c(GRP, "Total")) {
  dd <- if (g == "Total") M else M[g4 == g]
  add(g, "Participants, n", "", "n", as.character(nrow(dd)), NA_real_, nrow(dd))
}
add_num("Age, years", "_age", 1)
add_bin("Female", "female")
add_cat("Education (harmonized level)", "raeducel", c(1, 2, 3), c("< high school", "High school", "College or above"))
add_bin("Household income below cohort median", "ses_inc_median")
add_bin("Lowest wealth quartile/quintile", "low_wealth")
add_cat("CKM stage at baseline", "ckm_stage_asn", c(0, 1, 2, 3), c("Stage 0", "Stage 1", "Stage 2", "Stage 3"))
# ---- lifestyle / comorbidity ----
add_bin("Current smoker", "currsmoke")
add_bin("Current drinker", "current_drinker")
add_bin("Physically inactive", "pa_inactive")
add_bin("Cancer", "cancre")
add_bin("Chronic lung disease", "lunge")
add_bin("Arthritis", "arthre")
add_bin("Antihypertensive medication", "on_bp_med")
# ---- follow-up / events ----
for (g in c(GRP, "Total")) {
  dd <- if (g == "Total") M else M[g4 == g]
  f <- dd$fu_d[!is.na(dd$fu_d)]
  add(g, "Follow-up (death), months", "", "median (IQR)",
      sprintf("%.0f (%.0f\u2013%.0f)", median(f), quantile(f, .25), quantile(f, .75)), NA_real_, length(f))
}
for (g in c(GRP, "Total")) {
  dd <- if (g == "Total") M else M[g4 == g]
  add(g, "Deaths (main pool)", "", "n (%)", fmt_np(sum(dd$ev_d == 1, na.rm = TRUE), nrow(dd)), NA_real_, nrow(dd))
  add(g, "Incident stroke (first report)", "", "n (%)", fmt_np(sum(dd$ev_s == 1, na.rm = TRUE), nrow(dd)), NA_real_, nrow(dd))
}

T1 <- rbindlist(rows, fill = TRUE)
fwrite(T1, file.path(D_OUT, "Table1_基线特征_四类暴露_v1.0.csv"))

# ---- Part B: per cohort x exposure ----
B <- M[, .(n = .N, deaths = sum(ev_d == 1, na.rm = TRUE), strokes = sum(ev_s == 1, na.rm = TRUE),
           med_fu = as.numeric(median(fu_d, na.rm = TRUE))), by = .(cohort, g4)]
Bw <- dcast(B, cohort ~ g4, value.var = c("n", "deaths"))
B2 <- B[, .(n_total = sum(n), deaths = sum(deaths), strokes = sum(strokes),
            med_fu = as.numeric(median(M[cohort == .BY$cohort]$fu_d, na.rm = TRUE))), by = cohort]
setnames(B2, "cohort", "cohort")
Bw2 <- merge(B2, dcast(B, cohort ~ g4, value.var = "n"), by = "cohort", suffixes = c("", "_g"))
fwrite(Bw2, file.path(D_OUT, "Table1_逐库_四类暴露_v1.0.csv"))
lg("part B:\n", paste(capture.output(print(Bw2)), collapse = "\n"))

# ---- missingness for footnote ----
miss <- M[, .(n_total = .N,
  age = sum(is.na(`_age`)), female = sum(is.na(ragender)), edu = sum(is.na(raeducel)),
  income = sum(is.na(ses_inc_median)), wealth = sum(is.na(ses_wealth_q)),
  smoke = sum(is.na(currsmoke)), drink = sum(is.na(current_drinker)), pa = sum(is.na(pa_inactive)),
  cancer = sum(is.na(cancre)), lung = sum(is.na(lunge)), arthritis = sum(is.na(arthre)),
  bpmed = sum(is.na(on_bp_med)), stage = sum(is.na(ckm_stage_asn)))]
fwrite(miss, file.path(D_OUT, "Table1_缺失率_v1.0.csv"))

# ---- stage distribution by cohort (main pool, exposure-complete) ----
S4 <- M[, .N, by = .(cohort, ckm_stage_asn)][order(cohort, ckm_stage_asn)]
S4w <- dcast(S4, cohort ~ ckm_stage_asn, value.var = "N", fill = 0)
S4w[, n_total := rowSums(.SD), .SDcols = as.character(0:3)]
fwrite(S4w, file.path(D_OUT, "TableS4_分期分布_逐库_v1.0.csv"))
lg("S4 stage by cohort:
", paste(capture.output(print(S4w)), collapse = "
"))
lg("missingness:\n", paste(capture.output(print(miss)), collapse = "\n"))
lg("done. rows =", nrow(T1))
