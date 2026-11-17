# =====================================================================
# 20_step04_cohort_models.R  (v1.0, 2026-11-06)
# Topic 6 (Sensory impairment x CKM) - STEP04 in-cohort models
#
# Frozen design: v1.0 co-author doc + v0.3 SAP; deviations registered in
#   01_方案设计/变更登记_20261106_STEP04前置口径.md (RC-window = person exclusion,
#   SES = raeducel 3-level + within-cohort wealth quantiles) - both parameterized.
#
# Exposure : ONE model with three mutually-exclusive dummies vs "neither":
#            e_vi (VI-only) / e_hi (HI-only) / e_ds (DSI).
#            MHAS: e_vi only (no general hearing item).
#            HAALSI: HI binary construct -> rows kept, note = hi_binary.
# Outcomes : death  -> Cox (main, unweighted, calendar months)
#            stroke -> Fine-Gray (main; death competing) + cause-specific Cox (secondary)
# Models   : M1 age_y + ragender; M2 + raeducel + ses_inc_median + ses_wealth_q;
#            M3 = M2 + currsmoke + current_drinker + pa_inactive + cancre + lunge
#                 + arthre + on_bp_med
#            HAALSI (reduced): M2 = M1 + raeducel + ses_wealth_q (no income);
#            M3 = M2 + currsmoke + current_drinker + pa_inactive + hearte + has_dm
#            HRS: wave (8/9) forced into every model; hrs_half == 1 only.
# MI       : m=20 per cohort, covariates only; predictors = covariates + exposure
#            dummies + death outcome + Nelson-Aalen cumhaz + stroke indicator.
#            MI = main caliber; complete-case (cc) rows in same output for S6.
# PH       : cox.zph on complete-case M3 (death) and M3 cause-specific (stroke).
# Strata   : CKM stage 0-1 / 2 / 3 under M3 (MI). HR reported only if stratum
#            events >= 10 AND exposure-cell events >= 5 (frozen sparse rules).
# Dose-res : ordinal per-1 HR + RCS (3 knots, splines::ns) on CC, death endpoint.
#            Hearing skipped where absent/binary; vision = far, else near/overall.
# RC window: persons whose OWN-endpoint event occurs < RC_WINDOW months are removed
#            (death: early deaths dropped; stroke: early strokes dropped, early
#            deaths kept as competing events).
# Seed     : 20260904 (+ per-cohort offset for mice chains)
# ASCII-only source (GBK console); all outputs UTF-8-safe CSVs.
# =====================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(survival); library(cmprsk); library(mice); library(splines)
}))

SEED <- 20260904
M_IMP <- 20
RC_WINDOW <- 24
SPARSE_STRATUM <- 10
SPARSE_CELL <- 5
set.seed(SEED)

ROOT <- normalizePath(".", winslash = "/")
D_IN <- file.path(ROOT, "03_数据"); D_OUT <- file.path(ROOT, "04_结果"); D_QC <- file.path(ROOT, "05_质控")
LOGF <- file.path(D_QC, "STEP04_run_log_20261106.txt")
logf <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "
", file = LOGF, append = TRUE), silent = TRUE)
unlink(LOGF)
cat("STEP04 log ->", LOGF, "\n")

MASTER <- "D:/数据/国际老年数据库/00_跨库资源/10_主题研究/06_生物标志物Biomarker/002_10年心血管风险计算&CKM/03_数据/ALL_DB_CKM_Master_v1.6.csv"
COV46  <- "D:/数据/国际老年数据库/00_跨库资源/10_主题研究/03_通用协变量CommonCovariates/03_输出数据/ALL_DB_Covariates_Long_v4.6.csv"
logf("R:", R.version.string)
logf("survival", as.character(packageVersion("survival")), "| cmprsk",
     as.character(packageVersion("cmprsk")), "| mice", as.character(packageVersion("mice")),
     "| data.table", as.character(packageVersion("data.table")))
logf("params: M_IMP=", M_IMP, " RC_WINDOW=", RC_WINDOW, " SPARSE=", SPARSE_STRATUM, "/", SPARSE_CELL, " SEED=", SEED)

norm_id <- function(x) sub("^0+", "", sub("\\.0$", "", tolower(as.character(x))))

# ---------------------------------------------------------------
# 1. Data assembly
# ---------------------------------------------------------------
lm <- fread(file.path(D_IN, "中介_界标表_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(nrow(lm) == 50183)

# death time v2.1 (2026-11-06: CHARLS wave-midpoint correction, see 21_rebuild_outcome_times.py)
dt <- fread(file.path(D_IN, "死亡时间_v2.3.csv"), colClasses = list(character = "id"))
# duplicated keys only differ in base_wave (verified 2026-11-06); fu/event are identical
stopifnot(dt[, max(uniqueN(paste(fu_months, event))), by = .(cohort, id)]$V1 == 1)
dt <- unique(dt, by = c("cohort", "id"))
stopifnot(!any(duplicated(dt[, .(cohort, id)])))
lm <- merge(lm, dt[, .(cohort, id, fu_d = fu_months, ev_d = event)],
            by = c("cohort", "id"), all.x = TRUE, sort = FALSE)

# stroke time v1.0 (2026-11-06: interval-censored rebuild per frozen rule; HAALSI passthrough)
sk <- fread(file.path(D_IN, "卒中时间_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(!any(duplicated(sk[, .(cohort, id)])))
stopifnot(nrow(sk) == 50183)
lm <- merge(lm, sk[, .(cohort, id, fu_s = fu_s_months, ev_s = stroke_event,
                       stroke_prev = stroke_prevalent)],
            by = c("cohort", "id"), all.x = TRUE, sort = FALSE)

BASE_WAVE <- c(CHARLS = 1, ELSA = 2, SHARE = 6, KLoSA = 1, MHAS = 3, HAALSI = 1)
mv <- fread(MASTER, select = c("id", "study", "wave", "baseline_prevalent_hearte", "has_dm"),
            colClasses = list(character = "id"))
mv <- mv[study != "HRS_VBS"]
mv[, base := (study == "HRS_DBS" & wave %in% c(8, 9)) |
            (study != "HRS_DBS" & wave == BASE_WAVE[study])]
mv <- mv[base == TRUE][, cohort := fifelse(study == "HRS_DBS", "HRS", study)][, id := norm_id(id)]
stopifnot(!any(duplicated(mv[, .(cohort, id)])))
mv <- mv[, .(cohort, id, heart_prev = baseline_prevalent_hearte, has_dm)]
lm <- merge(lm, mv, by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
stopifnot(sum(is.na(lm$ev_s)) == 0)

cv <- fread(COV46, select = c("Sourcedatabase", "pid_unique", "wave", "hearte", "diabe"),
            colClasses = list(character = "pid_unique"))
cv[, cohort := Sourcedatabase][, id := norm_id(sub("^[^_]+_", "", pid_unique))]
cv <- cv[cohort %in% c("CHARLS", "ELSA", "SHARE", "KLoSA", "MHAS", "HAALSI") & wave == BASE_WAVE[cohort]]
cv <- unique(cv[, .(cohort, id, hearte46 = hearte)])
stopifnot(!any(duplicated(cv[, .(cohort, id)])))
lm <- merge(lm, cv, by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
# heart disease source for HAALSI reduced M3: Master prevalent flag, fallback v4.6 hearte
lm[, heart_d := fifelse(cohort == "HAALSI",
                        fifelse(!is.na(heart_prev), heart_prev, hearte46), NA_real_)]
ha_cov <- lm[cohort == "HAALSI", .(heart_cov = mean(!is.na(heart_d)), dm_cov = mean(!is.na(has_dm)))]
logf("[join] HAALSI hearte cov=", round(ha_cov$heart_cov, 4), " has_dm cov=", round(ha_cov$dm_cov, 4))
stopifnot(ha_cov$dm_cov > 0.99)

lm[, main := (cohort != "HRS") | (hrs_half == 1)]
lm[cohort == "HRS", wave := base_wave]   # landmark table carries base_wave (8/9)
stopifnot(lm[main == TRUE, .N] == 45949)
stopifnot(lm[main == TRUE, sum(ev_d == 1, na.rm = TRUE)] == 9438)  # v2.3 CHARLS 2026-11-12 + HAALSI W4 2026-11-13
logf("[1] universe ok: 50,183 | main pool 45,949 | deaths(main) = 9,438")

# ---------------------------------------------------------------
# 2. Exposure dummies (checked vs exp4)
# ---------------------------------------------------------------
lm[, e_vi := fifelse(is.na(vi), NA_real_, as.numeric(vi == 1 & (is.na(hi) | hi == 0)))]
lm[, e_hi := fifelse(is.na(hi), NA_real_, as.numeric(hi == 1 & (is.na(vi) | vi == 0)))]
lm[, e_ds := fifelse(is.na(vi) | is.na(hi), NA_real_, as.numeric(vi == 1 & hi == 1))]
lm[cohort == "MHAS", `:=`(e_hi = NA_real_, e_ds = NA_real_)]   # no hearing item
lm[, exp4n := suppressWarnings(as.numeric(exp4))]
chk <- lm[!is.na(exp4n)]
mism <- sum((chk$e_vi == 1 & chk$exp4n != 1) | (chk$e_hi == 1 & chk$exp4n != 2) |
            (chk$e_ds == 1 & chk$exp4n != 3) |
            (chk$e_vi == 0 & chk$e_hi == 0 & chk$e_ds == 0 & chk$exp4n != 0), na.rm = TRUE)
logf("[2] exposure vs exp4 mismatches = ", mism)
stopifnot(mism == 0)

# ---------------------------------------------------------------
# 3. Model specs
# ---------------------------------------------------------------
M1V <- c("age_y", "ragender")
M2V_FULL <- c(M1V, "raeducel", "ses_inc_median", "ses_wealth_q")
M3V_FULL <- c(M2V_FULL, "currsmoke", "current_drinker", "pa_inactive", "cancre", "lunge", "arthre", "on_bp_med")
M2V_HAALSI <- c(M1V, "raeducel", "ses_wealth_q")
M3V_HAALSI <- c(M2V_HAALSI, "currsmoke", "current_drinker", "pa_inactive", "heart_d", "has_dm")
FACTORV <- c("ragender", "raeducel", "ses_wealth_q", "wave")

get_vars <- function(co, model) {
  v <- switch(model, M1 = M1V,
              M2 = if (co == "HAALSI") M2V_HAALSI else M2V_FULL,
              M3 = if (co == "HAALSI") M3V_HAALSI else M3V_FULL)
  if (co == "HRS") v <- c(v, "wave")
  v
}
model_set_of <- function(co, model) if (co == "HAALSI" && model == "M3") "reduced" else "full"
exp_terms <- function(co) if (co == "MHAS") "e_vi" else "e_vi + e_hi + e_ds"
exp_dums  <- function(co) if (co == "MHAS") c("e_vi") else c("e_vi", "e_hi", "e_ds")

fact <- function(d, vars) { for (v in intersect(FACTORV, vars)) d[[v]] <- droplevels(factor(d[[v]])); d }
cc_keep <- function(d, vars) Reduce(`&`, lapply(vars, function(v) !is.na(d[[v]])))

pool_rubin <- function(betas, vcovs, dfcom, tag = "") {
  m0 <- length(betas)
  # per-imp: drop coefficients that are NA/Inf (rank-deficient terms, e.g. constant
  # covariates in strata); imps that lose a term only occasionally are dropped entirely
  clean <- lapply(seq_len(m0), function(i) {
    b <- betas[[i]]; v <- vcovs[[i]]
    okn <- is.finite(b) & is.finite(diag(v))
    list(b = b[okn], v = v[okn, okn, drop = FALSE],
         dropped = names(b)[!okn])
  })
  keep <- vapply(clean, function(x) length(x$b) > 0 && nrow(x$v) == length(x$b), TRUE)
  if (sum(keep) < m0)
    logf("  WARN pool [", tag, "]: drop ", m0 - sum(keep), "/", m0,
         " imp(s) with unusable fits")
  if (sum(keep) < max(5, ceiling(m0 / 2)))
    stop("pool_rubin [", tag, "]: too few usable imputations (", sum(keep), "/", m0, ")")
  clean <- clean[keep]
  m <- length(clean)
  all_drop <- unique(unlist(lapply(clean, `[[`, "dropped")))
  if (length(all_drop))
    logf("  NOTE pool [", tag, "]: degenerate term(s) dropped from pooling: ",
         paste(all_drop, collapse = ", "))
  betas <- lapply(clean, `[[`, "b"); vcovs <- lapply(clean, `[[`, "v")
  nm0 <- names(betas[[1]])
  nms <- Reduce(intersect, lapply(betas, names))
  if (length(nms) < length(nm0)) logf("  WARN: pooling drops ",
      length(nm0) - length(nms), " coefficient(s) not present in all imps: ",
      paste(setdiff(nm0, nms), collapse = ", "))
  b  <- do.call(cbind, lapply(betas, function(x) x[nms]))
  vs <- lapply(vcovs, function(v) v[nms, nms, drop = FALSE])
  stopifnot(all(is.finite(b)), all(vapply(vs, function(v) all(is.finite(v)), TRUE)))
  bbar <- rowMeans(b); names(bbar) <- nms
  Ubar <- Reduce(`+`, vs) / m
  Bm <- cov(t(b)); Tm <- Ubar + (1 + 1/m) * Bm
  lam <- pmin(pmax((1 + 1/m) * diag(Bm) / diag(Tm), 0), 0.999999)
  dfold <- (m - 1) / lam^2
  dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  df <- 1 / (1/dfold + 1/dfobs)
  se <- sqrt(diag(Tm)); names(se) <- nms
  z <- bbar / se
  pv <- 2 * pt(-abs(z), df); names(pv) <- nms
  list(beta = bbar, se = se, p = pv)
}
row_for <- function(e, co, oc, me, mo, cal, st, ms, n, ne, hr, lo, hi, pp, note) data.table(
  cohort = co, exposure = e, outcome = oc, method = me, model = mo, caliber = cal,
  substratum = st, model_set = ms, n = n, n_events = ne, HR = hr, lcl = lo, ucl = hi,
  p = pp, ph_p = NA_real_, rc_window = RC_WINDOW, note = note)

fit_fg <- function(d, vars, tag = "") {
  d <- as.data.table(d)
  cf <- fact(as.data.frame(d[, .SD, .SDcols = vars]), vars)
  ok <- cc_keep(cf, vars); cf <- cf[ok, ]
  X <- model.matrix(~ ., cf)[, -1, drop = FALSE]
  if (ncol(X) > 0) {
    sdv <- suppressWarnings(apply(X, 2, sd))
    cst <- !is.finite(sdv) | sdv == 0
    if (any(cst)) {
      logf("  NOTE FG [", tag, "]: constant design column(s) dropped: ",
           paste(colnames(X)[cst], collapse = ", "))
      X <- X[, !cst, drop = FALSE]
    }
  }
  dd <- d[ok, ]
  f <- tryCatch(crr(dd$t, dd$fstatus, X, failcode = 1, cencode = 0), error = function(e) e)
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f)))
  dimnames(f$var) <- list(names(f$coef), names(f$coef))
  list(ok = TRUE, coef = f$coef, var = f$var, n = nrow(X), nevent = sum(dd$fstatus == 1), keep = ok)
}
fit_coxf <- function(form, d) {
  warns <- character(0)
  f <- tryCatch(withCallingHandlers(
         coxph(form, data = d, ties = "efron", model = TRUE),
         warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
         error = function(e) e)
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f)))
  list(ok = TRUE, fit = f, warn = warns)
}
lrt_p <- function(f0, f1) {   # LR test p, nested cox models on identical data
  ch <- 2 * (f1$loglik[2] - f0$loglik[2])
  dfd <- length(coef(f1)) - length(coef(f0))
  if (is.na(ch) || ch < 0 || dfd <= 0) return(NA_real_)
  pchisq(ch, df = dfd, lower.tail = FALSE)
}

# ---------------------------------------------------------------
# 4. Per-cohort driver
# ---------------------------------------------------------------
cohs <- c("CHARLS", "ELSA", "SHARE", "KLoSA", "MHAS", "HAALSI", "HRS")
res_all <- list(); qc_rows <- list(); ph_rows <- list(); dose_rows <- list()

for (co in cohs) {
  t0 <- Sys.time(); logf("=== ", co, " ===")
  uni <- lm[cohort == co & main == TRUE]

  # --- endpoint frames ---
  dset <- uni[!is.na(fu_d) & !is.na(ev_d) & !is.na(age_y)]
  drop_d <- dset[ev_d == 1 & fu_d < RC_WINDOW, .N]
  dset <- dset[!(ev_d == 1 & fu_d < RC_WINDOW)]
  dset[, `:=`(t = fu_d, ev = as.integer(ev_d == 1))]
  sset0 <- uni[!is.na(fu_s) & !is.na(age_y) & !is.na(ev_s) & stroke_prev == 0]
  drop_s <- sset0[ev_s == 1 & fu_s < RC_WINDOW, .N]
  sset <- sset0[!(ev_s == 1 & fu_s < RC_WINDOW)]
  sset[, `:=`(t = fu_s, ev = as.integer(ev_s == 1),
              fstatus = fifelse(ev_s == 1, 1L, fifelse(!is.na(ev_d) & ev_d == 1, 2L, 0L)))]
  logf("  death n=", nrow(dset), " ev=", sum(dset$ev), " (early dropped ", drop_d, ") | stroke n=",
       nrow(sset), " ev=", sum(sset$ev), " comp=", sum(sset$fstatus == 2), " (early dropped ", drop_s, ")")

  qc_rows[[co]] <- data.table(
    cohort = co, n_pool = uni[, .N],
    n_death = nrow(dset), ev_death_used = sum(dset$ev), ev_death_dropped_early = drop_d,
    n_stroke = nrow(sset), ev_stroke_used = sum(sset$ev), ev_stroke_dropped_early = drop_s,
    n_stroke_baseline_prev = uni[, sum(stroke_prev == 1, na.rm = TRUE)],
    competing_death_in_stroke_set = sum(sset$fstatus == 2),
    med_fu_death = round(median(dset$t), 1), med_fu_stroke = round(median(sset$t), 1),
    n_vi_missing = uni[, sum(is.na(vi))],
    n_hi_missing = fifelse(co == "MHAS", NA_integer_, uni[, sum(is.na(hi))]))

  # --- MI (union frame = uni rows) ---
  m3v <- get_vars(co, "M3"); dums <- exp_dums(co)
  mi_f <- as.data.frame(uni[, .SD, .SDcols = unique(c("id", m3v, dums))])
  mi_f <- fact(mi_f, m3v)
  ph_t <- uni$fu_d; ph_e <- uni$ev_d
  ph_t[is.na(ph_t)] <- median(ph_t, na.rm = TRUE); ph_e[is.na(ph_e)] <- 0
  mi_f$t <- ph_t; mi_f$ev <- ph_e
  mi_f$cumhaz <- mice::nelsonaalen(mi_f, t, ev)
  mi_f$ev_s_pred <- fifelse(is.na(uni$ev_s), 0, uni$ev_s)
  imp_vars <- names(mi_f)[vapply(mi_f, function(x) anyNA(x) && !all(is.na(x)), TRUE)]
  imp_vars <- setdiff(imp_vars, c("id", "t", "ev", "cumhaz", "ev_s_pred"))
  imp_vars <- setdiff(imp_vars, dums)   # exposure is NEVER imputed (frozen rule)
  logf("  MI imputes: ", paste(imp_vars, collapse = ", "))
  set.seed(SEED + match(co, cohs) * 1000)
  imp <- mice(mi_f[, setdiff(names(mi_f), "id"), drop = FALSE], m = M_IMP, maxit = 5, printFlag = FALSE)
  comps <- lapply(seq_len(M_IMP), function(i) {
    ci <- complete(imp, i); ci$id <- mi_f$id; ci })
  logf("  MI done")

  # helper: imputed copy of a frame
  imp_copy <- function(fr, i) {
    ci <- comps[[i]]
    d2 <- copy(fr)
    for (v in imp_vars) d2[[v]] <- ci[[v]][match(d2$id, ci$id)]
    fact(d2, m3v)
  }

  # --- MI grid: M1/M2/M3 x (death Cox, stroke FG, stroke CS) ---
  for (mod in c("M1", "M2", "M3")) {
    vars <- get_vars(co, mod); mset <- model_set_of(co, mod)
    # death
    fits <- lapply(seq_len(M_IMP), function(i) {
      dd <- imp_copy(dset, i)
      dd <- dd[cc_keep(dd, c(dums, vars)), ]
      fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), dd)
    })
    okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
    if (sum(okf) > 0) {
      na1 <- names(which(!is.finite(stats::coef(fits[[which(okf)[1]]]$fit))))
      if (length(na1)) logf("  DEBUG [", co, "|", mod, "|death|mi] imp", which(okf)[1],
                            " NA coefs: ", paste(na1, collapse = ", "))
    } else logf("  DEBUG [", co, "|", mod, "|death|mi] all imps errored: ",
                utils::head(fits[[1]]$msg, 1))
    if (sum(okf) == M_IMP) {
      p <- length(stats::coef(fits[[1]]$fit))
      nn <- sum(cc_keep(as.data.frame(imp_copy(dset, 1)[, .SD, .SDcols = c(dums, vars)]), c(dums, vars)))
      ne <- sum(dset$ev == 1)
      pf <- pool_rubin(lapply(fits, function(x) stats::coef(x$fit)),
                       lapply(fits, function(x) vcov(x$fit)), dfcom = nn - p,
                       tag = paste0(co, "|", mod, "|death|Cox|mi|overall"))
      for (e in dums) {
        res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", mod, "mi", "overall",
          mset, nn, ne, exp(pf$beta[e]), exp(pf$beta[e] - 1.96 * pf$se[e]),
          exp(pf$beta[e] + 1.96 * pf$se[e]), pf$p[e], NA_character_)
      }
    } else logf("  WARN: MI death ", mod, " failed in ", M_IMP - sum(okf), " imps")
    # stroke FG
    fits <- lapply(seq_len(M_IMP), function(i) fit_fg(imp_copy(sset, i), c(dums, vars)))
    okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
    if (sum(okf) > 0) {
      na1 <- names(which(!is.finite(fits[[which(okf)[1]]]$coef)))
      if (length(na1)) logf("  DEBUG [", co, "|", mod, "|FG|mi] imp", which(okf)[1],
                            " NA coefs: ", paste(na1, collapse = ", "))
    } else logf("  DEBUG [", co, "|", mod, "|FG|mi] all imps errored: ",
                utils::head(fits[[1]]$msg, 1))
    if (sum(okf) == M_IMP) {
      p <- length(fits[[1]]$coef); nn <- fits[[1]]$n; ne <- fits[[1]]$nevent
      pf <- pool_rubin(lapply(fits, `[[`, "coef"), lapply(fits, `[[`, "var"), dfcom = ne - p,
                       tag = paste0(co, "|", mod, "|stroke|FG|mi|overall"))
      nm <- names(fits[[1]]$coef)
      for (e in dums) {
        j <- which(nm == e); if (length(j) != 1) next
        res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", mod, "mi", "overall",
          mset, nn, ne, exp(pf$beta[j]), exp(pf$beta[j] - 1.96 * pf$se[j]),
          exp(pf$beta[j] + 1.96 * pf$se[j]), pf$p[j], NA_character_)
      }
    } else logf("  WARN: MI FG ", mod, " failed in ", M_IMP - sum(okf), " imps")
    # stroke CS Cox
    fits <- lapply(seq_len(M_IMP), function(i) {
      ss <- imp_copy(sset, i); ss <- ss[cc_keep(ss, c(dums, vars)), ]
      fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), ss)
    })
    okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
    if (sum(okf) == M_IMP) {
      p <- length(stats::coef(fits[[1]]$fit))
      nn <- sum(cc_keep(as.data.frame(imp_copy(sset, 1)[, .SD, .SDcols = c(dums, vars)]), c(dums, vars)))
      ne <- sum(sset$ev == 1)
      pf <- pool_rubin(lapply(fits, function(x) stats::coef(x$fit)),
                       lapply(fits, function(x) vcov(x$fit)), dfcom = nn - p,
                       tag = paste0(co, "|", mod, "|stroke|CS|mi|overall"))
      for (e in dums) {
        res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "CauseSpecificCox", mod, "mi", "overall",
          mset, nn, ne, exp(pf$beta[e]), exp(pf$beta[e] - 1.96 * pf$se[e]),
          exp(pf$beta[e] + 1.96 * pf$se[e]), pf$p[e], NA_character_)
      }
    } else logf("  WARN: MI CS ", mod, " failed in ", M_IMP - sum(okf), " imps")
  }

  # --- CC rows (M1/M3 death; M3 FG; M3 CS) ---
  for (mod in c("M1", "M3")) {
    vars <- get_vars(co, mod); mset <- model_set_of(co, mod)
    dd <- fact(as.data.frame(dset[, .SD, .SDcols = c("t", "ev", "id", dums, vars)]), vars)
    dd <- dd[cc_keep(dd, c(dums, vars)), ]
    fc <- fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), dd)
    if (isTRUE(fc$ok)) {
      s <- summary(fc$fit)$coefficients; ci <- summary(fc$fit)$conf.int
      for (e in dums) res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", mod, "cc",
        "overall", mset, nrow(dd), sum(dd$ev), ci[e, "exp(coef)"], ci[e, "lower .95"],
        ci[e, "upper .95"], s[e, "Pr(>|z|)"], NA_character_)
    }
    if (mod == "M3") {
      sf <- fact(as.data.frame(sset[, .SD, .SDcols = c("t", "fstatus", dums, vars)]), vars)
      sf <- sf[cc_keep(sf, c(dums, vars)), ]
      ff <- fit_fg(sf, c(dums, vars))
      if (isTRUE(ff$ok)) {
        se <- sqrt(diag(ff$var)); z <- ff$coef / se
        for (e in dums) {
          j <- which(names(ff$coef) == e); if (length(j) != 1) next
          res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", mod, "cc",
            "overall", mset, ff$n, ff$nevent, exp(ff$coef[j]), exp(ff$coef[j] - 1.96 * se[j]),
            exp(ff$coef[j] + 1.96 * se[j]), 2 * pnorm(-abs(z[j])), NA_character_)
        }
      } else logf("  WARN: CC FG M3 failed: ", ff$msg)
      ss <- fact(as.data.frame(sset[, .SD, .SDcols = c("t", "ev", dums, vars)]), vars)
      ss <- ss[cc_keep(ss, c(dums, vars)), ]
      fc2 <- fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), ss)
      if (isTRUE(fc2$ok)) {
        s <- summary(fc2$fit)$coefficients; ci <- summary(fc2$fit)$conf.int
        for (e in dums) res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "CauseSpecificCox",
          mod, "cc", "overall", mset, nrow(ss), sum(ss$ev), ci[e, "exp(coef)"], ci[e, "lower .95"],
          ci[e, "upper .95"], s[e, "Pr(>|z|)"], NA_character_)
        # PH
        zp <- tryCatch(cox.zph(fc2$fit), error = function(e) e)
        if (!inherits(zp, "error")) {
          tt <- as.data.frame(zp$table)
          ph_rows[[paste0(co, "_stroke")]] <- data.table(cohort = co, outcome = "stroke",
            term = rownames(tt), chisq = tt$chisq, df = tt$df, p = tt$p,
            n_cc = nrow(ss), n_events_cc = sum(ss$ev))
        } else logf("  NOTE: cox.zph stroke failed for ", co, ": ",
                    conditionMessage(zp))
      }
    }
  }
  # PH death (CC M3)
  vars <- get_vars(co, "M3")
  dd <- fact(as.data.frame(dset[, .SD, .SDcols = c("t", "ev", dums, vars)]), vars)
  dd <- dd[cc_keep(dd, c(dums, vars)), ]
  fc <- fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), dd)
  if (isTRUE(fc$ok)) {
    zp <- tryCatch(cox.zph(fc$fit), error = function(e) e)
    if (!inherits(zp, "error")) {
      tt <- as.data.frame(zp$table)
      ph_rows[[paste0(co, "_death")]] <- data.table(cohort = co, outcome = "death",
        term = rownames(tt), chisq = tt$chisq, df = tt$df, p = tt$p,
        n_cc = nrow(dd), n_events_cc = sum(dd$ev))
    } else logf("  NOTE: cox.zph death failed for ", co, ": ", conditionMessage(zp))
  }

  # --- CKM strata (M3, MI) ---
  dset[, stg := fifelse(ckm_stage_asn %in% c(0, 1), "ckm0_1", fifelse(ckm_stage_asn == 2, "ckm2", "ckm3"))]
  sset[, stg := fifelse(ckm_stage_asn %in% c(0, 1), "ckm0_1", fifelse(ckm_stage_asn == 2, "ckm2", "ckm3"))]
  vars <- get_vars(co, "M3"); mset <- model_set_of(co, "M3")
  for (stg in c("ckm0_1", "ckm2", "ckm3")) {
    stg_v <- stg
    ds <- dset[stg == stg_v]; ss <- sset[stg == stg_v]
    # death
    fits <- lapply(seq_len(M_IMP), function(i) {
      dd <- imp_copy(ds, i); dd <- dd[cc_keep(dd, c(dums, vars)), ]
      if (sum(dd$ev) < SPARSE_STRATUM) return(NULL)
      fit_coxf(as.formula(paste("Surv(t, ev) ~", exp_terms(co), "+", paste(vars, collapse = " + "))), dd)
    })
    okf <- vapply(fits, function(x) !is.null(x) && isTRUE(x$ok), TRUE)
    for (e in dums) {
      ne_cell <- sum(ds[get(e) == 1, ev], na.rm = TRUE)
      ne_model <- sum(ds$ev, na.rm = TRUE); nn <- sum(!is.na(ds[[e]]))
      if (sum(okf) == M_IMP) {
        if (ne_cell >= SPARSE_CELL) {
          pf <- tryCatch({
            p <- length(stats::coef(fits[[1]]$fit))
            pool_rubin(lapply(fits, function(x) stats::coef(x$fit)),
                       lapply(fits, function(x) vcov(x$fit)), dfcom = nn - p,
                       tag = paste0(co, "|M3|death|Cox|mi|", stg))
          }, error = function(e) e)
          if (inherits(pf, "error")) {
            res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", "M3", "mi", stg,
              mset, nn, ne_model, NA_real_, NA_real_, NA_real_, NA_real_, "degenerate_fit")
          } else res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", "M3", "mi", stg,
            mset, nn, ne_model, exp(pf$beta[e]), exp(pf$beta[e] - 1.96 * pf$se[e]),
            exp(pf$beta[e] + 1.96 * pf$se[e]), pf$p[e], NA_character_)
        } else res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", "M3", "mi", stg,
          mset, nn, ne_cell, NA_real_, NA_real_, NA_real_, NA_real_, "sparse_cell_lt5")
      } else res_all[[length(res_all) + 1]] <- row_for(e, co, "death", "Cox", "M3", "mi", stg,
        mset, nn, ne_model, NA_real_, NA_real_, NA_real_, NA_real_,
        if (ne_model < SPARSE_STRATUM) "sparse_stratum_lt10" else "fit_failed")
    }
    # stroke FG
    fits <- lapply(seq_len(M_IMP), function(i) {
      s2 <- imp_copy(ss, i)
      if (sum(s2$fstatus == 1) < SPARSE_STRATUM) return(NULL)
      fit_fg(s2, c(dums, vars))
    })
    okf <- vapply(fits, function(x) !is.null(x) && isTRUE(x$ok), TRUE)
    for (e in dums) {
      ne_cell <- sum(ss[get(e) == 1, ev], na.rm = TRUE)
      ne_model <- sum(ss$ev, na.rm = TRUE); nn <- sum(!is.na(ss[[e]]))
      if (sum(okf) == M_IMP) {
        if (ne_cell >= SPARSE_CELL) {
          pf <- tryCatch({
            p <- length(fits[[1]]$coef)
            pool_rubin(lapply(fits, `[[`, "coef"), lapply(fits, `[[`, "var"), dfcom = ne_cell - p,
                       tag = paste0(co, "|M3|stroke|FG|mi|", stg))
          }, error = function(e) e)
          j <- which(names(fits[[1]]$coef) == e)
          if (inherits(pf, "error")) {
            res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", "M3", "mi", stg,
              mset, nn, ne_model, NA_real_, NA_real_, NA_real_, NA_real_, "degenerate_fit")
          } else res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", "M3", "mi", stg,
            mset, nn, ne_model, exp(pf$beta[j]), exp(pf$beta[j] - 1.96 * pf$se[j]),
            exp(pf$beta[j] + 1.96 * pf$se[j]), pf$p[j], NA_character_)
        } else res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", "M3", "mi", stg,
          mset, nn, ne_cell, NA_real_, NA_real_, NA_real_, NA_real_, "sparse_cell_lt5")
      } else res_all[[length(res_all) + 1]] <- row_for(e, co, "stroke", "FineGray", "M3", "mi", stg,
        mset, nn, ne_model, NA_real_, NA_real_, NA_real_, NA_real_,
        if (ne_model < SPARSE_STRATUM) "sparse_stratum_lt10" else "fit_failed")
    }
  }
  dset[, stg := NULL]; sset[, stg := NULL]

  # --- dose-response (CC, death) ---
  dose_vars <- c("hearing", "vision_far", "vision_near", "vision_overall")
  done_v <- character(0)
  for (v0 in c("hearing", "vision")) {
    cand <- if (v0 == "hearing") "hearing" else c("vision_far", "vision_near", "vision_overall")
    if (v0 == "hearing" && co == "MHAS") next
    for (vv in cand) {
      if (!vv %in% names(dset)) next
      vals <- unique(na.omit(dset[[vv]])); vals <- vals[vals %in% 1:5]
      if (length(vals) < 3) next
      v <- vv; done_v <- c(done_v, vv)
      for (mod in c("M1", "M3")) {
        vars <- get_vars(co, mod)
        dd <- as.data.frame(dset[, .SD, .SDcols = c("t", "ev", vv, vars)])
        dd <- fact(dd, vars)
        dd[[vv]][!dd[[vv]] %in% 1:5] <- NA
        dd <- dd[cc_keep(dd, c(vv, vars)), ]
        if (length(unique(dd[[vv]])) < 3) next
        fc <- fit_coxf(as.formula(paste("Surv(t, ev) ~", vv, "+", paste(vars, collapse = " + "))), dd)
        kn <- quantile(dd[[vv]], c(.10, .50, .90)); pn <- NA_real_
        if (length(unique(kn)) == 3) {
          fc2 <- fit_coxf(as.formula(paste("Surv(t, ev) ~ ns(", vv, ", knots = c(",
                                           paste(kn, collapse = ","), ")) +", paste(vars, collapse = " + "))), dd)
          if (isTRUE(fc$ok) && isTRUE(fc2$ok)) {
            pn <- lrt_p(fc$fit, fc2$fit)
          }
        }
        if (isTRUE(fc$ok)) {
          s <- summary(fc$fit)$coefficients; ci <- summary(fc$fit)$conf.int
          dose_rows[[length(dose_rows) + 1]] <- data.table(cohort = co, dose_var = vv, model = mod,
            n = nrow(dd), n_events = sum(dd$ev), per1_HR = ci[vv, "exp(coef)"],
            per1_lcl = ci[vv, "lower .95"], per1_ucl = ci[vv, "upper .95"],
            per1_p = s[vv, "Pr(>|z|)"], p_nonlinear = pn,
            note = if (v0 == "hearing" && co == "HAALSI") "hi_binary" else NA_character_)
        }
      }
      break  # first usable vision column only
    }
  }
  logf("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

# ---------------------------------------------------------------
# 5. Outputs
# ---------------------------------------------------------------
res <- rbindlist(res_all, use.names = TRUE, fill = TRUE)
res <- res[order(cohort, outcome, method, substratum, exposure, model, caliber)]
fwrite(res, file.path(D_OUT, "STEP04_库内模型_v1.0.csv"))
qc <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
fwrite(qc, file.path(D_QC, "STEP04_分析集与排除_v1.0.csv"))
ph <- rbindlist(ph_rows, use.names = TRUE, fill = TRUE)
fwrite(ph, file.path(D_QC, "STEP04_PH检验_v1.0.csv"))
dose <- rbindlist(dose_rows, use.names = TRUE, fill = TRUE)
fwrite(dose, file.path(D_OUT, "STEP04_剂量反应_v1.0.csv"))
logf("[DONE] result rows = ", nrow(res), " | dose rows = ", nrow(dose))
cat("STEP04 DONE\n")
