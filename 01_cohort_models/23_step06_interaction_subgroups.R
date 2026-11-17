# =====================================================================
# 23_step06_interaction_subgroups.R  (v1.0, 2026-11-06)
# Topic 6 - STEP06 interaction (multiplicative + additive) and subgroups
#
# Frozen design: v0.3 8.4 / v1.0 10.2 (H3 verdict step)
#   1) Multiplicative interaction: exposure dummy x CKM stage (0-1 ref)
#      per cohort (Wald P) -> pooled interaction HR (REML+KH); MI m=20
#      Death = Cox; stroke = Fine-Gray (main) + cause-specific Cox
#   2) Additive interaction: RERI/AP/SI on HR scale (delta-method CI),
#      modifier = CKM stage >=2 vs 0-1; death endpoint; pooled RERI via REML
#   3) Subgroups (death endpoint, MI): age <65/>=65, sex, SES tertile
#      (within-cohort wealth tertiles), country income (cohort-level:
#      UMIC = CHARLS/MHAS/HAALSI, HIC = ELSA/SHARE/KLoSA/HRS - re-pool only);
#      per-cohort exp x subgroup Wald P, pooled per stratum
#   MI machinery identical to 20_step04 (same seeds -> identical imputations)
# Output: 04_结果/STEP06_交互与亚组_v1.0.csv
# ASCII-only source (GBK console).
# =====================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(survival); library(cmprsk); library(mice); library(metafor)
}))

SEED <- 20260904
M_IMP <- 20
RC_WINDOW <- 24
set.seed(SEED)

ROOT <- normalizePath(".", winslash = "/")
D_IN <- file.path(ROOT, "03_数据"); D_OUT <- file.path(ROOT, "04_结果"); D_QC <- file.path(ROOT, "05_质控")
LOGF <- file.path(D_QC, "STEP06_run_log_20261106.txt")
logf <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "
", file = LOGF, append = TRUE), silent = TRUE)
unlink(LOGF)
cat("STEP06 log ->", LOGF, "\n")
logf("metafor", as.character(packageVersion("metafor")), "| survival", as.character(packageVersion("survival")))

MASTER <- "D:/数据/国际老年数据库/00_跨库资源/10_主题研究/06_生物标志物Biomarker/002_10年心血管风险计算&CKM/03_数据/ALL_DB_CKM_Master_v1.6.csv"
COV46  <- "D:/数据/国际老年数据库/00_跨库资源/10_主题研究/03_通用协变量CommonCovariates/03_输出数据/ALL_DB_Covariates_Long_v4.6.csv"

norm_id <- function(x) sub("^0+", "", sub("\\.0$", "", tolower(as.character(x))))

# ---------------------------------------------------------------
# 1. Data assembly (identical to 20_step04, incl. outcome times v2.1/v1.0)
# ---------------------------------------------------------------
lm <- fread(file.path(D_IN, "中介_界标表_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(nrow(lm) == 50183)
dt <- fread(file.path(D_IN, "死亡时间_v2.3.csv"), colClasses = list(character = "id"))
stopifnot(dt[, max(uniqueN(paste(fu_months, event))), by = .(cohort, id)]$V1 == 1)
dt <- unique(dt, by = c("cohort", "id"))
lm <- merge(lm, dt[, .(cohort, id, fu_d = fu_months, ev_d = event)],
            by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
sk <- fread(file.path(D_IN, "卒中时间_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(!any(duplicated(sk[, .(cohort, id)])), nrow(sk) == 50183)
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
cv <- fread(COV46, select = c("Sourcedatabase", "pid_unique", "wave", "hearte"),
            colClasses = list(character = "pid_unique"))
cv[, cohort := Sourcedatabase][, id := norm_id(sub("^[^_]+_", "", pid_unique))]
cv <- cv[cohort %in% c("CHARLS", "ELSA", "SHARE", "KLoSA", "MHAS", "HAALSI") & wave == BASE_WAVE[cohort]]
cv <- unique(cv[, .(cohort, id, hearte46 = hearte)])
lm <- merge(lm, cv, by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
lm[, heart_d := fifelse(cohort == "HAALSI",
                        fifelse(!is.na(heart_prev), heart_prev, hearte46), NA_real_)]
lm[, main := (cohort != "HRS") | (hrs_half == 1)]
stopifnot(lm[main == TRUE, .N] == 45949)
stopifnot(lm[main == TRUE, sum(ev_d == 1, na.rm = TRUE)] == 9438)  # v2.3 CHARLS 2026-11-12 + HAALSI W4 2026-11-13
lm[cohort == "HRS", wave := base_wave]

# exposure dummies + CKM stage indicators
lm[, e_vi := fifelse(is.na(vi), NA_real_, as.numeric(vi == 1 & (is.na(hi) | hi == 0)))]
lm[, e_hi := fifelse(is.na(hi), NA_real_, as.numeric(hi == 1 & (is.na(vi) | vi == 0)))]
lm[, e_ds := fifelse(is.na(vi) | is.na(hi), NA_real_, as.numeric(vi == 1 & hi == 1))]
lm[cohort == "MHAS", `:=`(e_hi = NA_real_, e_ds = NA_real_)]
lm[, stg2 := as.numeric(ckm_stage_asn == 2)]
lm[, stg3 := as.numeric(ckm_stage_asn == 3)]
lm[, stg23 := as.numeric(ckm_stage_asn %in% c(2, 3))]     # additive-interaction modifier
lm[, age65 := as.numeric(age_y >= 65)]
lm[, sex2 := as.numeric(ragender == 2)]

# SES tertile: within-cohort equal-frequency tertiles of the wealth ordinal
# (rank-based, robust to tied ordinal values)
lm[main == TRUE, w_tert := {
  q <- ses_wealth_q
  out <- rep(NA_integer_, .N)
  ok <- !is.na(q)
  if (any(ok)) {
    r <- rank(q[ok], ties.method = "first")
    out[ok] <- as.integer(cut(r, breaks = quantile(r, probs = c(0, 1/3, 2/3, 1)),
                              include.lowest = TRUE, labels = FALSE))
  }
  out
}, by = cohort]
logf("[1] assembly ok; SES tertile dist: ",
     paste(lm[main == TRUE, sort(unique(w_tert))], collapse = ","))

# ---------------------------------------------------------------
# 2. Model specs (same as STEP04 + user rulings)
# ---------------------------------------------------------------
M1V <- c("age_y", "ragender")
M2V_FULL <- c(M1V, "raeducel", "ses_inc_median", "ses_wealth_q")
M3V_FULL <- c(M2V_FULL, "currsmoke", "current_drinker", "pa_inactive", "cancre", "lunge", "arthre", "on_bp_med")
M2V_HAALSI <- c(M1V, "raeducel", "ses_wealth_q")
M3V_HAALSI <- c(M2V_HAALSI, "currsmoke", "current_drinker", "pa_inactive", "heart_d", "has_dm")
FACTORV <- c("ragender", "raeducel", "ses_wealth_q", "wave")
get_vars <- function(co, model = "M3") {
  v <- switch(model, M1 = M1V, M2 = if (co == "HAALSI") M2V_HAALSI else M2V_FULL,
              M3 = if (co == "HAALSI") M3V_HAALSI else M3V_FULL)
  if (co == "HRS") v <- c(v, "wave")
  v
}
exp_dums <- function(co) if (co == "MHAS") c("e_vi") else c("e_vi", "e_hi", "e_ds")
fact <- function(d, vars) { for (v in intersect(FACTORV, vars)) d[[v]] <- droplevels(factor(d[[v]])); d }
cc_keep <- function(d, vars) Reduce(`&`, lapply(vars, function(v) !is.na(d[[v]])))

pool_rubin <- function(betas, vcovs, dfcom, tag = "") {
  m0 <- length(betas)
  clean <- lapply(seq_len(m0), function(i) {
    b <- betas[[i]]; v <- vcovs[[i]]
    okn <- is.finite(b) & is.finite(diag(v))
    list(b = b[okn], v = v[okn, okn, drop = FALSE])
  })
  keep <- vapply(clean, function(x) length(x$b) > 0 && nrow(x$v) == length(x$b), TRUE)
  if (sum(keep) < max(5, ceiling(m0 / 2)))
    stop("pool_rubin [", tag, "]: too few usable imputations (", sum(keep), "/", m0, ")")
  clean <- clean[keep]
  m <- length(clean)
  betas <- lapply(clean, `[[`, "b"); vcovs <- lapply(clean, `[[`, "v")
  nms <- Reduce(intersect, lapply(betas, names))
  b  <- do.call(cbind, lapply(betas, function(x) x[nms]))
  vs <- lapply(vcovs, function(v) v[nms, nms, drop = FALSE])
  bbar <- rowMeans(b); names(bbar) <- nms
  Ubar <- Reduce(`+`, vs) / m
  Bm <- cov(t(b)); Tm <- Ubar + (1 + 1/m) * Bm
  lam <- pmin(pmax((1 + 1/m) * diag(Bm) / diag(Tm), 0), 0.999999)
  dfold <- (m - 1) / lam^2
  dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  df <- 1 / (1/dfold + 1/dfobs)
  se <- sqrt(diag(Tm)); names(se) <- nms
  pv <- 2 * pt(-abs(bbar / se), df); names(pv) <- nms
  list(beta = bbar, se = se, p = pv, df = df)
}
fit_coxf <- function(form, d) {
  warns <- character(0)
  f <- tryCatch(withCallingHandlers(coxph(form, data = d, ties = "efron", model = TRUE),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f)))
  list(ok = TRUE, fit = f)
}
fit_fg_X <- function(d, xvars) {   # d: data.table with t/fstatus; design from xvars columns
  cf <- fact(as.data.frame(d[, .SD, .SDcols = xvars]), xvars)
  ok <- cc_keep(cf, xvars); cf <- cf[ok, ]
  X <- model.matrix(~ ., cf)[, -1, drop = FALSE]
  sdv <- suppressWarnings(apply(X, 2, sd)); cst <- !is.finite(sdv) | sdv == 0
  if (any(cst)) X <- X[, !cst, drop = FALSE]
  dd <- d[ok, ]
  f <- tryCatch(crr(dd$t, dd$fstatus, X, failcode = 1, cencode = 0), error = function(e) e)
  if (inherits(f, "error")) return(list(ok = FALSE, msg = conditionMessage(f)))
  dimnames(f$var) <- list(colnames(X), colnames(X))
  list(ok = TRUE, coef = f$coef, var = f$var, n = nrow(X), nevent = sum(dd$fstatus == 1))
}

cohs <- c("CHARLS", "ELSA", "SHARE", "KLoSA", "MHAS", "HAALSI", "HRS")
OUTC <- c("cohort","block","outcome","method","exposure","level","n","n_events",
          "est","lcl","ucl","p","I2","tau2","pi_lo","pi_hi","k","note")
res <- list()

# ---------------------------------------------------------------
# 3. Per-cohort: MI + interaction/RERI/subgroup fits
# ---------------------------------------------------------------
COHORT_CSV <- file.path(D_OUT, "STEP06_逐库_交互与亚组_v1.0.csv")
REUSE <- Sys.getenv("STEP06_REUSE") == "1" && file.exists(COHORT_CSV)

if (!REUSE) for (co in cohs) {
  t0 <- Sys.time(); logf("=== ", co, " ===")
  uni <- lm[cohort == co & main == TRUE]
  dset <- uni[!is.na(fu_d) & !is.na(ev_d) & !is.na(age_y)]
  drop_d <- dset[ev_d == 1 & fu_d < RC_WINDOW, .N]
  dset <- dset[!(ev_d == 1 & fu_d < RC_WINDOW)]
  dset[, `:=`(t = fu_d, ev = as.integer(ev_d == 1))]
  sset0 <- uni[!is.na(fu_s) & !is.na(age_y) & !is.na(ev_s) & stroke_prev == 0]
  drop_s <- sset0[ev_s == 1 & fu_s < RC_WINDOW, .N]
  sset <- sset0[!(ev_s == 1 & fu_s < RC_WINDOW)]
  sset[, `:=`(t = fu_s, ev = as.integer(ev_s == 1),
              fstatus = fifelse(ev_s == 1, 1L, fifelse(!is.na(ev_d) & ev_d == 1, 2L, 0L)))]
  logf("  death n=", nrow(dset), " ev=", sum(dset$ev), " (early ", drop_d, ") | stroke n=",
       nrow(sset), " ev=", sum(sset$ev), " (early ", drop_s, ")")

  m3v <- get_vars(co); dums <- exp_dums(co)
  mi_f <- as.data.frame(uni[, .SD, .SDcols = unique(c("id", m3v, dums))])
  mi_f <- fact(mi_f, m3v)
  ph_t <- uni$fu_d; ph_e <- uni$ev_d
  ph_t[is.na(ph_t)] <- median(ph_t, na.rm = TRUE); ph_e[is.na(ph_e)] <- 0
  mi_f$t <- ph_t; mi_f$ev <- ph_e
  mi_f$cumhaz <- mice::nelsonaalen(mi_f, t, ev)
  mi_f$ev_s_pred <- fifelse(is.na(uni$ev_s), 0, uni$ev_s)
  imp_vars <- setdiff(names(mi_f)[vapply(mi_f, function(x) anyNA(x) && !all(is.na(x)), TRUE)],
                      c("id", "t", "ev", "cumhaz", "ev_s_pred", dums))
  set.seed(SEED + match(co, cohs) * 1000)
  imp <- mice(mi_f[, setdiff(names(mi_f), "id"), drop = FALSE], m = M_IMP, maxit = 5, printFlag = FALSE)
  comps <- lapply(seq_len(M_IMP), function(i) { ci <- complete(imp, i); ci$id <- mi_f$id; ci })
  imp_copy <- function(fr, i) {
    ci <- comps[[i]]; d2 <- copy(fr)
    for (v in imp_vars) d2[[v]] <- ci[[v]][match(d2$id, ci$id)]
    fact(d2, m3v)
  }
  logf("  MI done (", paste(imp_vars, collapse = ","), ")")

  # ---- 3a. multiplicative interaction: exp x stage (0-1 ref) ----
  for (oc in c("death", "stroke")) {
    # death -> Cox only (no competing risk); stroke -> FineGray main + CS Cox secondary
    for (me in if (oc == "death") c("Cox") else c("FineGray", "CauseSpecificCox")) {
      fr0 <- if (oc == "death") dset else sset
      fits <- lapply(seq_len(M_IMP), function(i) {
        dd <- imp_copy(fr0, i)
        need <- c(dums, m3v); ok <- cc_keep(dd, need)
        dd <- dd[ok, ]
        if (oc == "stroke" && me == "FineGray") {
          xv <- c(dums, "stg2", "stg3", paste0(rep(dums, 2), "_",
                 rep(c("stg2", "stg3"), each = length(dums))))
          dd[, e_vi_stg2 := e_vi * stg2]; dd[, e_hi_stg2 := e_hi * stg2]
          dd[, e_ds_stg2 := e_ds * stg2]
          dd[, e_vi_stg3 := e_vi * stg3]; dd[, e_hi_stg3 := e_hi * stg3]
          dd[, e_ds_stg3 := e_ds * stg3]
          fit_fg_X(dd, xv)
        } else {
          interms <- paste(vapply(dums, function(e) paste0(e, ":stg2 + ", e, ":stg3"), ""), collapse = " + ")
          form <- as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "), "+ stg2 + stg3 +",
                                   interms, "+", paste(m3v, collapse = " + ")))
          f <- fit_coxf(form, dd)
          if (isTRUE(f$ok)) list(ok = TRUE, coef = stats::coef(f$fit), var = vcov(f$fit),
                                 n = nrow(dd), nevent = sum(dd$ev)) else f
        }
      })
      okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
      if (sum(okf) < M_IMP / 2) { logf("  WARN int ", oc, me, " failed"); next }
      fl <- fits[okf]
      # cohort-level Rubin pool of interaction coefficients
      nms0 <- names(fl[[1]]$coef)
      ints <- nms0[grepl("[:_]stg", nms0)]
      for (tm in ints) {
        pf <- tryCatch(pool_rubin(lapply(fl, function(x) x$coef[tm]), lapply(fl, function(x) x$var[tm, tm, drop = FALSE]),
                         dfcom = fl[[1]]$nevent - length(fl[[1]]$coef), tag = paste(co, oc, me, tm)),
                       error = function(e) { logf("  NOTE: term not estimable: ", co, oc, me, tm, " (", conditionMessage(e), ")"); NULL })
        if (is.null(pf)) next
        res[[length(res) + 1]] <- data.table(cohort = co, block = "interaction", outcome = oc,
          method = me, exposure = sub(":stg.*", "", tm), level = sub(".*:", "", tm),
          n = fl[[1]]$n, n_events = fl[[1]]$nevent, est = exp(pf$beta[1]),
          lcl = exp(pf$beta[1] - 1.96 * pf$se[1]), ucl = exp(pf$beta[1] + 1.96 * pf$se[1]),
          p = pf$p[1], I2 = NA_real_, tau2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
          k = NA_integer_, note = NA_character_)
      }
    }
  }

  # ---- 3b. additive interaction (RERI/AP/SI), death, modifier = stage>=2 ----
  for (e in dums) {
    fits <- lapply(seq_len(M_IMP), function(i) {
      dd <- imp_copy(dset, i)
      dd <- dd[cc_keep(dd, c(e, m3v)), ]
      dd[, a_Eonly := get(e) * (1 - stg23)]
      dd[, b_Monly := (1 - get(e)) * stg23]
      dd[, c_both := get(e) * stg23]
      f <- fit_coxf(as.formula(paste("Surv(t, ev) ~ a_Eonly + b_Monly + c_both +",
                                     paste(m3v, collapse = " + "))), dd)
      if (isTRUE(f$ok)) list(ok = TRUE, coef = stats::coef(f$fit), var = vcov(f$fit),
                             n = nrow(dd), nevent = sum(dd$ev)) else f
    })
    okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
    if (sum(okf) < M_IMP / 2) { logf("  WARN RERI ", e, " failed"); next }
    fl <- fits[okf]
    tms <- c("a_Eonly", "b_Monly", "c_both")
    pf <- tryCatch(pool_rubin(lapply(fl, function(x) x$coef[tms]), lapply(fl, function(x) x$var[tms, tms, drop = FALSE]),
                     dfcom = fl[[1]]$nevent - length(fl[[1]]$coef), tag = paste(co, "RERI", e)),
                   error = function(e) { logf("  NOTE: RERI not estimable: ", co, e, " (", conditionMessage(e), ")"); NULL })
    if (is.null(pf)) next
    a <- pf$beta["a_Eonly"]; b <- pf$beta["b_Monly"]; c0 <- pf$beta["c_both"]
    # proper covariance from Rubin: recompute T matrix directly
    m <- length(fl)
    bs <- do.call(cbind, lapply(fl, function(x) x$coef[tms]))
    vs <- lapply(fl, function(x) x$var[tms, tms, drop = FALSE])
    Ubar <- Reduce(`+`, vs) / m; Bm <- cov(t(bs)); Tm <- Ubar + (1 + 1/m) * Bm
    ea <- exp(a); eb <- exp(b); ec <- exp(c0)
    reri <- ec - ea - eb + 1
    g <- c(-ea, -eb, ec)
    var_reri <- as.numeric(t(g) %*% Tm %*% g)
    ap <- reri / ec
    den <- (ea - 1) + (eb - 1)
    si <- if (den > 0) (ec - 1) / den else NA_real_
    res[[length(res) + 1]] <- data.table(cohort = co, block = "RERI", outcome = "death",
      method = "Cox", exposure = e, level = "stage23_vs_01",
      n = fl[[1]]$n, n_events = fl[[1]]$nevent, est = reri,
      lcl = reri - 1.96 * sqrt(var_reri), ucl = reri + 1.96 * sqrt(var_reri),
      p = 2 * pnorm(-abs(reri / sqrt(var_reri))),
      I2 = NA_real_, tau2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
      k = NA_integer_, note = paste0("AP=", round(ap, 3), ";SI=", ifelse(is.na(si), "NA", round(si, 3))))
  }

  # ---- 3c. subgroups (death, M3): age65, sex, SES tertile ----
  for (sg in c("age65", "sex2", "w_tert")) {
    levs <- if (sg == "w_tert") 1:3 else c(0, 1)
    for (lv in levs) {
      fits <- lapply(seq_len(M_IMP), function(i) {
        dd <- imp_copy(dset, i)
        dd <- dd[cc_keep(dd, c(dums, m3v)) & get(sg) == lv, ]
        f <- fit_coxf(as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "), "+",
                                       paste(m3v, collapse = " + "))), dd)
        if (isTRUE(f$ok)) list(ok = TRUE, coef = stats::coef(f$fit), var = vcov(f$fit),
                               n = nrow(dd), nevent = sum(dd$ev)) else f
      })
      okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
      if (sum(okf) < M_IMP / 2) next
      fl <- fits[okf]
      pf <- tryCatch(pool_rubin(lapply(fl, function(x) x$coef[dums]),
                       lapply(fl, function(x) x$var[dums, dums, drop = FALSE]),
                       dfcom = fl[[1]]$nevent - length(fl[[1]]$coef), tag = paste(co, sg, lv)),
                     error = function(e) { logf("  NOTE: subgroup not estimable: ", co, sg, lv, " (", conditionMessage(e), ")"); NULL })
      if (is.null(pf)) next
      for (e in dums) {
        res[[length(res) + 1]] <- data.table(cohort = co, block = paste0("subgroup_", sg),
          outcome = "death", method = "Cox", exposure = e, level = as.character(lv),
          n = fl[[1]]$n, n_events = fl[[1]]$nevent, est = exp(pf$beta[e]),
          lcl = exp(pf$beta[e] - 1.96 * pf$se[e]), ucl = exp(pf$beta[e] + 1.96 * pf$se[e]),
          p = pf$p[e], I2 = NA_real_, tau2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
          k = NA_integer_, note = NA_character_)
      }
    }
    # per-cohort multiplicative interaction P (exp x subgroup)
    for (lv_pair in list(c(0, 1))) {
      fits <- lapply(seq_len(M_IMP), function(i) {
        dd <- imp_copy(dset, i)
        dd <- dd[cc_keep(dd, c(dums, m3v, sg)), ]
        if (sg == "w_tert") {
          dd[, t2 := as.numeric(get(sg) == 2)]; dd[, t3 := as.numeric(get(sg) == 3)]
          iv <- paste(vapply(dums, function(e) paste0(e, ":t2 + ", e, ":t3"), ""), collapse = " + ")
          form <- as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "),
                                   "+ t2 + t3 +", iv, "+", paste(m3v, collapse = " + ")))
        } else {
          iv <- paste(paste0(dums, ":", sg), collapse = " + ")
          form <- as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "), "+", sg, "+",
                                   iv, "+", paste(m3v, collapse = " + ")))
        }
        f <- fit_coxf(form, dd)
        if (isTRUE(f$ok)) list(ok = TRUE, coef = stats::coef(f$fit), var = vcov(f$fit),
                               n = nrow(dd), nevent = sum(dd$ev)) else f
      })
      okf <- vapply(fits, function(x) isTRUE(x$ok), TRUE)
      if (sum(okf) < M_IMP / 2) next
      fl <- fits[okf]
      nms0 <- names(fl[[1]]$coef)
      ints <- nms0[grepl(":", nms0)]
      for (tm in ints) {
        pf <- tryCatch(pool_rubin(lapply(fl, function(x) x$coef[tm]), lapply(fl, function(x) x$var[tm, tm, drop = FALSE]),
                         dfcom = fl[[1]]$nevent - length(fl[[1]]$coef), tag = paste(co, sg, tm)),
                       error = function(e) { logf("  NOTE: term not estimable: ", co, sg, tm, " (", conditionMessage(e), ")"); NULL })
        if (is.null(pf)) next
        res[[length(res) + 1]] <- data.table(cohort = co, block = paste0("intP_", sg),
          outcome = "death", method = "Cox", exposure = sub(":.*", "", tm),
          level = sub(".*:", "", tm), n = fl[[1]]$n, n_events = fl[[1]]$nevent,
          est = exp(pf$beta[1]), lcl = exp(pf$beta[1] - 1.96 * pf$se[1]),
          ucl = exp(pf$beta[1] + 1.96 * pf$se[1]), p = pf$p[1],
          I2 = NA_real_, tau2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
          k = NA_integer_, note = NA_character_)
      }
    }
  }
  logf("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

# ---------------------------------------------------------------
# 4. Stage-2 pooling (REML+KH+PI) for interaction, RERI, subgroups
#    - yi is a real column (metafor cannot take an expression)
#    - RERI pooled on the LINEAR scale (delta-method se); HR-scale blocks on log scale
#    - one summary row per cell; cohort rows written separately
# ---------------------------------------------------------------
if (REUSE) {
  dt_all <- fread(COHORT_CSV)
  logf("[reuse] cohort rows loaded from ", basename(COHORT_CSV), " = ", nrow(dt_all))
} else {
  dt_all <- rbindlist(res, use.names = TRUE, fill = TRUE)
  fwrite(dt_all, COHORT_CSV)
}

# A3 pooling-policy alignment (2026-11-11/13): pooled HI/DSI cells follow the
# STEP05 pooling set - HAALSI (binary hearing construct) stays in the
# per-cohort file as a descriptive row but does not enter pooled e_hi/e_ds
# estimates (interaction, RERI, subgroups; k = 5 as in STEP05).
# Fine-Gray interaction rows carry compound exposure names (e_hi_stg2 etc.);
# normalize them first so the filter and the pooling keys are clean.
dt_all[block == "interaction" & grepl("_stg", exposure),
       `:=`(level = sub(".*_stg", "stg", exposure), exposure = sub("_stg.*", "", exposure))]
dt_all <- dt_all[!(cohort == "HAALSI" & exposure %in% c("e_hi", "e_ds"))]

pool_stage2 <- function(d, scale = c("log", "linear")) {
  scale <- match.arg(scale)
  d <- copy(d)[!is.na(est) & !is.na(lcl) & !is.na(ucl)]
  d[, se := if (scale == "log") (log(ucl) - log(lcl)) / (2 * 1.96)
             else (ucl - lcl) / (2 * 1.96)]
  d <- d[is.finite(se) & se > 0 & is.finite(est)]
  d[, yi := if (scale == "log") log(est) else est]
  k <- nrow(d)
  base <- data.table(block = d$block[1], outcome = d$outcome[1], method = d$method[1],
                     exposure = d$exposure[1], level = d$level[1], k_pool = k,
                     cohorts = paste(sort(d$cohort), collapse = "+"),
                     pooled_est = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
                     p_pool = NA_real_, I2_p = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
                     n_total = sum(d$n, na.rm = TRUE), ev_total = sum(d$n_events, na.rm = TRUE),
                     note_pool = NA_character_)
  if (k < 3) { base[, note_pool := "k_lt3_descriptive"]
               base[, `:=`(pooled_est = if (scale == "log") NA_real_ else mean(d$est),
                           note_pool = "k_lt3_descriptive_range_only")]
               return(base) }
  r <- tryCatch(rma(yi = yi, sei = se, data = as.data.frame(d), method = "REML", test = "knha"),
                error = function(e) e)
  if (inherits(r, "error")) {
    base[, note_pool := paste0("rma_failed:", conditionMessage(r))]
    return(base)
  }
  pi <- predict(r)
  tr <- if (scale == "log") exp else identity
  base[, `:=`(pooled_est = tr(r$b[1]), pooled_lcl = tr(r$ci.lb[1]), pooled_ucl = tr(r$ci.ub[1]),
              p_pool = r$pval[1], I2_p = r$I2, pi_lo = tr(pi$pi.lb), pi_hi = tr(pi$pi.ub),
              note_pool = NA_character_)]
  base
}

final <- list()
for (bl in unique(dt_all$block)) {
  b <- dt_all[block == bl]
  sc <- if (bl == "RERI") "linear" else "log"
  if (grepl("^subgroup_|^intP_", bl)) {
    for (ex in unique(b$exposure)) {
      for (lv in unique(b[exposure == ex, level])) {
        # subgroup blocks are per (exposure, level); intP blocks are the interaction terms
        final[[length(final) + 1]] <- pool_stage2(b[exposure == ex & level == lv], sc)
      }
    }
  } else {
    for (oc in unique(b$outcome)) {
      for (me in unique(b[b$outcome == oc, method])) {
        for (ex in unique(b[b$outcome == oc & method == me, exposure])) {
          sub <- b[outcome == oc & method == me & exposure == ex]
          if (bl == "interaction") {
            for (lv in unique(sub$level)) final[[length(final) + 1]] <- pool_stage2(sub[level == lv], sc)
          } else final[[length(final) + 1]] <- pool_stage2(sub, sc)
        }
      }
    }
  }
}
fin <- rbindlist(final, use.names = TRUE, fill = TRUE)
setorder(fin, block, outcome, method, exposure, level)
fwrite(fin, file.path(D_OUT, "STEP06_合并_交互与亚组_v1.0.csv"))
logf("[DONE] cohort rows = ", nrow(dt_all), "; pooled rows = ", nrow(fin),
     "; pooled estimated = ", sum(!is.na(fin$pooled_est)),
     "; k<3 descriptive = ", sum(fin$note_pool == "k_lt3_descriptive_range_only", na.rm = TRUE))
cat("STEP06 DONE\n")
