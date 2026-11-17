# =====================================================================
# 29_step04b_cvd_endpoint.R  (v1.0, 2026-11-07)
# Topic 6 - SECONDARY ANALYSIS (user question 2026-11-07): outcome switched to
#           incident CVD composite (new heart disease or stroke).
#           Registered as secondary; does NOT alter the frozen main endpoints.
#
# Data   : landmark v1.1 (main-pool flags) + CVD时间_v1.0.csv (28_rebuild_cvd_time.py)
# Rules  : same as STEP04/05 (RC window 24m person-exclusion, MI m=20 covariates-only,
#          HRS half-sample + wave covariate, M1/M2/M3, REML+KH+PI pooling);
#          HAALSI composite = hearte-only -> in-cohort rows kept with note, EXCLUDED
#          from pooling (composite definition not comparable); MHAS has no HI/DSI.
# Output : 04_结果/STEP04b_库内模型_CVD终点_v1.0.csv
#          04_结果/STEP05b_合并_CVD终点_v1.0.csv
# =====================================================================

suppressWarnings(suppressMessages({
  library(data.table); library(survival); library(mice); library(metafor)
}))
SEED <- 20260904; M_IMP <- 20; RC_WINDOW <- 24
set.seed(SEED)
ROOT <- normalizePath(".", winslash = "/")
D_IN <- file.path(ROOT, "03_数据"); D_OUT <- file.path(ROOT, "04_结果"); D_QC <- file.path(ROOT, "05_质控")
LOGF <- file.path(D_QC, "STEP04b_run_log_20261107.txt")
logf <- function(...) try(cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n", file = LOGF, append = TRUE), silent = TRUE)
unlink(LOGF)
cat("STEP04b log ->", LOGF, "\n")
logf("metafor", as.character(packageVersion("metafor")), "| survival", as.character(packageVersion("survival")))

norm_id <- function(x) sub("^0+", "", sub("\\.0$", "", tolower(as.character(x))))

# ---------------- assembly ----------------
lm11 <- fread(file.path(D_IN, "中介_界标表_v1.1.csv"), colClasses = list(character = "id"))
stopifnot(nrow(lm11) == 50183)
lm11[, main := (cohort != "HRS") | (hrs_half == 1)]
lm11[cohort == "HRS", wave := base_wave]
cvd <- fread(file.path(D_IN, "CVD时间_v1.0.csv"), colClasses = list(character = "id"))
stopifnot(!any(duplicated(cvd[, .(cohort, id)])))
lm11 <- merge(lm11, cvd[, .(cohort, id, fu_cvd = fu_cvd_months, ev_cvd = cvd_event,
                            cvd_prev = cvd_prevalent, cvd_source)],
              by = c("cohort", "id"), all.x = TRUE, sort = FALSE)
lm11[, e_vi := fifelse(is.na(vi), NA_real_, as.numeric(vi == 1 & (is.na(hi) | hi == 0)))]
lm11[, e_hi := fifelse(is.na(hi), NA_real_, as.numeric(hi == 1 & (is.na(vi) | vi == 0)))]
lm11[, e_ds := fifelse(is.na(vi) | is.na(hi), NA_real_, as.numeric(vi == 1 & hi == 1))]
lm11[cohort == "MHAS", `:=`(e_hi = NA_real_, e_ds = NA_real_)]

M1V <- c("age_y", "ragender")
M2V_FULL <- c(M1V, "raeducel", "ses_inc_median", "ses_wealth_q")
M3V_FULL <- c(M2V_FULL, "currsmoke", "current_drinker", "pa_inactive", "cancre", "lunge", "arthre", "on_bp_med")
FACTORV <- c("ragender", "raeducel", "ses_wealth_q", "wave")
get_vars <- function(co, model) {
  v <- switch(model, M1 = M1V, M2 = M2V_FULL, M3 = M3V_FULL)
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
  clean <- clean[keep]; m <- length(clean)
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
  list(beta = bbar, se = se, p = pv)
}

cohs <- c("CHARLS", "ELSA", "SHARE", "KLoSA", "MHAS", "HAALSI", "HRS")
res_all <- list(); qc_rows <- list(); ph_rows <- list()
EXPLAB <- c(e_vi = "VI-only", e_hi = "HI-only", e_ds = "DSI")

for (co in cohs) {
  t0 <- Sys.time(); logf("=== ", co, " ===")
  uni <- lm11[cohort == co & main == TRUE]
  dset <- uni[!is.na(fu_cvd) & !is.na(ev_cvd) & !is.na(age_y) & cvd_prev == 0]
  drop_d <- dset[ev_cvd == 1 & fu_cvd < RC_WINDOW, .N]
  dset <- dset[!(ev_cvd == 1 & fu_cvd < RC_WINDOW)]
  dset[, `:=`(t = fu_cvd, ev = as.integer(ev_cvd == 1))]
  logf("  cvd n=", nrow(dset), " ev=", sum(dset$ev), " (early dropped ", drop_d, ")")
  qc_rows[[co]] <- data.table(cohort = co, n = nrow(dset), ev_pre = sum(dset$ev) + drop_d,
                              ev_dropped_early = drop_d, ev_used = sum(dset$ev),
                              cvd_source = dset$cvd_source[1])

  m3v <- get_vars(co, "M3"); dums <- exp_dums(co)
  mi_f <- as.data.frame(uni[, .SD, .SDcols = unique(c("id", m3v, dums))])
  mi_f <- fact(mi_f, m3v)
  ph_t <- uni$fu_cvd; ph_e <- uni$ev_cvd
  ph_t[is.na(ph_t)] <- median(ph_t, na.rm = TRUE); ph_e[is.na(ph_e)] <- 0
  mi_f$t <- ph_t; mi_f$ev <- ph_e
  mi_f$cumhaz <- mice::nelsonaalen(mi_f, t, ev)
  mi_f$aux_death <- fifelse(is.na(uni$death_after), 0, uni$death_after)
  imp_vars <- setdiff(names(mi_f)[vapply(mi_f, function(x) anyNA(x) && !all(is.na(x)), TRUE)],
                      c("id", "t", "ev", "cumhaz", "aux_death", dums))
  logf("  MI imputes: ", paste(imp_vars, collapse = ","))
  set.seed(SEED + match(co, cohs) * 1000)
  imp <- mice(mi_f[, setdiff(names(mi_f), "id"), drop = FALSE], m = M_IMP, maxit = 5, printFlag = FALSE)
  comps <- lapply(seq_len(M_IMP), function(i) { ci <- complete(imp, i); ci$id <- mi_f$id; ci })
  imp_copy <- function(fr, i) {
    ci <- comps[[i]]; d2 <- copy(fr)
    for (v in imp_vars) d2[[v]] <- ci[[v]][match(d2$id, ci$id)]
    fact(d2, m3v)
  }

  for (mod in c("M1", "M2", "M3")) {
    vars <- get_vars(co, mod)
    fits <- lapply(seq_len(M_IMP), function(i) {
      dd <- imp_copy(dset, i); dd <- dd[cc_keep(dd, c(dums, vars)), ]
      tryCatch(coxph(as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "), "+",
                                      paste(vars, collapse = " + "))), data = dd, ties = "efron", model = TRUE),
               error = function(e) e)
    })
    okf <- vapply(fits, function(x) !inherits(x, "error"), TRUE)
    if (sum(okf) < M_IMP / 2) { logf("  WARN ", mod, " failed"); next }
    fl <- fits[okf]
    p <- length(stats::coef(fl[[1]]))
    nn <- nrow(dset); ne <- sum(dset$ev)
    pf <- tryCatch(pool_rubin(lapply(fl, function(x) stats::coef(x)),
                              lapply(fl, function(x) vcov(x)), dfcom = nn - p,
                              tag = paste(co, mod)), error = function(e) e)
    if (inherits(pf, "error")) { logf("  WARN pool ", co, mod, ": ", conditionMessage(pf)); next }
    for (e in dums) {
      res_all[[length(res_all) + 1]] <- data.table(
        cohort = co, exposure = EXPLAB[e], outcome = "incident_CVD", method = "Cox",
        model = mod, caliber = "mi", substratum = "overall",
        model_set = "full", n = nn, n_events = ne,
        HR = exp(pf$beta[e]), lcl = exp(pf$beta[e] - 1.96 * pf$se[e]),
        ucl = exp(pf$beta[e] + 1.96 * pf$se[e]), p = pf$p[e],
        ph_p = NA_real_, rc_window = RC_WINDOW,
        note = if (co == "HAALSI") "cvd_hearte_only" else NA_character_)
    }
  }
  # PH (CC M3)
  vars <- get_vars(co, "M3")
  dd <- fact(as.data.frame(dset[, .SD, .SDcols = c("t", "ev", dums, vars)]), vars)
  dd <- dd[cc_keep(dd, c(dums, vars)), ]
  fc <- tryCatch(coxph(as.formula(paste("Surv(t, ev) ~", paste(dums, collapse = " + "), "+",
                                        paste(vars, collapse = " + "))), data = dd, ties = "efron", model = TRUE),
                 error = function(e) e)
  if (!inherits(fc, "error")) {
    zp <- tryCatch(cox.zph(fc), error = function(e) e)
    if (!inherits(zp, "error")) {
      tt <- as.data.frame(zp$table)
      ph_rows[[co]] <- data.table(cohort = co, term = rownames(tt), chisq = tt$chisq,
                                  df = tt$df, p = tt$p, n_cc = nrow(dd), n_events_cc = sum(dd$ev))
    }
  }
  logf("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

res <- rbindlist(res_all, use.names = TRUE, fill = TRUE)
fwrite(res, file.path(D_OUT, "STEP04b_库内模型_CVD终点_v1.0.csv"))
qc <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
fwrite(qc, file.path(D_QC, "STEP04b_分析集_CVD终点_v1.0.csv"))
ph <- rbindlist(ph_rows, use.names = TRUE, fill = TRUE)
fwrite(ph, file.path(D_QC, "STEP04b_PH_CVD终点_v1.0.csv"))

# ---------------- pooling (REML+KH+PI) ----------------
mi <- res[!is.na(HR)]
mi[, se := (log(ucl) - log(lcl)) / (2 * 1.96)]
pool_rows <- list()
for (ex in c("VI-only", "HI-only", "DSI")) {
  # composite-comparability: HAALSI (hearte-only) excluded from all pools
  d <- mi[exposure == ex & model == "M3" & cohort != "HAALSI"]
  k <- nrow(d)
  row <- data.table(exposure = ex, k = k, cohorts = paste(sort(d$cohort), collapse = "+"),
                    pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
                    p = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                    pi_lo = NA_real_, pi_hi = NA_real_,
                    loo_min = NA_real_, loo_max = NA_real_, note = NA_character_)
  if (k < 5) {
    row[, `:=`(desc_min = min(d$HR), desc_max = max(d$HR), note = "k_lt5_descriptive")]
    setnames(row, c("desc_min", "desc_max"), c("desc_min_HR", "desc_max_HR"))
    pool_rows[[length(pool_rows) + 1]] <- row; next
  }
  r <- tryCatch(rma(yi = log(HR), sei = se, data = as.data.frame(d), method = "REML", test = "knha"),
                error = function(e) e)
  if (inherits(r, "error")) { row[, note := paste0("rma_failed:", conditionMessage(r))];
    pool_rows[[length(pool_rows) + 1]] <- row; next }
  pi <- predict(r)
  row[, `:=`(pooled_HR = exp(r$b[1]), pooled_lcl = exp(r$ci.lb[1]), pooled_ucl = exp(r$ci.ub[1]),
             p = r$pval[1], I2 = r$I2, tau2 = r$tau2, pi_lo = exp(pi$pi.lb), pi_hi = exp(pi$pi.ub))]
  loo <- vapply(seq_len(k), function(i) {
    tryCatch(exp(rma(yi = log(HR), sei = se, data = as.data.frame(d[-i]),
                     method = "REML", test = "knha")$b[1]), error = function(e) NA_real_)
  }, numeric(1))
  row[, `:=`(loo_min = min(loo, na.rm = TRUE), loo_max = max(loo, na.rm = TRUE))]
  pool_rows[[length(pool_rows) + 1]] <- row
}
fin <- rbindlist(pool_rows, use.names = TRUE, fill = TRUE)
fwrite(fin, file.path(D_OUT, "STEP05b_合并_CVD终点_v1.0.csv"))
logf("[DONE] in-cohort rows = ", nrow(res), "; pooled rows = ", nrow(fin))
cat("STEP04b/05b DONE\n")
