# 31_step07_repool.R - re-pool STEP07 from the (label-repaired) per-cohort CSV
suppressWarnings(suppressMessages({ library(data.table); library(metafor) }))
ROOT <- normalizePath(".", winslash = "/")
D_OUT <- file.path(ROOT, "04_结果")
res <- fread(file.path(D_OUT, "STEP07_敏感性核心_逐库_v1.0.csv"))
res[, se := (log(ucl) - log(lcl)) / (2 * 1.96)]

pool_rows <- list()
do_pool <- function(d, label) {
  for (oc in c("death", "stroke")) {
    for (ex in c("VI-only", "HI-only", "DSI")) {
      excl <- if (ex %in% c("HI-only", "DSI")) c("HAALSI", "MHAS") else character(0)
      sel <- d[block == label & outcome == oc & exposure == ex & model == "M3" &
               !cohort %in% excl & !is.na(HR)]
      if (nrow(sel) == 0) next
      k <- nrow(sel)
      row <- data.table(block = label, exposure = ex, outcome = oc, k = k,
        cohorts = paste(sort(sel$cohort), collapse = "+"),
        pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
        p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
        desc_min_HR = NA_real_, desc_max_HR = NA_real_, note = NA_character_)
      if (k < 5) {
        dmin <- min(sel$HR); dmax <- max(sel$HR)
        row[, `:=`(desc_min_HR = dmin, desc_max_HR = dmax, note = "k_lt5_descriptive")]
      } else {
        r <- tryCatch(rma(yi = log(HR), sei = se, data = as.data.frame(sel),
                          method = "REML", test = "knha"), error = function(e) e)
        if (!inherits(r, "error")) {
          pi <- predict(r)
          row[, `:=`(pooled_HR = exp(r$b[1]), pooled_lcl = exp(r$ci.lb[1]),
                     pooled_ucl = exp(r$ci.ub[1]), p = r$pval[1], I2 = r$I2,
                     pi_lo = exp(pi$pi.lb), pi_hi = exp(pi$pi.ub))]
        } else row[, note := "rma_failed"]
      }
      pool_rows[[length(pool_rows) + 1]] <<- row
    }
  }
}
for (bl in c("S1_cut5", "S2_w12", "S2_w36")) do_pool(res[!is.na(block) & block == bl], bl)

# S3: CS stroke M3 from STEP04 (substratum = overall only)
st4 <- fread(file.path(D_OUT, "STEP04_库内模型_v1.0.csv"))
cs <- st4[caliber == "mi" & model == "M3" & outcome == "stroke" & method == "CauseSpecificCox" &
          substratum == "overall" & !is.na(HR)]
cs[, se := (log(ucl) - log(lcl)) / (2 * 1.96)]
for (ex in c("VI-only", "HI-only", "DSI")) {
  excl <- if (ex %in% c("HI-only", "DSI")) c("HAALSI", "MHAS") else character(0)
  sel <- cs[exposure == c("VI-only"="e_vi", "HI-only"="e_hi", "DSI"="e_ds")[[ex]] & cohort %in% setdiff(c("CHARLS","ELSA","SHARE","KLoSA","MHAS","HAALSI","HRS"), excl)]
  k <- nrow(sel)
  row <- data.table(block = "S3_cs_swap", exposure = ex, outcome = "stroke", k = k,
    cohorts = paste(sort(sel$cohort), collapse = "+"),
    pooled_HR = NA_real_, pooled_lcl = NA_real_, pooled_ucl = NA_real_,
    p = NA_real_, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
    desc_min_HR = NA_real_, desc_max_HR = NA_real_, note = NA_character_)
  if (k < 5) { dmin <- min(sel$HR); dmax <- max(sel$HR)
    row[, `:=`(desc_min_HR = dmin, desc_max_HR = dmax, note = "k_lt5_descriptive")] }
  else {
    r <- tryCatch(rma(yi = log(HR), sei = se, data = as.data.frame(sel), method = "REML", test = "knha"),
                  error = function(e) e)
    if (!inherits(r, "error")) {
      pi <- predict(r)
      row[, `:=`(pooled_HR = exp(r$b[1]), pooled_lcl = exp(r$ci.lb[1]), pooled_ucl = exp(r$ci.ub[1]),
                 p = r$pval[1], I2 = r$I2, pi_lo = exp(pi$pi.lb), pi_hi = exp(pi$pi.ub))]
    } else row[, note := "rma_failed"]
  }
  pool_rows[[length(pool_rows) + 1]] <- row
}

# S19 rows (in-cohort, no pooling) appended as reference rows
s19 <- res[block == "S19_hrs_full"]
s19p <- s19[, .(block, exposure, outcome, k = 1L,
                 cohorts = cohort, pooled_HR = HR, pooled_lcl = lcl, pooled_ucl = ucl,
                 p = p, I2 = NA_real_, pi_lo = NA_real_, pi_hi = NA_real_,
                 desc_min_HR = NA_real_, desc_max_HR = NA_real_, note = "hrs_full_incohort")]
fin <- rbindlist(c(pool_rows, list(s19p)), use.names = TRUE, fill = TRUE)
fin <- fin[order(block, outcome, exposure)]
fwrite(fin, file.path(D_OUT, "STEP07_敏感性核心_池化_v1.0.csv"))
cat("repooled rows =", nrow(fin), "\n")
