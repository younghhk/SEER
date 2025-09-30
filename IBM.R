# Author: Grace Hong
#
# Purpose:
#   Provides a function to compute Incidence-Based Mortality (IBM) rates
#   and rate ratios (RR) from registry data.
#
# Features:
#   - Computes age-adjusted IBM rates using standardized weights.
#   - Confidence intervals (CIs) for rates:
#       • Fay–Feuer method (recommended for rare events, small counts)
#       • Tiwari modification (alternative small-count adjustment)
#   - Rate ratios (RR) between groups:
#       • Confidence intervals estimated using the Delta log-normal method.

compute_dsr_and_rr_for_subset <- function(df,
                                          idx1, idx2,
                                          label_group1,
                                          label_group2,
                                          subset_label,
                                          count_col   = "Count",
                                          pop_col     = "Population",
                                          weight_col  = "weight",
                                          alpha = 0.05,
                                          scale = 1e5,
                                          normalize_weights = TRUE,
                                          ci_method = c("fayfeuer","tiwari")) {
  ci_method <- match.arg(ci_method)
  
  # ---------- Helper: Directly Standardized Rate (DSR) and variance ----------
  dsr_components <- function(Count, Pop, w, normalize_weights = TRUE) {
    stopifnot(length(Count) == length(Pop), length(Pop) == length(w))
    if (length(Count) == 0) stop("Selected strata are empty (length 0).")
    if (normalize_weights) {
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) stop("Weights must sum to a positive finite value.")
      w <- w / sw
    }
    dsr     <- sum((Count / Pop) * w)
    var_dsr <- sum((Count / (Pop^2)) * (w^2))
    list(dsr = dsr, var = var_dsr, se = sqrt(var_dsr), w = w)
  }
  
  # ---------- Helper: SEER-style Confidence Intervals ----------
  # Based on Fay–Feuer or Tiwari methods
  seer_ci <- function(Count, Pop, std_w_norm, alpha = 0.05,
                      method = c("fayfeuer","tiwari")) {
    method <- match.arg(method)
    
    if (length(Count) == 0) return(c(lower = NA_real_, upper = NA_real_))
    
    wi <- std_w_norm / Pop                               # SEER definition
    rate <- sum((Count / Pop) * std_w_norm)              # DSR (unscaled)
    v <- sum((wi^2) * Count)                             # variance component
    
    ok <- (Pop > 0) | (Count > 0)                        # valid strata
    
    if (method == "fayfeuer") {
      wm <- max(wi, na.rm = TRUE)
      z  <- wm^2
    } else {
      wm <- mean(wi[ok], na.rm = TRUE)
      z  <- mean((wi^2)[ok], na.rm = TRUE)
    }
    
    if (!is.finite(rate) || rate <= 0 || !is.finite(v) || v <= 0 ||
        !is.finite(wm) || !is.finite(z)) {
      return(c(lower = NA_real_, upper = NA_real_))
    }
    
    # Lower CI
    df_lo  <- (2 * rate^2) / v
    mul_lo <- v / (2 * rate)
    lo <- mul_lo * qchisq(alpha/2, df = df_lo)
    
    # Upper CI
    df_hi  <- (2 * (rate + wm)^2) / (v + z)
    mul_hi <- (v + z) / (2 * (rate + wm))
    hi <- mul_hi * qchisq(1 - alpha/2, df = df_hi)
    
    c(lower = lo, upper = hi)
  }
  
  # ---------- Extract data for the two groups ----------
  Count1 <- df[[count_col]][idx1]
  Pop1   <- df[[pop_col]][idx1]
  w1     <- df[[weight_col]][idx1]
  Count2 <- df[[count_col]][idx2]
  Pop2   <- df[[pop_col]][idx2]
  w2     <- df[[weight_col]][idx2]
  
  if (length(Count1) == 0) stop("idx1 selects 0 rows for this subset.")
  if (length(Count2) == 0) stop("idx2 selects 0 rows for this subset.")
  if (length(Count1) != length(Count2)) stop("Groups must have the same number of strata.")
  if (any(!is.finite(w1)) || any(!is.finite(w2))) stop("Non-finite weights.")
  
  # ---------- Normalize weights ----------
  if (max(abs(w1 - w2)) > 1e-12) {
    warning("Standard weights differ between groups; using group 1's weights.")
  }
  std_w_norm <- if (normalize_weights) w1 / sum(w1) else w1
  
  # ---------- Compute DSR estimates ----------
  g1 <- dsr_components(Count1, Pop1, std_w_norm, normalize_weights = FALSE)
  g2 <- dsr_components(Count2, Pop2, std_w_norm, normalize_weights = FALSE)
  
  rate1 <- g1$dsr * scale
  rate2 <- g2$dsr * scale
  
  # ---------- Confidence Intervals for rates ----------
  ci1 <- seer_ci(Count1, Pop1, std_w_norm, alpha, method = ci_method) * scale
  ci2 <- seer_ci(Count2, Pop2, std_w_norm, alpha, method = ci_method) * scale
  ci_label <- if (ci_method == "tiwari") "SEER (Tiwari modification)" else "SEER (Fay–Feuer)"
  
  # ---------- Rate Ratio (RR) with Delta log-normal CI ----------
  RR       <- g2$dsr / g1$dsr
  z        <- qnorm(1 - alpha/2)
  se_logRR <- sqrt(g1$var / (g1$dsr^2) + g2$var / (g2$dsr^2))
  RR_lo    <- exp(log(RR) - z * se_logRR)
  RR_hi    <- exp(log(RR) + z * se_logRR)
  
  # ---------- Final outputs ----------
  rates <- data.frame(
    Group        = c(label_group1, label_group2),
    Strata_n     = c(length(Count1), length(Count2)),
    DSR          = c(g1$dsr, g2$dsr),
    Var_DSR      = c(g1$var, g2$var),
    Rate_per1e5  = c(rate1, rate2),
    Rate_CI_low  = c(ci1[1], ci2[1]),
    Rate_CI_high = c(ci1[2], ci2[2]),
    CI_method    = ci_label,
    row.names = NULL
  )
  
  rr_label <- paste0("RR of DSRs for ", subset_label)
  
  rr_of_dsrs <- data.frame(
    Measure   = rr_label,
    Estimate  = RR,
    CI_low    = RR_lo,
    CI_high   = RR_hi,
    CI_method = "Delta log-normal",
    row.names = NULL
  )
  
  list(rates = rates, rr_of_dsrs = rr_of_dsrs)
}

