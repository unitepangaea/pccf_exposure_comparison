# =============================================================================
# SW3  Postal-code values for methods A (weighted mean) and B (weighted centroid)
# Spatial-weight protocol, Sections 2.4-2.5 (Equations 3-5 and 7)
#
# Inputs: data/derived/sw2/cand_exposures.rds, centroids.rds, msdi_cutpoints.csv
#
# Weight sets (one row per postal code and weight set):
#   pop    primary: w_pop (PCCF+ WCF weights for C, 5, 4; DA population share for D)
#   pccf   PCCF+ 8D selection probabilities: equal weights for D (sensitivity 1), w_pop
#          otherwise. Under these weights method A is the expected value of method C for
#          every link source, so the MSE decomposition (Eq. 7) holds exactly.
#
# Method A, for each EVI measure (pix, b1, b2 x 2025, 2026, mean), weights renormalized over
# candidates with a non-missing value (Section 2.4):
#   A_evi_<m>    weighted mean                 A_V_<m>     allocation variance V_p (Eq. 7)
#   A_miss_<m>   share of weight with missing values
# MSDI (mat, soc): A_S_<d> weighted mean factor score; A_q_<d> its quintile from the national
#   cut points (Eq. 4); A_pi_<d>_1..5 weighted candidate quintile distribution (Eq. 5);
#   A_V_S_<d> weighted variance of candidate scores; A_miss_S_<d>; A_nq_<d> number of distinct
#   candidate quintiles; A_pimax_<d> largest pi.
# Method B: EVI, DAuid, scores and quintiles at the weighted centroid (sw2), for the same
#   weight set (pop: variant pop; pccf: variant equal for Link_Source D, pop otherwise);
#   B_ll_* = centroid from latitude/longitude (sensitivity 2), on pop rows only.
#
# Output (data/derived/sw3/, or sw3/test_<pr>/ when SW_PROV is set):
#   postal_values.rds, sw3_summary.csv
# Checks (stop if any fails):
#   1. weights renormalize correctly: A values lie within the range of candidate values,
#      V_p >= 0, pi sums to 1 wherever a candidate has a quintile
#   2. Link_Source F: method A and method B equal method C exactly (EVI, all windows and
#      years; MSDI scores and quintiles) for every Link_Source F address
#   3. one row per postal code and weight set; every postal code of SW1 present
#
# Run: SW_PROV=11 Rscript scripts/sw3_postal_values.R   (test)
#      Rscript scripts/sw3_postal_values.R             (national)
# Log: data/derived/sw3_postal_values.log
# =============================================================================

suppressPackageStartupMessages(library(data.table))

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW    <- DERIVED
IN    <- file.path(SW, "sw2")
PROVS <- if (nzchar(Sys.getenv("SW_PROV"))) strsplit(Sys.getenv("SW_PROV"), ",")[[1]] else NULL
OUT   <- if (is.null(PROVS)) file.path(SW, "sw3") else file.path(SW, "sw3", paste0("test_", paste(PROVS, collapse = "_")))
LOG   <- file.path(SW, "sw3_postal_values.log")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TOL <- 1e-12

msg <- function(...) {
  line <- gsub("  +", " ", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(...)))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK  " else "CHECK FAIL", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}
fmt <- function(n) format(n, big.mark = ",")
same <- function(u, v) (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & abs(u - v) < TOL)

EVI_M <- c(as.vector(outer(c("pix", "b1", "b2"), c(2025, 2026), paste, sep = "_")), paste0(c("pix", "b1", "b2"), "_mean"))
DIMS  <- c(mat = "MAT", soc = "SOC")
SCORE <- c(mat = "SCOREMAT", soc = "SCORESOC"); QUINT <- c(mat = "QuintMat", soc = "QuintSoc")

msg("SW3 start", if (is.null(PROVS)) "(national)" else paste0("(test: ", paste(PROVS, collapse = ","), ")"))
cand <- readRDS(file.path(IN, "cand_exposures.rds"))
cen  <- readRDS(file.path(IN, "centroids.rds"))
cuts <- fread(file.path(IN, "msdi_cutpoints.csv"))
if (!is.null(PROVS)) { cand <- cand[pr %chin% PROVS]; cen <- cen[pr %chin% PROVS] }
quint_of <- function(score, cut) ifelse(is.na(score), NA_integer_, 1L + rowSums(outer(score, cut, ">")))
msg("Candidates:", fmt(nrow(cand)), "| postal codes:", fmt(uniqueN(cand$PCODE)))

# ---- Method A ---------------------------------------------------------------------------------
# Weighted sums use data.table's grouped sum (GForce); the variance is a second pass around
# the weighted mean. Weights are renormalized over candidates with a value.
method_A <- function(w) {
  d <- cand[, .(pr, PCODE, Link_Source, n_cand, w = get(w))]
  out <- unique(d[, .(pr, PCODE, Link_Source, n_cand)])
  wmean_var <- function(e, name) {
    ok <- !is.na(e)
    d[, `:=`(wo = w * ok, e0 = fifelse(ok, e, 0))]
    s <- d[, .(sw = sum(wo), swe = sum(wo * e0), lo = min(fifelse(wo > 0, e0, Inf)), hi = max(fifelse(wo > 0, e0, -Inf))), by = PCODE]
    s[, mu := fifelse(sw > 0, swe / sw, NA_real_)]
    d[s, mu := i.mu, on = "PCODE"]
    v <- d[, .(sv = sum(wo * (e0 - fcoalesce(mu, 0))^2)), by = PCODE]
    s[v, V := fifelse(sw > 0, i.sv / sw, NA_real_), on = "PCODE"]
    out[s, (paste0("A_", name)) := i.mu, on = "PCODE"]
    out[s, (paste0("A_V_", name)) := i.V, on = "PCODE"]
    out[s, (paste0("A_miss_", name)) := 1 - i.sw, on = "PCODE"]
    s[, .(PCODE, mu, V, lo, hi)]
  }
  rng <- list()
  for (m in EVI_M) rng[[m]] <- wmean_var(cand[[paste0("evi_", m)]], paste0("evi_", m))
  for (dk in names(DIMS)) {
    rng[[dk]] <- wmean_var(cand[[SCORE[[dk]]]], paste0("S_", dk))
    cutv <- cuts[dimension == DIMS[[dk]]]$cut
    out[, (paste0("A_q_", dk)) := quint_of(get(paste0("A_S_", dk)), cutv)]
    q <- cand[[QUINT[[dk]]]]; okq <- !is.na(q)
    d[, `:=`(wq = w * okq, qq = fifelse(okq, q, 0L))]
    p <- d[, .(sw = sum(wq), p1 = sum(wq * (qq == 1L)), p2 = sum(wq * (qq == 2L)), p3 = sum(wq * (qq == 3L)),
               p4 = sum(wq * (qq == 4L)), p5 = sum(wq * (qq == 5L)), nq = uniqueN(qq[wq > 0])), by = PCODE]
    for (j in 1:5) p[, (paste0("p", j)) := fifelse(sw > 0, get(paste0("p", j)) / sw, NA_real_)]
    setnames(p, c(paste0("p", 1:5), "nq"), c(paste0("A_pi_", dk, "_", 1:5), paste0("A_nq_", dk)))
    out <- merge(out, p[, !"sw"], by = "PCODE", all.x = TRUE, sort = FALSE)
    out[, (paste0("A_pimax_", dk)) := do.call(pmax, .SD), .SDcols = paste0("A_pi_", dk, "_", 1:5)]
  }
  attr(out, "rng") <- rng
  out
}
A <- list(pop = method_A("w_pop"))
cand[, w_pccf := fifelse(Link_Source == "D", w_equal, w_pop)]
A$pccf <- method_A("w_pccf")
msg("Method A computed for", fmt(nrow(A$pop)), "postal codes x 2 weight sets")

# check 1: A within candidate range, V >= 0, pi sums to 1
for (ws in names(A)) {
  rng <- attr(A[[ws]], "rng")
  bad <- sum(vapply(rng, function(r) r[!is.na(mu), sum(mu < lo - 1e-9 | mu > hi + 1e-9 | V < -1e-12)], 1))
  pis <- c(A[[ws]][!is.na(A_pi_mat_1), rowSums(.SD), .SDcols = paste0("A_pi_mat_", 1:5)],
           A[[ws]][!is.na(A_pi_soc_1), rowSums(.SD), .SDcols = paste0("A_pi_soc_", 1:5)])
  check(bad == 0 && all(abs(pis - 1) < 1e-9),
        sprintf("%s weights: A within candidate range, V_p >= 0, pi sums to 1 (%d violations)", ws, bad))
  attr(A[[ws]], "rng") <- NULL
}

# ---- Method B -----------------------------------------------------------------------------------
bcols <- c(paste0("evi_", EVI_M), "DAuid", "da_method", "SCOREMAT", "SCORESOC", "QuintMat", "QuintSoc", "x", "y", "lat", "lon")
take_B <- function(cen_rows, prefix) {
  b <- cen_rows[, c("PCODE", bcols), with = FALSE]
  setnames(b, bcols, paste0(prefix, c(paste0("evi_", EVI_M), "DAuid", "da_method", "S_mat", "S_soc", "q_mat", "q_soc", "x", "y", "lat", "lon")))
  b
}
B_pop  <- take_B(cen[variant == "pop"], "B_")
B_pccf <- take_B(cen[(variant == "equal" & Link_Source == "D") | (variant == "pop" & Link_Source != "D")], "B_")
B_ll   <- take_B(cen[variant == "latlon"], "B_ll_")

pv <- rbind(
  merge(merge(A$pop, B_pop, by = "PCODE"), B_ll, by = "PCODE")[, wset := "pop"],
  merge(A$pccf, B_pccf, by = "PCODE")[, wset := "pccf"], fill = TRUE)
setcolorder(pv, c("pr", "PCODE", "Link_Source", "wset", "n_cand"))
setorder(pv, pr, PCODE, wset)
check(pv[, .N, by = .(PCODE, wset)][, all(N == 1)] && pv[wset == "pop", .N] == uniqueN(cand$PCODE) &&
        pv[wset == "pccf", .N] == uniqueN(cand$PCODE),
      sprintf("one row per postal code and weight set (%s postal codes)", fmt(uniqueN(cand$PCODE))))

# ---- Check 2: Link_Source F, A = B = C ----------------------------------------------------------
lf <- list.files(file.path(VAL_DERIVED, "linked"), "^linked_[0-9]+[.]rds$", full.names = TRUE)
if (!is.null(PROVS)) lf <- lf[sub("^linked_([0-9]+)[.]rds$", "\\1", basename(lf)) %in% PROVS]
a <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE & pccf_Link_Source == "F", .(ADDR_GUID, PCODE)]))
p2 <- readRDS(file.path(VAL_DERIVED, "p2", "p2_analysis.rds"))
p2 <- p2[, c("ADDR_GUID", paste0("evi_pccf_", EVI_M), "pccf_SCOREMAT", "pccf_SCORESOC", "pccf_QuintMat", "pccf_QuintSoc"), with = FALSE]
a <- p2[a, on = "ADDR_GUID", nomatch = NULL]; rm(p2); invisible(gc())
for (ws in c("pop", "pccf")) {
  f <- pv[wset == ws & Link_Source == "F"][a, on = "PCODE", nomatch = NULL]
  ok <- Reduce(`&`, c(
    lapply(EVI_M, function(m) same(f[[paste0("A_evi_", m)]], f[[paste0("evi_pccf_", m)]]) &
                              same(f[[paste0("B_evi_", m)]], f[[paste0("evi_pccf_", m)]])),
    list(same(f$A_S_mat, f$pccf_SCOREMAT), same(f$A_S_soc, f$pccf_SCORESOC),
         same(f$B_S_mat, f$pccf_SCOREMAT), same(f$B_S_soc, f$pccf_SCORESOC),
         same(f$A_q_mat, f$pccf_QuintMat), same(f$A_q_soc, f$pccf_QuintSoc),
         same(f$B_q_mat, f$pccf_QuintMat), same(f$B_q_soc, f$pccf_QuintSoc))))
  check(nrow(f) == nrow(a) && all(ok), sprintf("Link_Source F (%s weights): A = B = C for all %s addresses (%s mismatches)",
                                               ws, fmt(nrow(a)), fmt(sum(!ok))))
}

saveRDS(pv, file.path(OUT, "postal_values.rds"))
msg("Saved", file.path(OUT, "postal_values.rds"), "(", fmt(nrow(pv)), "rows,", ncol(pv), "columns )")

# ---- Summary ----------------------------------------------------------------------------------------
summ <- pv[, .(postal_codes = .N,
               mean_sd_evi_b1 = round(mean(sqrt(A_V_evi_b1_mean), na.rm = TRUE), 4),
               pct_V_evi_b1_pos = round(100 * mean(A_V_evi_b1_mean > 1e-12, na.rm = TRUE), 1),
               mean_miss_evi_b1 = round(mean(A_miss_evi_b1_mean), 5),
               mean_miss_S_mat = round(mean(A_miss_S_mat), 4),
               pct_multi_quint_mat = round(100 * mean(A_nq_mat > 1, na.rm = TRUE), 1),
               mean_pimax_mat = round(mean(A_pimax_mat, na.rm = TRUE), 3),
               pct_qA_ne_qB_mat = round(100 * mean(A_q_mat != B_q_mat, na.rm = TRUE), 1),
               mean_abs_evi_A_minus_B = round(mean(abs(A_evi_b1_mean - B_evi_b1_mean), na.rm = TRUE), 4)),
           keyby = .(wset, Link_Source)]
fwrite(summ, file.path(OUT, "sw3_summary.csv")); print(summ)
msg("SW3 done")
