# =============================================================================
# SW5  Modal quintile rule for method A and the six sensitivity analyses
# Spatial-weight protocol, Sections 2.4 (amendment of 14 October 2026) and 2.7
#
# Inputs: VAL_DERIVED/p2/p2_analysis.rds, VAL_DERIVED/linked/linked_<pr>.rds,
#         data/derived/sw3/postal_values.rds
#
# M   Modal quintile of method A: the quintile q maximizing pi_p(q); ties go to the quintile
#     closest to the mean-score quintile, then the lower one. Compared with the mean-score
#     quintile (A), method B and method C, overall and by stratum, with Bayesian bootstrap
#     intervals (postal codes, 200 replicates) for differences in misclassification and kappa.
# S1  Alternative weights for Link_Source D: the weight set that is not primary (with the default
#     primary pccf = PCCF+ selection probabilities, S1 uses the DA population-share weights, pop)
# S2  Method B centroid from latitude/longitude instead of projected coordinates (computed with the
#     pop weights, for which the lat/long centroid is available; identical to pccf except Link_Source D)
# S3  EVI windows (pixel, 3x3, 5x5) and years (2025, 2026, mean)
# S4  Addresses with a NAR building reference point only
# S5  Including urban letter-carrier addresses with a blockface reference point
# S6  MSE decomposition (Eq. 7) for every EVI measure and both MSDI scores, PCCF+ weights
#
# Outputs (data/derived/sw5/, or sw5/test_<pr>/): modal_quintiles.rds, m_modal.csv,
#   m_modal_bootstrap.csv, s1_alt_weights.csv, s1_expected_misclassification.csv,
#   s2_latlon_centroid.csv, s3_evi_measures.csv, s4_building_only.csv, s5_with_urban_blockface.csv,
#   s6_decomposition.csv, sw5_results.xlsx
# Checks (stop if any fails):
#   1. the modal quintile never has a larger expected misclassification than the mean-score
#      quintile (sum over postal codes of 1 - pi at the assigned quintile)
#   2. S1: for Link_Source D under equal weights, the expected misclassification of random
#      allocation, mean(1 - pi_p(t_i)), equals the observed misclassification of method C within
#      3.29 bootstrap standard errors (material and social)
#
# Run: SW_PROV=11 Rscript scripts/sw5_sensitivity.R   (test)
#      Rscript scripts/sw5_sensitivity.R             (national)
# Log: data/derived/sw5_sensitivity.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(writexl) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW    <- DERIVED
PROVS <- if (nzchar(Sys.getenv("SW_PROV"))) strsplit(Sys.getenv("SW_PROV"), ",")[[1]] else NULL
OUT   <- if (is.null(PROVS)) file.path(SW, "sw5") else file.path(SW, "sw5", paste0("test_", paste(PROVS, collapse = "_")))
LOG   <- file.path(SW, "sw5_sensitivity.log")
REPS  <- as.integer(Sys.getenv("SW_REPS", "200"))
PRIMARY <- Sys.getenv("SW_PRIMARY", "pccf"); stopifnot(PRIMARY %in% c("pccf", "pop"))
ALT <- setdiff(c("pccf", "pop"), PRIMARY)
KEEP  <- c("F", "C", "D", "5", "4")
N_MAIN <- 12410638L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261015)

msg <- function(...) {
  line <- gsub("  +", " ", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(...)))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK  " else "CHECK FAIL", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}
fmt <- function(n) format(n, big.mark = ",")
use <- function(file, objs) {
  for (e in parse(file.path(VALIDATION, "scripts", file), keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% objs)
      eval(e, globalenv())
  stopifnot(all(objs %in% ls(globalenv())))
}
use("p2_variables.R", c("R_EARTH", "haversine"))
use("p3_analysis.R", c("dm_levels", "rur_levels", "wkappa"))

EVI_M <- c(as.vector(outer(c("pix", "b1", "b2"), c(2025, 2026), paste, sep = "_")), paste0(c("pix", "b1", "b2"), "_mean"))
msg("SW5 start", if (is.null(PROVS)) "(national)" else paste0("(test: ", paste(PROVS, collapse = ","), ")"), "| reps", REPS, "| primary weights", PRIMARY)

# ---- 1. Address table: main sample plus urban letter-carrier blockface addresses (S5) ------------------
lf <- list.files(file.path(VAL_DERIVED, "linked"), "^linked_[0-9]+[.]rds$", full.names = TRUE)
if (!is.null(PROVS)) lf <- lf[sub("^linked_([0-9]+)[.]rds$", "\\1", basename(lf)) %in% PROVS]
pc_of <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE, .(ADDR_GUID, PCODE)]))
cols <- c("ADDR_GUID", "PROV_CODE", "ref_source", "delivery_mode", "rurality_sac", "pccf_Link_Source", "pccf_coded",
          "ref_lat", "ref_lon", "ref_DAUID", "dist_m", paste0("evi_ref_", EVI_M), paste0("evi_pccf_", EVI_M),
          "ref_SCOREMAT", "ref_SCORESOC", "ref_QuintMat", "ref_QuintSoc", "pccf_SCOREMAT", "pccf_SCORESOC", "pccf_QuintMat", "pccf_QuintSoc")
x <- readRDS(file.path(VAL_DERIVED, "p2", "p2_analysis.rds"))[pccf_coded == TRUE & pccf_Link_Source %chin% KEEP, cols, with = FALSE]
x <- x[pc_of, on = "ADDR_GUID", nomatch = NULL]; rm(pc_of); invisible(gc())
setnames(x, "pccf_Link_Source", "Link_Source")
x[, delivery_mode := factor(delivery_mode, dm_levels)]
x[, rurality := factor(fcase(rurality_sac %in% c("CA (tracted)", "CA (non-tracted)"), "CA",
                             rurality_sac == "Territories (outside CA)", "Territories",
                             default = as.character(rurality_sac)), rur_levels)]
x[, `:=`(province = factor(PROV_CODE), Link_Source = factor(Link_Source, KEEP), rurality_sac = NULL, pccf_coded = NULL)]
x[, main := !(ref_source == "blockface" & delivery_mode %in% dm_levels[1:3])]
x[, pidx := as.integer(factor(PCODE))]
npc <- max(x$pidx)
msg("Addresses:", fmt(nrow(x)), "| main sample:", fmt(x[, sum(main)]), "| urban letter-carrier blockface (S5):", fmt(x[, sum(!main)]))
if (is.null(PROVS)) check(x[, sum(main)] == N_MAIN, sprintf("main sample = %s addresses", fmt(N_MAIN)))

# ---- 2. Postal-code values and the modal quintile -------------------------------------------------------
pv <- readRDS(file.path(SW, "sw3", "postal_values.rds"))
if (!is.null(PROVS)) pv <- pv[pr %chin% PROVS]
mode_q <- function(P, qa) {            # argmax pi; ties -> closest to the mean-score quintile, then lower
  ok <- !is.na(P[, 1])
  mx <- apply(P, 1, max)
  S <- ifelse(P >= mx - 1e-9, abs(col(P) - ifelse(is.na(qa), 3L, qa)) * 10 + col(P), Inf)
  out <- max.col(-S, ties.method = "first"); out[!ok] <- NA_integer_; out
}
for (k in c("mat", "soc")) {
  P <- as.matrix(pv[, paste0("A_pi_", k, "_", 1:5), with = FALSE])
  pv[, (paste0("A_mode_", k)) := mode_q(P, get(paste0("A_q_", k)))]
  pm <- P[cbind(seq_len(nrow(P)), pmax(1L, pv[[paste0("A_mode_", k)]]))]
  pa <- P[cbind(seq_len(nrow(P)), pmax(1L, pv[[paste0("A_q_", k)]]))]
  ok <- !is.na(pm) & !is.na(pa)
  check(all(pm[ok] >= pa[ok] - 1e-12), sprintf("%s: modal quintile has the highest pi in every postal code (%s postal code x weight-set rows); %.2f%% differ from the mean-score quintile",
                                               k, fmt(sum(ok)), 100 * mean(pv[[paste0("A_mode_", k)]][ok] != pv[[paste0("A_q_", k)]][ok])))
}
rm(P); saveRDS(pv[, .(pr, PCODE, Link_Source, wset, A_q_mat, A_mode_mat, A_q_soc, A_mode_soc)], file.path(OUT, "modal_quintiles.rds"))

# attach columns from one weight set to x (by postal code), with an optional prefix
add_cols <- function(ws, vars, prefix = "") {
  d <- pv[wset == ws]; m <- match(x$PCODE, d$PCODE)
  for (v in vars) set(x, j = paste0(prefix, v), value = d[[v]][m])
  invisible(paste0(prefix, vars))
}
drop_cols <- function(v) x[, (intersect(v, names(x))) := NULL]

# ---- 3. Metric helpers -------------------------------------------------------------------------------------
GROUPS <- list(overall = character(0), Link_Source = "Link_Source", delivery_mode = "delivery_mode",
               rurality = "rurality", province = "province")
by_groups <- function(fun, groups = names(GROUPS)) rbindlist(lapply(groups, function(g) {
  r <- fun(GROUPS[[g]]); r[, grouping := g]; setcolorder(r, "grouping"); r }), fill = TRUE)
# EVI agreement of named assigned columns with reference column ref, on rows `rows` (paired sample)
evi_tab <- function(rows, ref, meth, groups = names(GROUPS)) {
  cc <- c(ref, meth)
  ok <- rows & Reduce(`&`, lapply(cc, function(v) !is.na(x[[v]])))
  by_groups(function(by) rbindlist(lapply(names(meth), function(M) x[ok, .(method = M, N = .N,
    bias = round(mean(get(meth[[M]]) - get(ref)), 4), mae = round(mean(abs(get(meth[[M]]) - get(ref))), 4),
    rmse = round(sqrt(mean((get(meth[[M]]) - get(ref))^2)), 4), pearson_r = round(if (.N > 2) cor(get(meth[[M]]), get(ref)) else NA_real_, 3),
    mean_delta_vs_C = if ("C" %in% names(meth) && M != "C") round(mean(abs(get(meth[[M]]) - get(ref)) - abs(get(meth[["C"]]) - get(ref))), 4) else NA_real_),
    keyby = by])), groups)
}
msdi_tab <- function(rows, t, meth, groups = names(GROUPS)) {
  ok <- rows & Reduce(`&`, lapply(c(t, unlist(meth)), function(v) !is.na(x[[v]])))
  by_groups(function(by) rbindlist(lapply(names(meth), function(M) x[ok, .(method = M, N = .N,
    pct_misclassified = round(100 * mean(get(meth[[M]]) != get(t)), 2),
    pct_more_deprived = round(100 * mean(get(meth[[M]]) > get(t)), 2), pct_less_deprived = round(100 * mean(get(meth[[M]]) < get(t)), 2),
    pct_off_2plus = round(100 * mean(abs(get(meth[[M]]) - get(t)) >= 2), 2), kappa_linear = round(wkappa(get(t), get(meth[[M]])), 4),
    pct_assigned_q1 = round(100 * mean(get(meth[[M]]) == 1), 2), pct_assigned_q5 = round(100 * mean(get(meth[[M]]) == 5), 2),
    pct_ref_q1 = round(100 * mean(get(t) == 1), 2), pct_ref_q5 = round(100 * mean(get(t) == 5), 2)), keyby = by])), groups)
}
kappa_counts <- function(cnt, t, a) {
  m <- matrix(0, 5, 5); m[cbind(t, a)] <- cnt; n <- sum(m)
  w <- 1 - abs(outer(1:5, 1:5, "-")) / 4
  pe <- sum(w * outer(rowSums(m), colSums(m))) / n^2
  (sum(w * m) / n - pe) / (1 - pe)
}
gid_of  <- function(d, g) if (g == "overall") rep(1L, nrow(d)) else as.integer(d[[g]])
glab_of <- function(d, g) if (g == "overall") "all" else levels(d[[g]])
# bootstrap of misclassification and kappa differences between pairs of methods
boot_msdi <- function(rows, t, meth, pairs, groups) {
  ok <- rows & Reduce(`&`, lapply(c(t, unlist(meth)), function(v) !is.na(x[[v]])))
  cnt <- rbindlist(lapply(names(meth), function(M)
    x[ok, .(cnt = .N), by = .(pidx, Link_Source, delivery_mode, rurality, t = get(t), a = get(meth[[M]]))][, method := M]))
  q <- function(w) rbindlist(lapply(groups, function(g) {
    kk <- data.table(g = gid_of(cnt, g), method = cnt$method, t = cnt$t, a = cnt$a, wc = cnt$cnt * w[cnt$pidx])[, .(wc = sum(wc)), by = .(g, method, t, a)]
    st <- kk[, .(mis = sum(wc[t != a]) / sum(wc), kap = kappa_counts(wc, t, a)), by = .(g, method)]
    rbindlist(lapply(pairs, function(p) {
      a1 <- st[method == p[1]]; a2 <- st[method == p[2]][a1, on = "g"]
      data.table(grouping = g, group = glab_of(cnt, g)[a2$g], comparison = paste(p[1], "minus", p[2]),
                 misclass_diff_pct = 100 * (a2$i.mis - a2$mis), kappa_diff = a2$i.kap - a2$kap)
    }))
  }))
  est <- q(rep(1, npc))
  bt <- rbindlist(lapply(seq_len(REPS), function(b) q(rexp(npc))[, b := b]))
  bs <- bt[, .(mis_lo = quantile(misclass_diff_pct, .025, na.rm = TRUE), mis_hi = quantile(misclass_diff_pct, .975, na.rm = TRUE),
               kap_lo = quantile(kappa_diff, .025, na.rm = TRUE), kap_hi = quantile(kappa_diff, .975, na.rm = TRUE)), by = .(grouping, group, comparison)]
  merge(est, bs, by = c("grouping", "group", "comparison"))
}
out <- list()

# ---- M. Modal quintile rule ---------------------------------------------------------------------------------
v <- c(add_cols(PRIMARY, c("A_q_mat", "A_mode_mat", "B_q_mat", "A_q_soc", "A_mode_soc", "B_q_soc")))
for (k in c("mat", "soc")) {
  K <- tools::toTitleCase(k)
  meth <- list(A_mean = paste0("A_q_", k), A_mode = paste0("A_mode_", k), B = paste0("B_q_", k), C = paste0("pccf_Quint", K))
  out[[paste0("m_modal_", k)]] <- msdi_tab(x$main, paste0("ref_Quint", K), meth)[, dimension := k]
  t0 <- Sys.time()
  out[[paste0("m_modal_bootstrap_", k)]] <- boot_msdi(x$main, paste0("ref_Quint", K), meth,
    list(c("A_mode", "C"), c("A_mode", "A_mean"), c("A_mean", "C")), c("overall", "Link_Source", "delivery_mode", "rurality"))[, dimension := k]
  msg("M:", k, "done; bootstrap", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
}
m <- out$m_modal_mat[grouping == "overall"]
msg(sprintf("M: material misclassified overall: A mean %.2f%%, A mode %.2f%%, B %.2f%%, C %.2f%%; kappa A mean %.3f, A mode %.3f, C %.3f",
            m[method == "A_mean", pct_misclassified], m[method == "A_mode", pct_misclassified], m[method == "B", pct_misclassified],
            m[method == "C", pct_misclassified], m[method == "A_mean", kappa_linear], m[method == "A_mode", kappa_linear], m[method == "C", kappa_linear]))
drop_cols(v)

# ---- S1. Alternative weights for Link_Source D ---------------------------------------------------------------
v <- c(add_cols(ALT, c("A_evi_b1_mean", "B_evi_b1_mean", "A_q_mat", "A_mode_mat", "B_q_mat", "A_q_soc", "A_mode_soc", "B_q_soc"), prefix = "w_"),
       add_cols("pccf", c(paste0("A_pi_mat_", 1:5), paste0("A_pi_soc_", 1:5)), prefix = "e_"))
s1e <- evi_tab(x$main, "evi_ref_b1_mean", list(A = "w_A_evi_b1_mean", B = "w_B_evi_b1_mean", C = "evi_pccf_b1_mean"), c("overall", "Link_Source", "delivery_mode"))
s1m <- rbindlist(lapply(c("mat", "soc"), function(k) { K <- tools::toTitleCase(k)
  msdi_tab(x$main, paste0("ref_Quint", K), list(A_mean = paste0("w_A_q_", k), A_mode = paste0("w_A_mode_", k), B = paste0("w_B_q_", k),
           C = paste0("pccf_Quint", K)), c("overall", "Link_Source", "delivery_mode"))[, dimension := k] }))
out$s1_alt_weights <- rbind(s1e[, `:=`(measure = "evi_b1_mean", weights = ALT)], s1m[, `:=`(measure = "quintile", weights = ALT)], fill = TRUE)
# check: expected misclassification of random allocation = observed for method C (Link_Source D)
s1x <- rbindlist(lapply(c("mat", "soc"), function(k) {
  K <- tools::toTitleCase(k); t <- x[[paste0("ref_Quint", K)]]
  P <- as.matrix(x[, paste0("e_A_pi_", k, "_", 1:5), with = FALSE])
  em <- 1 - P[cbind(seq_len(nrow(P)), pmax(1L, ifelse(is.na(t), 1L, t)))]
  ok <- x$main & !is.na(t) & !is.na(em) & !is.na(x[[paste0("pccf_Quint", K)]])
  d <- data.table(pidx = x$pidx[ok], Link_Source = x$Link_Source[ok], em = em[ok], om = as.numeric(x[[paste0("pccf_Quint", K)]][ok] != t[ok]))
  a <- d[, .(n = .N, se = sum(em), so = sum(om)), by = .(pidx, Link_Source)]
  bt <- vapply(seq_len(REPS), function(b) { w <- rexp(npc)[a$pidx]
    r <- rowsum(cbind(a$n, a$se, a$so) * w, as.integer(a$Link_Source)); (r[, 3] - r[, 2]) / r[, 1] }, numeric(length(unique(a$Link_Source))))
  r <- a[, .(N = sum(n), expected_pct = 100 * sum(se) / sum(n), observed_C_pct = 100 * sum(so) / sum(n)), keyby = Link_Source]
  r[, `:=`(dimension = k, diff_pct = observed_C_pct - expected_pct, boot_se_pct = 100 * apply(matrix(bt, nrow = nrow(r)), 1, sd))]
  r[, z := diff_pct / boot_se_pct][]
}))
out$s1_expected_misclassification <- s1x
for (k in c("mat", "soc")) { r <- s1x[dimension == k & Link_Source == "D"]
  check(abs(r$z) <= 3.29, sprintf("S1 (%s): Link_Source D equal-weight expected misclassification %.2f%% vs observed C %.2f%% (z = %.2f)",
                                   k, r$expected_pct, r$observed_C_pct, r$z)) }
drop_cols(v)

# ---- S2. Centroid from latitude/longitude ---------------------------------------------------------------------
v <- c(add_cols("pop", c("B_lat", "B_lon", "B_evi_b1_mean", "B_q_mat", "B_q_soc", "B_DAuid", "B_ll_lat", "B_ll_lon", "B_ll_evi_b1_mean",
                       "B_ll_q_mat", "B_ll_q_soc", "B_ll_DAuid")))
x[, `:=`(dB = haversine(ref_lat, ref_lon, B_lat, B_lon), dBll = haversine(ref_lat, ref_lon, B_ll_lat, B_ll_lon))]
s2p <- by_groups(function(by) x[(main), .(N = .N, B_median_m = round(median(dB)), Bll_median_m = round(median(dBll)),
  B_mean_m = round(mean(dB)), Bll_mean_m = round(mean(dBll)), B_pct_same_da = round(100 * mean(B_DAuid == ref_DAUID), 2),
  Bll_pct_same_da = round(100 * mean(B_ll_DAuid == ref_DAUID), 2), pct_centroid_moved_gt_10m = round(100 * mean(haversine(B_lat, B_lon, B_ll_lat, B_ll_lon) > 10), 2)),
  keyby = by], c("overall", "Link_Source", "delivery_mode"))
s2e <- evi_tab(x$main, "evi_ref_b1_mean", list(B = "B_evi_b1_mean", B_latlon = "B_ll_evi_b1_mean"), c("overall", "Link_Source", "delivery_mode"))
s2m <- rbindlist(lapply(c("mat", "soc"), function(k) msdi_tab(x$main, paste0("ref_Quint", tools::toTitleCase(k)),
  list(B = paste0("B_q_", k), B_latlon = paste0("B_ll_q_", k)), c("overall", "Link_Source", "delivery_mode"))[, dimension := k]))
out$s2_latlon_centroid <- rbind(s2p[, part := "position"], s2e[, part := "evi_b1_mean"], s2m[, part := "quintile"], fill = TRUE)
drop_cols(c(v, "dB", "dBll"))

# ---- S3. EVI windows and years ---------------------------------------------------------------------------------
v <- c(add_cols(PRIMARY, c(paste0("A_evi_", EVI_M), paste0("B_evi_", EVI_M))))
out$s3_evi_measures <- rbindlist(lapply(EVI_M, function(mm)
  evi_tab(x$main, paste0("evi_ref_", mm), list(A = paste0("A_evi_", mm), B = paste0("B_evi_", mm), C = paste0("evi_pccf_", mm)),
          c("overall", "Link_Source", "delivery_mode"))[, measure := mm]))
drop_cols(v)

# ---- S4 and S5. Building-only sample; with urban letter-carrier blockface addresses ------------------------------
v <- c(add_cols(PRIMARY, c("A_evi_b1_mean", "B_evi_b1_mean", "A_q_mat", "A_mode_mat", "B_q_mat", "A_q_soc", "A_mode_soc", "B_q_soc")))
core <- function(rows) rbind(
  evi_tab(rows, "evi_ref_b1_mean", list(A = "A_evi_b1_mean", B = "B_evi_b1_mean", C = "evi_pccf_b1_mean"), c("overall", "Link_Source", "delivery_mode"))[, measure := "evi_b1_mean"],
  rbindlist(lapply(c("mat", "soc"), function(k) { K <- tools::toTitleCase(k)
    msdi_tab(rows, paste0("ref_Quint", K), list(A_mean = paste0("A_q_", k), A_mode = paste0("A_mode_", k), B = paste0("B_q_", k),
             C = paste0("pccf_Quint", K)), c("overall", "Link_Source", "delivery_mode"))[, measure := paste0("quintile_", k)] })), fill = TRUE)
out$s4_building_only <- core(x$main & x$ref_source == "building")
out$s5_with_urban_blockface <- core(rep(TRUE, nrow(x)))
drop_cols(v)
msg("S1-S5 done")

# ---- S6. MSE decomposition for every measure (PCCF+ weights) --------------------------------------------------
vars <- c(paste0("A_evi_", EVI_M), paste0("A_V_evi_", EVI_M), paste0("A_miss_evi_", EVI_M), "A_S_mat", "A_V_S_mat", "A_miss_S_mat",
          "A_S_soc", "A_V_S_soc", "A_miss_S_soc")
v <- add_cols("pccf", vars, prefix = "w_")
meas <- c(setNames(EVI_M, paste0("evi_", EVI_M)), score_mat = "S_mat", score_soc = "S_soc")
s6 <- rbindlist(lapply(names(meas), function(nm) {
  mm <- meas[[nm]]; isS <- grepl("^S_", mm)
  D <- if (isS) paste0("ref_SCORE", toupper(sub("S_", "", mm))) else paste0("evi_ref_", mm)
  C <- if (isS) paste0("pccf_SCORE", toupper(sub("S_", "", mm))) else paste0("evi_pccf_", mm)
  A <- paste0("w_A_", if (isS) mm else paste0("evi_", mm)); V <- paste0("w_A_V_", if (isS) mm else paste0("evi_", mm))
  Ms <- paste0("w_A_miss_", if (isS) mm else paste0("evi_", mm))
  ok <- x$main & !is.na(x[[D]]) & !is.na(x[[C]]) & !is.na(x[[A]]) & x[[Ms]] < 1e-12 & x$Link_Source != "F"
  d <- data.table(pidx = x$pidx[ok], grp = fifelse(x$Link_Source[ok] == "D", "D (PCCF+ weights)", "C/5/4"),
                  seC = (x[[C]][ok] - x[[D]][ok])^2, seA = (x[[A]][ok] - x[[D]][ok])^2, V = x[[V]][ok])
  a <- d[, .(n = .N, seC = sum(seC), seA = sum(seA), V = sum(V)), by = .(pidx, grp)]
  r <- a[, .(N = sum(n), mse_C = sum(seC) / sum(n), mse_A = sum(seA) / sum(n), mean_V = sum(V) / sum(n)), keyby = grp]
  bt <- vapply(seq_len(REPS), function(b) { w <- rexp(npc)[a$pidx]
    s <- rowsum(cbind(a$n, a$seC, a$seA, a$V) * w, a$grp); (s[, 2] - s[, 3] - s[, 4]) / s[, 1] }, numeric(nrow(r)))
  r[, `:=`(measure = nm, diff = mse_C - mse_A - mean_V, boot_se = apply(matrix(bt, nrow = nrow(r)), 1, sd))]
  r[, `:=`(rel_diff_pct = 100 * diff / mse_C, z = diff / boot_se, share_of_mse_C_from_V = mean_V / mse_C)][]
}))
out$s6_decomposition <- s6
big <- s6[abs(z) > 3.29]
msg("S6: decomposition |z| > 3.29 for", nrow(big), "of", nrow(s6), "measure x group rows", if (nrow(big)) paste0("(", paste(big$measure, big$grp, collapse = "; "), ")") else "")
drop_cols(v)

# ---- Write ----------------------------------------------------------------------------------------------------
tabs <- list(m_modal = rbind(out$m_modal_mat, out$m_modal_soc), m_modal_bootstrap = rbind(out$m_modal_bootstrap_mat, out$m_modal_bootstrap_soc),
             s1_alt_weights = out$s1_alt_weights, s1_expected_misclassification = out$s1_expected_misclassification,
             s2_latlon_centroid = out$s2_latlon_centroid, s3_evi_measures = out$s3_evi_measures,
             s4_building_only = out$s4_building_only, s5_with_urban_blockface = out$s5_with_urban_blockface, s6_decomposition = out$s6_decomposition)
for (n in names(tabs)) fwrite(tabs[[n]], file.path(OUT, paste0(n, ".csv")))
write_xlsx(lapply(tabs, as.data.frame), file.path(OUT, "sw5_results.xlsx"))
msg("SW5 done")
