# =============================================================================
# SW4  Comparison of methods A, B and C against the address reference (method D)
# Spatial-weight protocol, Sections 2.5-2.6 (Equations 7-10)
#
# Inputs: VAL_DERIVED/p2/p2_analysis.rds (validation study: reference D and PCCF+ draw C),
#         VAL_DERIVED/linked/linked_<pr>.rds (postal codes), data/derived/sw1 (candidates,
#         postal-code summaries), data/derived/sw3/postal_values.rds (A and B)
# Sample: validation main sample (PCCF+ coded; urban letter-carrier addresses with a blockface
#         reference excluded) restricted to Link_Source F, C, D, 5, 4 (protocol 2.2).
#         Primary weights: SW_PRIMARY (default pccf = PCCF+ 8D selection probabilities, equal
#         weights for Link_Source D; protocol amendment of 14 October 2026; "pop" gives the original
#         DA population-share weights for D). Greenness: 3x3 EVI, 2025-26 mean. MSDI: national quintiles.
#
# Address-level outcomes (sw4/address_outcomes.rds):
#   positional error   C: haversine distance of the validation study (dist_m); B: haversine
#                      distance to the weighted centroid; A: expected planar distance
#                      sum_j w_pj d(r_i, x_pj) (Eq. 9); planar B and C distances for the
#                      inequality check; expected distance under PCCF+ weights = E[error of C]
#   DA agreement       C, B: same DA as the reference; A: weight of candidates in the reference DA
#   greenness          e^M - e^D for M in A, B, C; paired differences Delta_AC (primary), Delta_BC,
#                      Delta_AB (Eq. 10)
#   deprivation        quintile of each method vs reference; expected misclassification of
#                      random allocation 1 - pi_p(t_i); score differences
# Tables (sw4/*.csv and sw4/sw4_results.xlsx), overall and by Link_Source, delivery mode,
#   rurality, delivery mode x rurality and province:
#   t4_candidates, t5_position, t6_evi, t7_msdi, t7_quintile_shares, decomposition,
#   bootstrap (200 Bayesian bootstrap replicates over postal codes: mean Delta, differences in
#   misclassification and weighted kappa between methods), model coefficients, tests and
#   adjusted predictions (OLS for Delta_AC; logistic for material misclassification with method x
#   delivery mode x rurality + province), as in the validation study
# Checks (stop if any fails):
#   1. main sample = 12,410,638 addresses; every address has postal-code values
#   2. d(r_i, x_B) <= expected distance of A for every address (triangle inequality, Eq. 9)
#   3. MSE decomposition (Eq. 7) under PCCF+ weights, Link_Source C/5/4 and D: MSE(C) - [MSE(A) +
#      mean V_p] within 3.29 bootstrap standard errors, for EVI (3x3) and the material score
#
# Run: SW_PROV=11 Rscript scripts/sw4_compare.R   (test; sample-size check skipped)
#      Rscript scripts/sw4_compare.R             (national)
#      SW_REPS=200 bootstrap replicates (default 200)
# Log: data/derived/sw4_compare.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(writexl) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW    <- DERIVED
PROVS <- if (nzchar(Sys.getenv("SW_PROV"))) strsplit(Sys.getenv("SW_PROV"), ",")[[1]] else NULL
OUT   <- if (is.null(PROVS)) file.path(SW, "sw4") else file.path(SW, "sw4", paste0("test_", paste(PROVS, collapse = "_")))
LOG   <- file.path(SW, "sw4_compare.log")
REPS  <- as.integer(Sys.getenv("SW_REPS", "200"))
PRIMARY <- Sys.getenv("SW_PRIMARY", "pccf"); stopifnot(PRIMARY %in% c("pccf", "pop"))
N_MAIN <- 12410638L
KEEP  <- c("F", "C", "D", "5", "4")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261014)

msg <- function(...) {
  line <- gsub("  +", " ", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(...)))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK  " else "CHECK FAIL", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}
fmt <- function(n) format(n, big.mark = ",")

# ---- Reuse validation functions (parsed; only these assignments are evaluated) ----------------
use <- function(file, objs) {
  for (e in parse(file.path(VALIDATION, "scripts", file), keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% objs)
      eval(e, globalenv())
  stopifnot(all(objs %in% ls(globalenv())))
}
use("p2_variables.R", c("R_EARTH", "haversine"))
use("p3_analysis.R", c("MIN_CELL", "dm_levels", "rur_levels", "wkappa", "drop_sparse", "model_cells",
                       "ols_cells", "coef_table", "adj_pred", "model_formulas"))

msg("SW4 start", if (is.null(PROVS)) "(national)" else paste0("(test: ", paste(PROVS, collapse = ","), ")"), "| reps", REPS, "| primary weights", PRIMARY)

# ---- 1. Address sample -----------------------------------------------------------------------------
lf <- list.files(file.path(VAL_DERIVED, "linked"), "^linked_[0-9]+[.]rds$", full.names = TRUE)
if (!is.null(PROVS)) lf <- lf[sub("^linked_([0-9]+)[.]rds$", "\\1", basename(lf)) %in% PROVS]
pc_of <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE, .(ADDR_GUID, PCODE)]))
cols <- c("ADDR_GUID", "PROV_CODE", "ref_source", "delivery_mode", "rurality_sac", "pccf_Link_Source", "pccf_coded",
          "ref_x", "ref_y", "ref_lat", "ref_lon", "pccf_x", "pccf_y", "dist_m", "ref_DAUID", "same_da",
          "evi_ref_b1_mean", "evi_pccf_b1_mean", "ref_SCOREMAT", "ref_SCORESOC", "ref_QuintMat", "ref_QuintSoc",
          "pccf_SCOREMAT", "pccf_SCORESOC", "pccf_QuintMat", "pccf_QuintSoc")
x <- readRDS(file.path(VAL_DERIVED, "p2", "p2_analysis.rds"))[pccf_coded == TRUE, cols, with = FALSE]
x <- x[pc_of, on = "ADDR_GUID", nomatch = NULL]; rm(pc_of); invisible(gc())
setnames(x, "pccf_Link_Source", "Link_Source")
x[, delivery_mode := factor(delivery_mode, dm_levels)]
x[, rurality := factor(fcase(rurality_sac %in% c("CA (tracted)", "CA (non-tracted)"), "CA",
                             rurality_sac == "Territories (outside CA)", "Territories",
                             default = as.character(rurality_sac)), rur_levels)]
x[, province := factor(PROV_CODE)]
urban_lc <- x$delivery_mode %in% dm_levels[1:3]
x <- x[!(ref_source == "blockface" & urban_lc) & Link_Source %chin% KEEP]
x[, Link_Source := factor(Link_Source, KEEP)]
x[, aid := .I]
msg("Main sample, Link_Source F/C/D/5/4:", fmt(nrow(x)))
if (is.null(PROVS)) check(nrow(x) == N_MAIN, sprintf("main sample = %s addresses (protocol 2.2)", fmt(N_MAIN)))

# postal-code values (A, B) and candidate summaries (SW1)
pvcols <- c("PCODE", "A_evi_b1_mean", "A_V_evi_b1_mean", "A_miss_evi_b1_mean", "A_S_mat", "A_S_soc", "A_q_mat", "A_q_soc",
            "A_V_S_mat", "A_miss_S_mat", paste0("A_pi_mat_", 1:5), paste0("A_pi_soc_", 1:5), "A_nq_mat", "A_pimax_mat",
            "B_evi_b1_mean", "B_DAuid", "B_S_mat", "B_S_soc", "B_q_mat", "B_q_soc", "B_x", "B_y", "B_lat", "B_lon")
pv <- readRDS(file.path(SW, "sw3", "postal_values.rds"))
pvw <- pv[wset == "pccf", .(PCODE, Aw_evi_b1_mean = A_evi_b1_mean, Aw_V_evi_b1_mean = A_V_evi_b1_mean,
                            Aw_miss_evi_b1_mean = A_miss_evi_b1_mean, Aw_S_mat = A_S_mat, Aw_V_S_mat = A_V_S_mat, Aw_miss_S_mat = A_miss_S_mat)]
pv <- pv[wset == PRIMARY, pvcols, with = FALSE]
x <- pv[x, on = "PCODE"]; x <- pvw[x, on = "PCODE"]; rm(pv, pvw); invisible(gc())
check(x[, !anyNA(B_DAuid)], "every address has postal-code values for methods A and B")
pcs <- readRDS(file.path(SW, "sw1", "postal_codes.rds"))[, .(PCODE, n_cand, n_eff, n_da)]
x <- pcs[x, on = "PCODE"]; rm(pcs)

# ---- 2. Positional error and DA agreement -------------------------------------------------------------
cand <- readRDS(file.path(SW, "sw1", "candidates.rds"))[, .(PCODE, DAuid, cx = x, cy = y, w_pop, w_pccf = fifelse(Link_Source == "D", w_equal, w_pop))]
if (!is.null(PROVS)) cand <- cand[PCODE %chin% x$PCODE]
pos <- rbindlist(lapply(split(x[, .(aid, PCODE, ref_x, ref_y, ref_DAUID, pr = PROV_CODE)], by = "pr"), function(d) {
  m <- cand[d, on = "PCODE", allow.cartesian = TRUE, nomatch = NULL]
  m[, dd := sqrt((cx - ref_x)^2 + (cy - ref_y)^2)]
  m[, wA := if (PRIMARY == "pccf") w_pccf else w_pop]
  m[, .(dA_exp = sum(wA * dd), dC_exp_pccf = sum(w_pccf * dd), A_w_refda = sum(wA * (DAuid == ref_DAUID))), by = aid]
}))
rm(cand); invisible(gc())
x[pos, on = "aid", `:=`(dA_exp = i.dA_exp, dC_exp_pccf = i.dC_exp_pccf, A_w_refda = i.A_w_refda)]; rm(pos)
x[, `:=`(dB = haversine(ref_lat, ref_lon, B_lat, B_lon), dC = dist_m,
         dB_planar = sqrt((B_x - ref_x)^2 + (B_y - ref_y)^2), dC_planar = sqrt((pccf_x - ref_x)^2 + (pccf_y - ref_y)^2),
         same_da_B = B_DAuid == ref_DAUID, same_da_C = same_da)]
viol <- x[, sum(dB_planar > dA_exp + 1e-6)]
check(viol == 0, sprintf("distance to the weighted centroid <= expected distance of A for all %s addresses (%s violations)", fmt(nrow(x)), fmt(viol)))

# ---- 3. Greenness and deprivation outcomes ---------------------------------------------------------------
x[, `:=`(eD = evi_ref_b1_mean, eA = A_evi_b1_mean, eB = B_evi_b1_mean, eC = evi_pccf_b1_mean)]
x[, evi_ok := !is.na(eD) & !is.na(eA) & !is.na(eB) & !is.na(eC)]
x[(evi_ok), `:=`(dAC = abs(eA - eD) - abs(eC - eD), dBC = abs(eB - eD) - abs(eC - eD), dAB = abs(eA - eD) - abs(eB - eD))]
for (k in c("mat", "soc")) {
  K <- tools::toTitleCase(k); t <- x[[paste0("ref_Quint", K)]]
  x[, (paste0("q_ok_", k)) := !is.na(t) & !is.na(get(paste0("A_q_", k))) & !is.na(get(paste0("B_q_", k))) &
                               !is.na(get(paste0("pccf_Quint", K)))]
  pim <- as.matrix(x[, paste0("A_pi_", k, "_", 1:5), with = FALSE])
  x[, (paste0("C_expmis_", k)) := fifelse(is.na(t), NA_real_, 1 - pim[cbind(seq_len(.N), pmax(1L, as.integer(t)))])]
}
x[, pidx := as.integer(factor(PCODE))]
npc <- max(x$pidx)
x[, c(paste0("A_pi_", rep(c("mat", "soc"), each = 5), "_", 1:5), "ref_x", "ref_y", "pccf_x", "pccf_y", "B_x", "B_y") := NULL]
invisible(gc())
keep_out <- c("aid", "ADDR_GUID", "PCODE", "pidx", "Link_Source", "delivery_mode", "rurality", "province", "ref_source",
              "n_cand", "n_eff", "n_da", "dA_exp", "dB", "dC", "dB_planar", "dC_planar", "dC_exp_pccf", "A_w_refda", "same_da_B", "same_da_C",
              "eD", "eA", "eB", "eC", "evi_ok", "dAC", "dBC", "dAB", "A_V_evi_b1_mean", "A_miss_evi_b1_mean",
              "ref_QuintMat", "A_q_mat", "B_q_mat", "pccf_QuintMat", "ref_QuintSoc", "A_q_soc", "B_q_soc", "pccf_QuintSoc",
              "ref_SCOREMAT", "A_S_mat", "B_S_mat", "pccf_SCOREMAT", "ref_SCORESOC", "A_S_soc", "B_S_soc", "pccf_SCORESOC",
              "q_ok_mat", "q_ok_soc", "C_expmis_mat", "C_expmis_soc", "A_nq_mat", "A_pimax_mat")
saveRDS(x[, keep_out, with = FALSE], file.path(OUT, "address_outcomes.rds"))
msg("Address outcomes saved:", fmt(nrow(x)), "rows | EVI paired:", fmt(x[, sum(evi_ok)]),
    "| material paired:", fmt(x[, sum(q_ok_mat)]), "| social paired:", fmt(x[, sum(q_ok_soc)]))

# compact analysis tables (avoid copying the full address table)
STRAT <- c("Link_Source", "delivery_mode", "rurality", "province")
ev <- x[(evi_ok), c("pidx", STRAT, "eD", "eA", "eB", "eC", "dAC", "dBC", "dAB"), with = FALSE]
mq <- lapply(c(mat = "mat", soc = "soc"), function(k) {
  K <- tools::toTitleCase(k); U <- toupper(k)
  x[get(paste0("q_ok_", k)) == TRUE, c(list(pidx = pidx), mget(STRAT), list(
    t = get(paste0("ref_Quint", K)), qA = get(paste0("A_q_", k)), qB = get(paste0("B_q_", k)), qC = get(paste0("pccf_Quint", K)),
    sD = get(paste0("ref_SCORE", U)), sA = get(paste0("A_S_", k)), sB = get(paste0("B_S_", k)), sC = get(paste0("pccf_SCORE", U)),
    expmis = get(paste0("C_expmis_", k))))]
})

# ---- 4. Descriptive tables ------------------------------------------------------------------------------
groupings <- list(overall = character(0), Link_Source = "Link_Source", delivery_mode = "delivery_mode",
                  rurality = "rurality", dm_x_rurality = c("delivery_mode", "rurality"), province = "province")
by_group <- function(fun) rbindlist(lapply(names(groupings), function(g) {
  r <- fun(groupings[[g]]); if (!nrow(r)) return(NULL)
  r[, grouping := g]; setcolorder(r, "grouping"); r }), fill = TRUE)
r4 <- function(v) round(v, 4)

t4 <- by_group(function(by) x[, .(addresses = .N, postal_codes = uniqueN(PCODE), mean_n_cand = round(mean(n_cand), 2),
                                  mean_n_eff = round(mean(n_eff), 2), mean_n_da = round(mean(n_da), 2),
                                  mean_sd_cand_evi = r4(mean(sqrt(A_V_evi_b1_mean), na.rm = TRUE)),
                                  pct_multi_quint_mat = round(100 * mean(A_nq_mat > 1, na.rm = TRUE), 1),
                                  mean_pimax_mat = round(mean(A_pimax_mat, na.rm = TRUE), 3)), keyby = by])
ps <- function(v, p) setNames(list(round(median(v)), round(mean(v)), round(quantile(v, .9)), round(100 * mean(v <= 100), 1),
                                   round(100 * mean(v > 5000), 1)),
                              paste0(p, c("median_m", "mean_m", "p90_m", "pct_le_100m", "pct_gt_5km")))
t5 <- by_group(function(by) x[, c(list(N = .N), ps(dA_exp, "A_exp_"), ps(dB, "B_"), ps(dC, "C_"),
  list(C_expected_planar_mean_m = round(mean(dC_exp_pccf)), C_realized_planar_mean_m = round(mean(dC_planar)),
       A_pct_weight_ref_da = round(100 * mean(A_w_refda), 1), B_pct_same_da = round(100 * mean(same_da_B), 1),
       C_pct_same_da = round(100 * mean(same_da_C), 1))), keyby = by])
t6 <- by_group(function(by) rbindlist(list(
  ev[, .(method = "A", N = .N, bias = r4(mean(eA - eD)), mae = r4(mean(abs(eA - eD))), rmse = r4(sqrt(mean((eA - eD)^2))),
         pearson_r = round(cor(eA, eD), 3), mean_delta_vs_C = r4(mean(dAC)), mean_delta_A_vs_B = r4(mean(dAB))), keyby = by],
  ev[, .(method = "B", N = .N, bias = r4(mean(eB - eD)), mae = r4(mean(abs(eB - eD))), rmse = r4(sqrt(mean((eB - eD)^2))),
         pearson_r = round(cor(eB, eD), 3), mean_delta_vs_C = r4(mean(dBC)), mean_delta_A_vs_B = r4(mean(dAB))), keyby = by],
  ev[, .(method = "C", N = .N, bias = r4(mean(eC - eD)), mae = r4(mean(abs(eC - eD))), rmse = r4(sqrt(mean((eC - eD)^2))),
         pearson_r = round(cor(eC, eD), 3), mean_delta_vs_C = NA_real_, mean_delta_A_vs_B = r4(mean(dAB))), keyby = by])))
msdi_stats <- function(d, by, M, k) {
  q <- paste0("q", M); s <- paste0("s", M)
  d[, .(dimension = k, method = M, N = .N,
        pct_misclassified = round(100 * mean(get(q) != t), 2), pct_more_deprived = round(100 * mean(get(q) > t), 2),
        pct_less_deprived = round(100 * mean(get(q) < t), 2), pct_off_2plus = round(100 * mean(abs(get(q) - t) >= 2), 2),
        kappa_linear = round(wkappa(t, get(q), "linear"), 4),
        pct_assigned_q1 = round(100 * mean(get(q) == 1), 2), pct_assigned_q5 = round(100 * mean(get(q) == 5), 2),
        pct_ref_q1 = round(100 * mean(t == 1), 2), pct_ref_q5 = round(100 * mean(t == 5), 2),
        C_expected_misclassified = if (M == "C") round(100 * mean(expmis), 2) else NA_real_,
        score_bias = signif(mean(get(s) - sD, na.rm = TRUE), 4), score_mae = signif(mean(abs(get(s) - sD), na.rm = TRUE), 4),
        score_rmse = signif(sqrt(mean((get(s) - sD)^2, na.rm = TRUE)), 4)), keyby = by]
}
t7 <- by_group(function(by) rbindlist(lapply(c("mat", "soc"), function(k)
  rbindlist(lapply(c("A", "B", "C"), function(M) msdi_stats(mq[[k]], by, M, k))))))
qshare <- rbindlist(lapply(c("mat", "soc"), function(k) rbindlist(lapply(c(D = "t", A = "qA", B = "qB", C = "qC"), function(v)
  mq[[k]][, .(N = .N), keyby = .(Link_Source, q = get(v))][, pct := round(100 * N / sum(N), 2), by = Link_Source]),
  idcol = "method")[, dimension := k]))
msg("Descriptive tables done")

# ---- 5. Bayesian bootstrap over postal codes ---------------------------------------------------------------
# Sums are aggregated to the finest cell (postal code x stratum variables); each replicate draws
# Exp(1) weights per postal code (Rubin 1981) and re-aggregates by stratum. The same function with
# unit weights gives the point estimates.
bgroups <- c("overall", "Link_Source", "delivery_mode", "rurality", "province")
fine <- c("pidx", STRAT)
fe <- ev[, .(n = .N, sAC = sum(dAC), sBC = sum(dBC), sAB = sum(dAB)), by = fine]
fm <- lapply(mq, function(d) d[, .(n = .N, mA = sum(qA != t), mB = sum(qB != t), mC = sum(qC != t)), by = fine])
fc <- lapply(mq, function(d) rbindlist(list(
  A = d[, .(cnt = .N), by = c(fine, "t", "qA")][, a := qA][, qA := NULL],
  B = d[, .(cnt = .N), by = c(fine, "t", "qB")][, a := qB][, qB := NULL],
  C = d[, .(cnt = .N), by = c(fine, "t", "qC")][, a := qC][, qC := NULL]), idcol = "method"))
gid_of  <- function(d, g) if (g == "overall") rep(1L, nrow(d)) else as.integer(d[[g]])
glab_of <- function(d, g) if (g == "overall") "all" else levels(d[[g]])
kappa_counts <- function(cnt, t, a) {    # linear weighted kappa from (weighted) counts
  m <- matrix(0, 5, 5); m[cbind(t, a)] <- cnt; n <- sum(m)
  w <- 1 - abs(outer(1:5, 1:5, "-")) / 4
  pe <- sum(w * outer(rowSums(m), colSums(m))) / n^2
  (sum(w * m) / n - pe) / (1 - pe)
}
# decomposition sums under PCCF+ weights (Link_Source C, 5, 4 checked; D reported)
x[, `:=`(dec_evi_ok = evi_ok & Aw_miss_evi_b1_mean < 1e-12,
         dec_S_ok = !is.na(ref_SCOREMAT) & !is.na(pccf_SCOREMAT) & !is.na(Aw_S_mat) & Aw_miss_S_mat < 1e-12)]
fd <- rbind(
  x[(dec_evi_ok), .(measure = "evi_b1", n = .N, seC = sum((eC - eD)^2), seA = sum((Aw_evi_b1_mean - eD)^2), sV = sum(Aw_V_evi_b1_mean)), by = .(pidx, Link_Source)],
  x[(dec_S_ok), .(measure = "score_mat", n = .N, seC = sum((pccf_SCOREMAT - ref_SCOREMAT)^2), seA = sum((Aw_S_mat - ref_SCOREMAT)^2), sV = sum(Aw_V_S_mat)), by = .(pidx, Link_Source)])
fd[, dgrp := fifelse(Link_Source %in% c("C", "5", "4"), "C/5/4", fifelse(Link_Source == "D", "D (PCCF+ weights)", "F"))]
fd[, dkey := paste(measure, dgrp, sep = ".")]

quantities <- function(w) {
  out <- list()
  for (g in bgroups) {
    se <- rowsum(cbind(fe$n, fe$sAC, fe$sBC, fe$sAB) * w[fe$pidx], gid_of(fe, g))
    lab <- glab_of(fe, g)[as.integer(rownames(se))]
    out[[length(out) + 1]] <- data.table(grouping = g, group = rep(lab, 3), dim = "evi",
      q = rep(c("delta_AC", "delta_BC", "delta_AB"), each = nrow(se)), value = as.vector(se[, 2:4] / se[, 1]))
    for (k in c("mat", "soc")) {
      f <- fm[[k]]
      sm <- rowsum(cbind(f$n, f$mA, f$mB, f$mC) * w[f$pidx], gid_of(f, g))
      rate <- sm[, 2:4, drop = FALSE] / sm[, 1]; lab <- glab_of(f, g)[as.integer(rownames(sm))]
      out[[length(out) + 1]] <- data.table(grouping = g, group = rep(lab, 3), dim = k,
        q = rep(c("misclass_A_minus_C", "misclass_B_minus_C", "misclass_A_minus_B"), each = nrow(sm)),
        value = 100 * c(rate[, 1] - rate[, 3], rate[, 2] - rate[, 3], rate[, 1] - rate[, 2]))
      cc <- fc[[k]]
      kk <- data.table(g = gid_of(cc, g), method = cc$method, t = cc$t, a = cc$a, wc = cc$cnt * w[cc$pidx])[
        , .(wc = sum(wc)), by = .(g, method, t, a)]
      ka <- dcast(kk[, .(kappa = kappa_counts(wc, t, a)), by = .(g, method)], g ~ method, value.var = "kappa")
      lab <- glab_of(cc, g)[ka$g]
      out[[length(out) + 1]] <- data.table(grouping = g, group = rep(lab, 3), dim = k,
        q = rep(c("kappa_A_minus_C", "kappa_B_minus_C", "kappa_A_minus_B"), each = nrow(ka)),
        value = c(ka$A - ka$C, ka$B - ka$C, ka$A - ka$B))
    }
  }
  sd <- rowsum(cbind(fd$n, fd$seC, fd$seA, fd$sV) * w[fd$pidx], fd$dkey)
  out[[length(out) + 1]] <- data.table(grouping = "decomposition", group = rownames(sd), dim = "dec",
    q = "mseC_minus_mseA_plus_V", value = (sd[, 2] - sd[, 3] - sd[, 4]) / sd[, 1])
  rbindlist(out)
}
pt <- quantities(rep(1, npc))
t0 <- Sys.time()
boot <- rbindlist(lapply(seq_len(REPS), function(b) {
  r <- quantities(rexp(npc))[, b := b]
  if (b %% 20 == 0) msg("bootstrap", b, "of", REPS, "|", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
  r
}))
bs <- boot[, .(boot_sd = sd(value), lo = quantile(value, .025, na.rm = TRUE), hi = quantile(value, .975, na.rm = TRUE)),
           by = .(grouping, group, dim, q)]
tb <- merge(pt[, .(grouping, group, dim, q, estimate = value)], bs, by = c("grouping", "group", "dim", "q"), all.x = TRUE)
msg("Bootstrap done:", REPS, "replicates in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

# ---- 6. MSE decomposition check (Eq. 7) --------------------------------------------------------------------
dec <- fd[, .(N = sum(n), mse_C = sum(seC) / sum(n), mse_A = sum(seA) / sum(n), mean_V = sum(sV) / sum(n)), by = .(measure, dgrp)]
dec[, `:=`(diff = mse_C - mse_A - mean_V, group = paste(measure, dgrp, sep = "."))]
dec[bs[grouping == "decomposition"], boot_se := i.boot_sd, on = "group"]
dec[, `:=`(rel_diff_pct = 100 * diff / mse_C, z = diff / boot_se, share_of_mse_C_from_V = mean_V / mse_C)]
print(dec)
for (i in which(dec$dgrp %in% c("C/5/4", "D (PCCF+ weights)")))
  check(abs(dec$z[i]) <= 3.29, sprintf("MSE decomposition, Link_Source %s, %s: MSE(C) %.6g vs MSE(A) + mean V %.6g (diff %.3g%%, z = %.2f)",
        dec$dgrp[i], dec$measure[i], dec$mse_C[i], dec$mse_A[i] + dec$mean_V[i], dec$rel_diff_pct[i], dec$z[i]))

# ---- 7. Cell-level models (validation design) ----------------------------------------------------------------
coefs <- list(); tests <- list(); preds <- list()
# OLS for Delta_AC: delivery mode x rurality + province
cells <- model_cells(ev, "dAC")
fmf <- model_formulas(cells, "ybar")
full <- ols_cells(cells, fmf$full); red <- ols_cells(cells, fmf$main)
est <- coef(full$fit); keep <- !is.na(est)
coefs$delta <- coef_table(est[keep], full$V, "delta_AC_evi_b1", "main")
df1 <- red$df - full$df; Fst <- ((red$rss - full$rss) / df1) / full$s2
tests$delta <- data.table(model = "delta_AC_evi_b1", test = "delivery mode x rurality", statistic = "F", value = Fst,
                          df1 = df1, df2 = full$df, p = pf(Fst, df1, full$df, lower.tail = FALSE), N = sum(cells$n))
preds$delta <- cbind(model = "delta_AC_evi_b1", adj_pred(full$fit, cells))
# logistic for material misclassification: method x delivery mode x rurality + province
mc <- drop_sparse(mq$mat)[, .(n = .N, A = sum(qA != t), B = sum(qB != t), C = sum(qC != t)), by = .(delivery_mode, rurality, province)]
mc <- melt(mc, id.vars = c("delivery_mode", "rurality", "province", "n"), variable.name = "method", value.name = "k")
mc[, `:=`(delivery_mode = droplevels(delivery_mode), rurality = droplevels(rurality), province = droplevels(province),
          method = factor(method, c("C", "A", "B")))]
# terms with a single level in the data are dropped, as in model_formulas() of the validation study
has <- function(v) nlevels(mc[[v]]) > 1L
core <- c("delivery_mode", "rurality")[c(has("delivery_mode"), has("rurality"))]
adj  <- if (has("province")) " + province" else ""
f_full <- as.formula(paste("cbind(k, n - k) ~", paste(c("method", core), collapse = " * "), adj))
f_red  <- as.formula(paste("cbind(k, n - k) ~ method +", if (length(core)) paste(core, collapse = " * ") else "1", adj))
gfull <- glm(f_full, family = binomial, data = mc); gred <- glm(f_red, family = binomial, data = mc)
est <- coef(gfull); keep <- !is.na(est)
ct <- coef_table(est[keep], vcov(gfull, complete = FALSE), "misclass_mat_by_method", "main", scale = "log-odds")
ct[, `:=`(odds_ratio = exp(estimate), or_low = exp(ci_low), or_high = exp(ci_high))]
coefs$misclass <- ct
lr <- gred$deviance - gfull$deviance; df1 <- gred$df.residual - gfull$df.residual
tests$misclass <- data.table(model = "misclass_mat_by_method", test = "method x (delivery mode, rurality) interactions",
                             statistic = "LR chi-square", value = lr, df1 = df1, df2 = NA_real_,
                             p = pchisq(lr, df1, lower.tail = FALSE), N = sum(mc$n) / 3)
pw <- mc[method == "C", .(w = sum(n)), by = province][, w := w / sum(w)]
grid <- CJ(method = levels(mc$method), delivery_mode = levels(mc$delivery_mode), rurality = levels(mc$rurality), province = pw$province)
for (v in c("method", "delivery_mode", "rurality", "province")) grid[, (v) := factor(get(v), levels(mc[[v]]))]
grid[, pred := suppressWarnings(predict(gfull, newdata = grid, type = "response"))]
grid <- merge(grid, pw, by = "province")
ap <- grid[, .(adjusted_prediction = sum(pred * w)), by = .(method, delivery_mode, rurality)]
ap <- merge(ap, mc[method == "C", .(N = sum(n)), by = .(delivery_mode, rurality)], by = c("delivery_mode", "rurality"), all.x = TRUE)
ap[is.na(N), `:=`(N = 0L, adjusted_prediction = NA_real_)]
preds$misclass <- cbind(model = "misclass_mat_by_method", ap)
msg("Models done")

# ---- 8. Write -------------------------------------------------------------------------------------------------
tabs <- list(t4_candidates = t4, t5_position = t5, t6_evi = t6, t7_msdi = t7, t7_quintile_shares = qshare,
             bootstrap = tb, decomposition = dec, model_coefficients = rbindlist(coefs, fill = TRUE),
             model_tests = rbindlist(tests, fill = TRUE), model_adjusted_predictions = rbindlist(preds, fill = TRUE))
for (n in names(tabs)) fwrite(tabs[[n]], file.path(OUT, paste0(n, ".csv")))
write_xlsx(lapply(tabs, as.data.frame), file.path(OUT, "sw4_results.xlsx"))
fwrite(boot, file.path(OUT, "bootstrap_replicates.csv.gz"))
o <- t6[grouping == "overall"]; m7 <- t7[grouping == "overall" & dimension == "mat"]
msg(sprintf("Overall EVI MAE: A %.4f, B %.4f, C %.4f | mean Delta_AC %.4f | material misclassified: A %.1f%%, B %.1f%%, C %.1f%%",
            o[method == "A", mae], o[method == "B", mae], o[method == "C", mae], o[method == "A", mean_delta_vs_C],
            m7[method == "A", pct_misclassified], m7[method == "B", pct_misclassified], m7[method == "C", pct_misclassified]))
msg("SW4 done")
