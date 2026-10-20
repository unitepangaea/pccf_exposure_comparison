# =============================================================================
# SW6  Manuscript tables, figures and supplementary tables
# Spatial-weight protocol, Section 3 (planned results)
#
# Inputs: data/derived/sw1-sw5 (primary weights: PCCF+ selection probabilities,
#         equal weights for Link_Source D; amendment of 14 October 2026), VAL_DERIVED/linked,
#         VAL_DERIVED/p2/p2_analysis.rds, VAL_DERIVED/p4/tables/Table1_sample.csv (validation study)
# Outputs (data/derived/sw6/):
#   tables/Table4_candidates.csv     candidates and weights by link source and delivery mode
#   tables/Table5_position.csv       positional error and DA agreement, methods A, B, C
#   tables/Table6_greenness.csv      greenness agreement and paired difference Delta_AC
#   tables/Table7_deprivation.csv    quintile agreement, both rules for method A
#   tables/manuscript_tables.xlsx, tables/supplementary_tables.xlsx (S1-S12)
#   figures/fig1_example_postal_code, fig2_delta_dm_rurality, fig3_decomposition,
#   figures/fig4_confusion_rural  (.png 300 dpi, .pdf, _ms.png without title)
# Check (stop if it fails): for every reference source, delivery mode and rurality category,
#   analysed addresses + addresses excluded for Link_Source 3 and 2 = validation Table 1
#
# Run: Rscript scripts/sw6_figures_tables.R
# Log: data/derived/sw6_figures_tables.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales)
                                 library(writexl); library(sf) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW   <- DERIVED
OUT  <- file.path(SW, "sw6"); FIG <- file.path(OUT, "figures"); TAB <- file.path(OUT, "tables")
LOG  <- file.path(SW, "sw6_figures_tables.log")
PRIMARY <- "pccf"
for (d in c(FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

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
# style of the validation study (P4): palette, theme and figure saving
use("p4_figures_tables.R", c("ink", "seq_blue", "seq_orange", "theme_paper", "save_fig", "dm_short", "rur_levels"))
use("p3_analysis.R", c("dm_levels"))
msg("SW6 start")

rd <- function(step, f) fread(file.path(SW, step, f))
t5 <- rd("sw4", "t5_position.csv"); t6 <- rd("sw4", "t6_evi.csv"); bt4 <- rd("sw4", "bootstrap.csv")
mm <- rd("sw5", "m_modal.csv"); mb <- rd("sw5", "m_modal_bootstrap.csv"); s6 <- rd("sw5", "s6_decomposition.csv")
LS <- c("F", "C", "D", "5", "4")
# fread returns "" (not NA) for empty character cells, so blank labels are set to NA before coalescing
lab_of <- function(d) {
  nz <- function(v) { v <- as.character(v); v[!is.na(v) & v == ""] <- NA_character_; v }
  d[, level := fcoalesce(nz(Link_Source), nz(delivery_mode), nz(rurality), nz(province), "All addresses")]
}

# ---- 1. Address outcomes and the check against validation Table 1 ---------------------------------------
x <- readRDS(file.path(SW, "sw4", "address_outcomes.rds"))
excl <- readRDS(file.path(VAL_DERIVED, "p2", "p2_analysis.rds"))[pccf_coded == TRUE & pccf_Link_Source %chin% c("3", "2", "1", "0"),
          .(ref_source, delivery_mode, rurality_sac)]
excl[, urban_lc := delivery_mode %in% dm_levels[1:3]]
excl <- excl[!(ref_source == "blockface" & urban_lc)]
excl[, rurality := fcase(rurality_sac %in% c("CA (tracted)", "CA (non-tracted)"), "CA",
                         rurality_sac == "Territories (outside CA)", "Territories", default = as.character(rurality_sac))]
t1 <- fread(file.path(VAL_DERIVED, "p4", "tables", "Table1_sample.csv"))
cmp <- rbind(
  merge(x[, .N, by = .(level = ref_source)], excl[, .(nx = .N), by = .(level = ref_source)], by = "level", all = TRUE)[, group := "Reference source"],
  merge(x[, .N, by = .(level = as.character(delivery_mode))], excl[, .(nx = .N), by = .(level = delivery_mode)], by = "level", all = TRUE)[, group := "Delivery mode"],
  merge(x[, .N, by = .(level = as.character(rurality))], excl[, .(nx = .N), by = .(level = rurality)], by = "level", all = TRUE)[, group := "Rurality (SAC)"])
cmp[is.na(N), N := 0L][is.na(nx), nx := 0L]
cmp <- merge(cmp, t1[, .(group, level, N_validation = N)], by = c("group", "level"), all = TRUE)
fwrite(cmp, file.path(TAB, "check_validation_table1.csv"))
check(cmp[, all(!is.na(N_validation) & N + nx == N_validation)] && x[, .N] + excl[, .N] == 12412300L,
      sprintf("analysed (%s) + Link_Source 3/2 excluded (%s) = validation Table 1 (12,412,300) in every category",
              fmt(nrow(x)), fmt(nrow(excl))))
rm(excl); invisible(gc())

# modal quintiles (primary weights)
mq <- readRDS(file.path(SW, "sw5", "modal_quintiles.rds"))[wset == PRIMARY, .(PCODE, A_mode_mat, A_mode_soc)]
x[mq, `:=`(A_mode_mat = i.A_mode_mat, A_mode_soc = i.A_mode_soc), on = "PCODE"]; rm(mq)

# ---- 2. Table 4: candidates and weights ----------------------------------------------------------------------
x[, n_eff_p := if (PRIMARY == "pccf") fifelse(Link_Source == "D", as.numeric(n_cand), n_eff) else n_eff]
cand_tab <- function(by) x[, .(postal_codes = uniqueN(PCODE), addresses = .N,
  mean_candidates = round(mean(n_cand), 2), mean_effective_candidates = round(mean(n_eff_p), 2),
  pct_multi_da = round(100 * mean(n_da > 1), 1), mean_sd_candidate_evi = round(mean(sqrt(A_V_evi_b1_mean), na.rm = TRUE), 4),
  pct_multi_quintile_mat = round(100 * mean(A_nq_mat > 1, na.rm = TRUE), 1), mean_max_pi_mat = round(mean(A_pimax_mat, na.rm = TRUE), 3)), keyby = by]
tab4 <- rbind(cand_tab(character(0))[, `:=`(block = "All", level = "All addresses")],
              cand_tab("Link_Source")[, `:=`(block = "Link source", level = as.character(Link_Source))][, Link_Source := NULL],
              cand_tab("delivery_mode")[, `:=`(block = "Delivery mode", level = dm_short[as.character(delivery_mode)])][, delivery_mode := NULL])
tab4[, pct_addresses := round(100 * addresses / x[, .N], 1)]
setcolorder(tab4, c("block", "level", "postal_codes", "addresses", "pct_addresses"))
# Note: postal codes and addresses here are addresses analysed (main sample); multiple-candidate
# statistics are address-weighted.

# ---- 3. Table 5: positional error ----------------------------------------------------------------------------
p5 <- lab_of(copy(t5)[grouping %in% c("overall", "delivery_mode", "rurality")])
check(!anyDuplicated(p5, by = c("grouping", "level")) && all(p5$level != ""), "Table 5: one row per grouping x level")
tab5 <- p5[, .(block = fcase(grouping == "overall", "All", grouping == "delivery_mode", "Delivery mode", default = "Rurality"),
               level = fifelse(grouping == "delivery_mode", dm_short[level], level), N,
               A_expected_median_m = A_exp_median_m, B_median_m, C_median_m,
               A_expected_mean_m = A_exp_mean_m, B_mean_m, C_mean_m,
               A_pct_weight_in_ref_da = A_pct_weight_ref_da, B_pct_same_da, C_pct_same_da)]

# ---- 4. Table 6: greenness -------------------------------------------------------------------------------------
g6 <- lab_of(copy(t6)[grouping %in% c("overall", "Link_Source", "delivery_mode")])
check(!anyDuplicated(g6, by = c("grouping", "level", "method")) && !anyNA(g6$level) && all(g6$level != ""),
      "Table 6: one row per grouping x level x method")
w6 <- dcast(g6, grouping + level ~ method, value.var = c("N", "bias", "mae", "rmse", "pearson_r"))
d6 <- bt4[dim == "evi" & q == "delta_AC" & grouping %in% c("overall", "Link_Source", "delivery_mode")]
d6[, level := fifelse(group == "all", "All addresses", group)]
w6 <- merge(w6, d6[, .(grouping, level, delta_AC = round(estimate, 4), delta_AC_lo = round(lo, 4), delta_AC_hi = round(hi, 4))],
            by = c("grouping", "level"), all.x = TRUE)
tab6 <- w6[, .(block = fcase(grouping == "overall", "All", grouping == "Link_Source", "Link source", default = "Delivery mode"),
               level = fifelse(grouping == "delivery_mode", dm_short[level], level), N = N_A,
               bias_A = bias_A, bias_B = bias_B, bias_C = bias_C, mae_A, mae_B, mae_C, rmse_A, rmse_B, rmse_C,
               r_A = pearson_r_A, r_B = pearson_r_B, r_C = pearson_r_C, delta_AC, delta_AC_lo, delta_AC_hi)]
ord <- function(d) d[order(match(block, c("All", "Link source", "Delivery mode", "Rurality")),
                           match(level, c("All addresses", LS, dm_short, rur_levels)))]
tab4 <- ord(tab4); tab5 <- ord(tab5); tab6 <- ord(tab6)

# ---- 5. Table 7: deprivation -----------------------------------------------------------------------------------
m7 <- lab_of(copy(mm)[grouping %in% c("overall", "Link_Source", "delivery_mode")])
check(!anyDuplicated(m7, by = c("dimension", "grouping", "level", "method")) && all(m7$level != ""),
      "Table 7: one row per dimension x grouping x level x method")
w7 <- dcast(m7, dimension + grouping + level ~ method, value.var = c("N", "pct_misclassified", "pct_off_2plus", "kappa_linear",
                                                                     "pct_assigned_q1", "pct_assigned_q5", "pct_ref_q1", "pct_ref_q5"))
b7 <- mb[comparison == "A_mode minus C" & grouping %in% c("overall", "Link_Source", "delivery_mode")]
b7[, level := fifelse(group == "all", "All addresses", group)]
w7 <- merge(w7, b7[, .(dimension, grouping, level, diff_mode_C = round(misclass_diff_pct, 2), diff_lo = round(mis_lo, 2),
                       diff_hi = round(mis_hi, 2))], by = c("dimension", "grouping", "level"), all.x = TRUE)
tab7 <- w7[, .(dimension = fifelse(dimension == "mat", "Material", "Social"),
               block = fcase(grouping == "overall", "All", grouping == "Link_Source", "Link source", default = "Delivery mode"),
               level = fifelse(grouping == "delivery_mode", dm_short[level], level), N = N_C,
               misclassified_A_mean = pct_misclassified_A_mean, misclassified_A_mode = pct_misclassified_A_mode,
               misclassified_B = pct_misclassified_B, misclassified_C = pct_misclassified_C,
               diff_A_mode_minus_C = diff_mode_C, diff_lo, diff_hi,
               off2_A_mean = pct_off_2plus_A_mean, off2_A_mode = pct_off_2plus_A_mode, off2_C = pct_off_2plus_C,
               kappa_A_mean = kappa_linear_A_mean, kappa_A_mode = kappa_linear_A_mode, kappa_B = kappa_linear_B, kappa_C = kappa_linear_C,
               q1_ref = pct_ref_q1_C, q1_A_mean = pct_assigned_q1_A_mean, q1_A_mode = pct_assigned_q1_A_mode, q1_C = pct_assigned_q1_C,
               q5_ref = pct_ref_q5_C, q5_A_mean = pct_assigned_q5_A_mean, q5_A_mode = pct_assigned_q5_A_mode, q5_C = pct_assigned_q5_C)]
tab7 <- tab7[order(dimension, match(block, c("All", "Link source", "Delivery mode")), match(level, c("All addresses", LS, dm_short)))]
msg("Tables 4-7 built")

# ---- 6. Figure 1: example rural postal code -----------------------------------------------------------------------
pv <- readRDS(file.path(SW, "sw3", "postal_values.rds"))[wset == PRIMARY]
xc <- x[Link_Source == "C" & delivery_mode == "Rural postal code", .(n_addr = .N, n_cand = n_cand[1], sd = sqrt(A_V_evi_b1_mean[1]),
                                                                    dAC = mean(dAC, na.rm = TRUE), province = province[1]), by = PCODE]
qs <- xc[sd > 0, quantile(sd, c(.4, .6))]
ex <- xc[n_cand %between% c(4, 8) & n_addr %between% c(150, 400) & sd %between% qs][order(-n_addr)][1]
msg("Figure 1 postal code:", ex$PCODE, "| candidates", ex$n_cand, "| addresses", ex$n_addr, "| SD of candidate EVI", round(ex$sd, 3))
ce <- readRDS(file.path(SW, "sw2", "cand_exposures.rds"))[PCODE == ex$PCODE]
ce[, w := if (PRIMARY == "pccf") fifelse(Link_Source == "D", w_equal, w_pop) else w_pop]
ex_pv <- pv[PCODE == ex$PCODE]
ad <- readRDS(file.path(VAL_DERIVED, "linked", paste0("linked_", ex$province, ".rds")))[PCODE == ex$PCODE & pccf_coded == TRUE,
        .(ADDR_GUID, ref_x, ref_y, pccf_x, pccf_y)]
ad <- ad[x[PCODE == ex$PCODE, .(ADDR_GUID, eD, eC)], on = "ADDR_GUID", nomatch = NULL]
drawn <- ad[, .N, by = .(pccf_x, pccf_y)]
das <- unique(c(ce$DAuid, ex_pv$B_DAuid))
da <- st_read(file.path(RAW, "lda_000b21a_e", "lda_000b21a_e.shp"), quiet = TRUE)
da <- st_transform(da[da$DAUID %in% das, "DAUID"], 3347)
# label offset grows with the symbol radius (scale_size_area: diameter proportional to sqrt(w), max 9 mm)
yr <- diff(st_bbox(da)[c("ymin", "ymax")])
ce[, ly := y + (4.5 * sqrt(w / max(w)) + 2.5) * yr / 120]
f1 <- ggplot() +
  geom_sf(data = da, fill = ink$na, colour = ink$axis, linewidth = 0.3) +
  geom_point(data = ad, aes(ref_x, ref_y, colour = "Address (reference location)"), size = 0.5, alpha = 0.6) +
  geom_point(data = ce, aes(x, y, size = w, fill = evi_b1_mean), shape = 21, colour = ink$primary, stroke = 0.4) +
  geom_point(data = ex_pv, aes(B_x, B_y, colour = "Weighted centroid (method B)"), shape = 4, size = 3.5, stroke = 1.2) +
  geom_text(data = ce, aes(x, ly, label = sprintf("%.2f", evi_b1_mean)), size = 2.6, colour = ink$primary) +
  scale_colour_manual(values = c(`Address (reference location)` = ink$secondary, `Weighted centroid (method B)` = "#eb6834"), name = NULL,
                      guide = guide_legend(override.aes = list(shape = c(16, 4), size = c(1.5, 3), alpha = 1, stroke = c(0.5, 1.2)), order = 3)) +
  scale_fill_gradientn(colours = c("#f5e9c8", "#9fcf8c", "#2f8f4e"), name = "Candidate EVI (3x3)", guide = guide_colourbar(order = 1)) +
  scale_size_area(max_size = 9, name = "Candidate weight", guide = guide_legend(order = 2)) +
  labs(title = paste("Candidate locations of rural postal code", ex$PCODE),
       subtitle = sprintf("Method A (weighted mean) EVI %.3f | method B (centroid, x) EVI %.3f | mean EVI at the %s addresses %.3f",
                          ex_pv$A_evi_b1_mean, ex_pv$B_evi_b1_mean, comma(nrow(ad)), mean(ad$eD, na.rm = TRUE)),
       x = NULL, y = NULL) +
  theme_paper() + theme(axis.text = element_blank(), axis.line.x = element_blank(), axis.ticks = element_blank(),
                        panel.grid.major = element_blank())
save_fig(f1, "fig1_example_postal_code", 7, 6)

# ---- 7. Figure 2: Delta_AC by delivery mode x rurality ------------------------------------------------------------
g2 <- copy(t6)[grouping == "dm_x_rurality" & method == "A" & N >= 100]
g2[, `:=`(dm = factor(dm_short[delivery_mode], rev(dm_short)), rur = factor(rurality, rur_levels))]
lim <- g2[, max(abs(mean_delta_vs_C))]
f2 <- ggplot(g2, aes(rur, dm, fill = mean_delta_vs_C)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sub("^-(0\\.0+)$", "\\1", sprintf("%.3f", mean_delta_vs_C)), colour = abs(mean_delta_vs_C) > lim * 0.55), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink$primary)) +
  scale_fill_gradient2(low = "#1c5cab", mid = "#f7f7f5", high = "#b84a1c", midpoint = 0, limits = c(-lim, lim),
                       name = "Mean paired\ndifference\n(EVI units)") +
  labs(title = "Paired difference in absolute EVI error, weighted mean (A) minus PCCF+ (C)",
       subtitle = "3x3 window, 2025-26; negative values favour method A; cells with at least 100 addresses",
       x = "Rurality of the reference location (SAC)", y = NULL) +
  theme_paper() + theme(panel.grid.major = element_blank())
save_fig(f2, "fig2_delta_dm_rurality", 7.5, 3.6)

# ---- 8. Figure 3: decomposition of the squared error of method C ---------------------------------------------------
g3 <- s6[measure %in% c("evi_b1_mean", "score_mat", "score_soc")]
g3[, meas := factor(c(evi_b1_mean = "Greenness (EVI, 3x3)", score_mat = "Material score", score_soc = "Social score")[measure],
                    c("Greenness (EVI, 3x3)", "Material score", "Social score"))]
g3[, grp := factor(fifelse(grp == "C/5/4", "Link source C, 5, 4", "Link source D"), c("Link source C, 5, 4", "Link source D"))]
g3l <- melt(g3[, .(meas, grp, `Error of the weighted mean (method A)` = 1 - share_of_mse_C_from_V,
                   `Allocation variance (removed by method A)` = share_of_mse_C_from_V)], id.vars = c("meas", "grp"), variable.name = "part")
f3 <- ggplot(g3l, aes(value, grp, fill = part)) +
  geom_col(width = 0.6, colour = "white") +
  geom_text(aes(label = percent(value, accuracy = 1), colour = part), position = position_stack(vjust = 0.5), size = 2.8,
            show.legend = FALSE) +
  scale_colour_manual(values = c(ink$primary, "white")) +
  facet_wrap(~ meas, ncol = 1) +
  scale_fill_manual(values = c("#cde2fb", "#1c5cab"), name = NULL) +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Decomposition of the mean squared error of PCCF+ random allocation (method C)",
       subtitle = "MSE(C) = MSE(A) + mean allocation variance (Equation 7), PCCF+ selection weights", x = "Share of MSE(C)", y = NULL) +
  theme_paper() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_fig(f3, "fig3_decomposition", 6.5, 4.6)

# ---- 9. Figure 4: confusion matrices, rural postal-code addresses --------------------------------------------------
r4 <- x[delivery_mode == "Rural postal code" & q_ok_mat == TRUE & !is.na(A_mode_mat)]
cm <- rbindlist(lapply(c("C. PCCF+ random allocation" = "pccf_QuintMat", "A. Weighted mean score" = "A_q_mat",
                         "A. Modal quintile" = "A_mode_mat"), function(v)
  r4[, .N, by = .(t = ref_QuintMat, a = get(v))][, pct := 100 * N / sum(N), by = t]), idcol = "method")
cm[, method := factor(method, c("C. PCCF+ random allocation", "A. Weighted mean score", "A. Modal quintile"))]
f4 <- ggplot(cm, aes(factor(a), factor(t, 5:1), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = round(pct), colour = pct > 45), size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink$primary)) +
  scale_fill_gradientn(colours = seq_blue, limits = c(0, 100), name = "Row %") +
  facet_wrap(~ method, nrow = 1) +
  labs(title = "Material deprivation quintile, rural postal-code addresses",
       subtitle = sprintf("Rows: quintile at the address; columns: assigned quintile (n = %s)", comma(nrow(r4))),
       x = "Assigned quintile", y = "Quintile at the address") +
  theme_paper() + theme(panel.grid.major = element_blank())
save_fig(f4, "fig4_confusion_rural", 8, 3.4)

# ---- 10. Supplementary tables -------------------------------------------------------------------------------------
cprov <- x[, .(postal_codes = uniqueN(PCODE), addresses = .N, mean_candidates = round(mean(n_cand), 2),
               mean_effective_candidates = round(mean(n_eff_p), 2)), keyby = .(province, Link_Source)]
supp <- list(
  S1_candidates_province_linksource = cprov,
  S2_position_province = t5[grouping == "province"],
  S3_greenness_province = t6[grouping == "province"],
  S4_deprivation_province = mm[grouping == "province"],
  S5_outcomes_link_source = rbind(t5[grouping == "Link_Source"][, table := "position"], t6[grouping == "Link_Source"][, table := "greenness"],
                                  mm[grouping == "Link_Source"][, table := "deprivation"], fill = TRUE),
  S6_sens1_alt_weights_D = rd("sw5", "s1_alt_weights.csv"),
  S6b_sens1_expected_misclassification = rd("sw5", "s1_expected_misclassification.csv"),
  S7_sens2_latlon_centroid = rd("sw5", "s2_latlon_centroid.csv"),
  S8_sens3_evi_windows_years = rd("sw5", "s3_evi_measures.csv"),
  S9_sens4_building_only = rd("sw5", "s4_building_only.csv"),
  S10_sens5_with_urban_blockface = rd("sw5", "s5_with_urban_blockface.csv"),
  S11_sens6_decomposition = s6,
  S12_models = rbind(rd("sw4", "model_coefficients.csv")[, part := "coefficients"], rd("sw4", "model_tests.csv")[, part := "tests"],
                     rd("sw4", "model_adjusted_predictions.csv")[, part := "adjusted predictions"], fill = TRUE))
readme <- data.table(sheet = names(supp), description = c(
  "Postal codes, addresses and mean (effective) number of candidates by province and PCCF+ link source",
  "Positional error and DA agreement by province (methods A expected distance, B, C)",
  "Greenness agreement by province (EVI 3x3, 2025-26), methods A, B, C",
  "Material and social quintile agreement by province, methods A (mean score, modal), B, C",
  "Positional error, greenness and deprivation outcomes by PCCF+ link source",
  "Sensitivity 1: DA population-share weights for Link_Source D (original protocol primary weighting)",
  "Sensitivity 1 check: expected misclassification of random allocation (1 - pi) vs observed, by link source",
  "Sensitivity 2: method B centroid from latitude/longitude",
  "Sensitivity 3: EVI windows (pixel, 3x3, 5x5) and years",
  "Sensitivity 4: addresses with a NAR building reference point",
  "Sensitivity 5: including urban letter-carrier addresses with a blockface reference",
  "Sensitivity 6: MSE decomposition for all EVI measures and MSDI scores",
  "Cell-level models: OLS for Delta_AC; logistic for material misclassification by method"))
mt <- list(Table4_candidates = tab4, Table5_position = tab5, Table6_greenness = tab6, Table7_deprivation = tab7)
for (n in names(mt)) fwrite(mt[[n]], file.path(TAB, paste0(n, ".csv")))
write_xlsx(lapply(mt, as.data.frame), file.path(TAB, "manuscript_tables.xlsx"))
write_xlsx(c(list(README = as.data.frame(readme)), lapply(supp, as.data.frame)), file.path(TAB, "supplementary_tables.xlsx"))
msg("Tables written:", TAB)
msg("SW6 done")
