# =============================================================================
# SW8  Manuscript tables as Markdown pipe tables for manuscript/manuscript_src.md
#
# Reads the step 6 tables (data/derived/sw6/tables/Table4-7) and writes
# data/derived/sw8/manuscript_tables.md with the four results tables of the preprint,
# numbered as in the manuscript (Tables 3-6; Tables 1-2 are the methods tables). The Markdown is
# pasted into manuscript_src.md, so every number in the tables comes from the analysis outputs.
#
# Run: Rscript scripts/sw8_manuscript_tables.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
TAB  <- file.path(DERIVED, "sw6", "tables")
OUT  <- file.path(DERIVED, "sw8"); dir.create(OUT, showWarnings = FALSE)

n0 <- function(v) formatC(v, format = "d", big.mark = ",")
f  <- function(v, d) {       # fixed decimals, true minus sign, no sign on values that round to zero
  s <- formatC(abs(v), format = "f", digits = d)
  s <- ifelse(!is.na(v) & v < 0 & grepl("[1-9]", s), paste0("−", s), s); s[is.na(v)] <- "—"; s
}
lab <- function(block, level) fifelse(block == "All", "All addresses", level)
# rows of one table with a panel heading row before each block (Link source, Delivery mode, Rurality)
with_panels <- function(d, k, row_fun) {
  rows <- list(); prev <- ""
  for (i in seq_len(nrow(d))) {
    r <- d[i]
    if (r$block != "All" && r$block != prev) rows[[length(rows) + 1]] <- panel(r$block, k)
    prev <- r$block
    rows[[length(rows) + 1]] <- row_fun(r)
  }
  rows
}
ci  <- function(e, lo, hi, d) sprintf("%s (%s, %s)", f(e, d), f(lo, d), f(hi, d))
md  <- function(header, rows) c(paste0("| ", paste(header, collapse = " | "), " |"),
                               paste0("|", paste(rep(" --- ", length(header)), collapse = "|"), "|"),
                               vapply(rows, function(r) paste0("| ", paste(r, collapse = " | "), " |"), ""))
panel <- function(name, k) c(paste0("*", name, "*"), rep("", k - 1))


out <- character()

# ---- Table 3: candidate locations ------------------------------------------------------------------
t4 <- fread(file.path(TAB, "Table4_candidates.csv"))
rows <- with_panels(t4, 10, function(r) with(r, c(lab(block, level), n0(postal_codes), n0(addresses), f(pct_addresses, 1),
  f(mean_candidates, 2), f(mean_effective_candidates, 2), f(pct_multi_da, 1), f(mean_sd_candidate_evi, 3),
  f(pct_multi_quintile_mat, 1), f(mean_max_pi_mat, 2))))
out <- c(out, "<!-- Table 3 -->", md(c("Group", "Postal codes", "Addresses", "%", "Cand.", "Eff. cand.",
                                       "Multi-DA (%)", "SD EVI", "Multi-quintile (%)", "Largest π"), rows), "")

# ---- Table 4: positional error -----------------------------------------------------------------------
t5 <- fread(file.path(TAB, "Table5_position.csv"))
rows <- with_panels(t5, 11, function(r) with(r, c(lab(block, level), n0(N), n0(A_expected_median_m), n0(B_median_m), n0(C_median_m),
  n0(A_expected_mean_m), n0(B_mean_m), n0(C_mean_m), f(A_pct_weight_in_ref_da, 1), f(B_pct_same_da, 1), f(C_pct_same_da, 1))))
out <- c(out, "<!-- Table 4 -->", md(c("Group", "Addresses", "Median A", "Median B", "Median C", "Mean A", "Mean B", "Mean C",
                                       "Same DA A", "Same DA B", "Same DA C"), rows), "")

# ---- Table 5: greenness ---------------------------------------------------------------------------------
t6 <- fread(file.path(TAB, "Table6_greenness.csv"))
rows <- with_panels(t6, 10, function(r) with(r, c(lab(block, level), n0(N), f(bias_A, 4), f(bias_C, 4), f(mae_A, 4), f(mae_B, 4),
  f(mae_C, 4), f(r_A, 3), f(r_C, 3), ci(delta_AC, delta_AC_lo, delta_AC_hi, 4))))
out <- c(out, "<!-- Table 5 -->", md(c("Group", "Addresses", "Bias A", "Bias C", "MAE A", "MAE B", "MAE C", "r A", "r C",
                                       "Δ~AC~ (95% interval)"), rows), "")

# ---- Table 6: deprivation ---------------------------------------------------------------------------------
t7 <- fread(file.path(TAB, "Table7_deprivation.csv"))
rows <- list()
for (dm in c("Material", "Social")) {
  rows[[length(rows) + 1]] <- panel(paste(dm, "deprivation"), 11)
  rows <- c(rows, with_panels(t7[dimension == dm], 11, function(r) with(r, c(lab(block, level), n0(N),
    f(misclassified_A_mean, 1), f(misclassified_A_mode, 1), f(misclassified_B, 1), f(misclassified_C, 1),
    ci(diff_A_mode_minus_C, diff_lo, diff_hi, 1), f(off2_A_mode, 1), f(off2_C, 1), f(kappa_A_mode, 3), f(kappa_C, 3)))))
}
out <- c(out, "<!-- Table 6 -->", md(c("Group", "Addresses", "A mean", "A modal", "B", "C",
                                       "A modal − C (95% interval)", "≥2 off A modal", "≥2 off C", "κ A modal", "κ C"), rows), "")

writeLines(out, file.path(OUT, "manuscript_tables.md"), useBytes = FALSE)
cat("Wrote", file.path(OUT, "manuscript_tables.md"), "\n")
