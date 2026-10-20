# =============================================================================
# Runs the analysis in order. Each step is a separate R process, so a failure stops the
# run at that step and later steps can be rerun individually. Every script ends with
# checks that stop it if they fail. Requires the pccf_nar_validation pipeline to have been
# run (see README.md). Approximate times on a 2024 laptop (16 GB): sw1 7 min, sw2 8 min,
# sw3 7 min, sw4 21 min, sw5 65 min, sw6 2 min, sw7 3 min, sw8 seconds.
# Run from the repository root: Rscript run_all.R
# =============================================================================

steps <- c(
  "scripts/sw1_candidates.R",           # candidate locations and weights per postal code
  "scripts/sw2_candidate_exposures.R",  # EVI and MSDI at candidates; weighted centroids
  "scripts/sw3_postal_values.R",        # methods A and B per postal code
  "scripts/sw4_compare.R",              # address-level comparison, bootstrap, models
  "scripts/sw5_sensitivity.R",          # modal quintile rule and sensitivity analyses 1-6
  "scripts/sw6_figures_tables.R",       # manuscript tables, figures, supplementary tables
  "scripts/sw7_pcweight_check.R",       # checks of the portable function R/pcweight.R
  "scripts/sw8_manuscript_tables.R")    # manuscript tables as Markdown

for (s in steps) {
  message("\n==== ", s, " (", format(Sys.time(), "%H:%M"), ")")
  rc <- system2(file.path(R.home("bin"), "Rscript"), shQuote(s))
  if (rc != 0) stop("Step failed: ", s)
}
message("\nAnalysis complete. Outputs are in data/derived; compare with outputs/.")
