# =============================================================================
# SW1  Candidate locations and population weights for each residential postal code
# Spatial-weight protocol, Section 2.3 (Table 1, Equations 1-2)
#
# For every postal code of a PCCF+-coded NAR address in the validation study
# (VAL_DERIVED/linked/linked_<pr>.rds), lists the PCCF+ 8D records from which the residential
# routine draws the assigned location, given the Link_Source PCCF+ reported:
#   F  pccfuniq   single record                         weight 1
#   C  wc6dups    DA representative points, by PCODE    PC6DAWT
#   D  pccfdups   blockface/block/DA points, by PCODE   DA population share (Eq. 2): DA21pop of
#                                                       the candidate DA / number of candidate
#                                                       rows of the postal code in that DA
#   5  wc5dups    DA points, by first 5 characters      PC5DAWT
#   4  wc4dups    DA points, by first 4 characters      PC4DAWT
#   3, 2, 1, 0 and uncoded records are excluded (decision 2026-10-14). PCCF+ 8D draws a
#   Link_Source 3 record from the 3-character WCF but then sets DA = 9999 and does not keep
#   the drawn point (engine Step 20; RPF 8), so like 2, 1 and 0 it has only partial geography.
# Raw weights are normalized within the postal code (Eq. 1); NA weights count as 0, and a
# postal code whose raw weights are all 0 or NA gets equal weights (flag w_fallback), as in
# the PCCF+ run of the validation study (VALIDATION/R/pccf_engine_patch.R). w_equal = 1/m is
# the PCCF+ 8D selection probability for Link_Source D and the sensitivity weight.
#
# Outputs (data/derived/sw1/):
#   cand_<pr>.rds      one row per postal code x candidate: pr, PCODE, Link_Source, cand_id,
#                      DAuid, lat, lon (NAD83), x, y (EPSG:3347), Rep_Pt_Type, w_raw, w_pop,
#                      w_equal, w_fallback
#   pcode_<pr>.rds     one row per postal code: Link_Source, n_addr, n_cand, n_eff (w_pop),
#                      n_eff_equal, n_da, w_fallback, realized-point check counts
#   candidates.rds, postal_codes.rds   all provinces (written when all provinces are done)
#   sw1_summary.csv    postal codes, addresses and candidates by Link_Source
# Checks (the script stops if any fails):
#   1. w_pop and w_equal sum to 1 in every postal code
#   2. Link_Source counts (addresses and postal codes) match P1; every postal code has
#      candidates; Link_Source F postal codes have exactly one candidate
#   3. the realized PCCF+ point (same DA, within 1 m) is in the candidate set for
#      >= 99.9% of Link_Source C, D, 5 and 4 addresses
#
# Run: Rscript scripts/sw1_candidates.R           (all provinces)
#      SW_PROV=11 Rscript scripts/sw1_candidates.R (one province, for testing)
#      SW_FORCE=1 ...  recomputes provinces that are already done
# Log: data/derived/sw1_candidates.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(sf) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
PCCF_R <- file.path(PCCF_DIR, "DATA_R")
LINKED <- file.path(VAL_DERIVED, "linked")
OUT    <- file.path(DERIVED, "sw1")
LOG    <- file.path(DERIVED, "sw1_candidates.log")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CRS_LAMBERT <- 3347     # as in scripts/p1_data_assembly.R
CRS_NAD83   <- 4269     # PCCF+ LAT/LONG
TOL_M       <- 1        # realized point = candidate if same DA and within 1 m
MIN_REALIZED <- 0.999
FORCE <- nzchar(Sys.getenv("SW_FORCE"))
PROVS <- if (nzchar(Sys.getenv("SW_PROV"))) strsplit(Sys.getenv("SW_PROV"), ",")[[1]] else
  sub("^linked_(\\d+)\\.rds$", "\\1", list.files(LINKED, "^linked_\\d+\\.rds$"))

msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste0(...))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK   " else "CHECK FAIL ", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}

# ---- PCCF+ candidate files ---------------------------------------------------------------

to_lambert <- function(lat, lon) {
  xy <- st_coordinates(st_transform(st_as_sf(data.table(lon = lon, lat = lat), coords = c("lon", "lat"),
                                             crs = CRS_NAD83), CRS_LAMBERT))
  list(xy[, 1], xy[, 2])
}

load_pccf <- function() {
  rd <- function(f) as.data.table(readRDS(file.path(PCCF_R, paste0(f, ".rds"))))
  std <- function(d, key, w) {
    d[, .(key = get(key), DAuid = paste0(PR, CD, DA), lat = as.numeric(LAT), lon = as.numeric(LONG),
          Rep_Pt_Type = if ("Rep_Pt_Type" %in% names(d)) as.character(Rep_Pt_Type) else NA_character_,
          w_raw = if (is.null(w)) NA_real_ else as.numeric(get(w)))]
  }
  src <- list(
    F = std(rd("pccf_pccfuniq"), "PCODE", NULL),
    C = std(rd("pccf_wc6dups"),  "PCODE", "PC6DAWT"),
    D = std(rd("pccf_pccfdups"), "PCODE", NULL),
    `5` = std(rd("pccf_wc5dups"), "PCODE5", "PC5DAWT"),
    `4` = std(rd("pccf_wc4dups"), "PCODE4", "PC4DAWT"))
  cand <- rbindlist(src, idcol = "Link_Source")
  # Link_Source D: DA population share split equally among the postal code's rows in that DA
  dapop <- rd("geo_dablkpnt")[, .(DAuid, DA21pop = as.numeric(DA21pop))]
  cand[dapop, DA21pop := i.DA21pop, on = "DAuid"]
  cand[Link_Source == "D", w_raw := DA21pop / .N, by = .(key, DAuid)]
  cand[Link_Source == "F", w_raw := 1]
  cand[, DA21pop := NULL]
  ok <- cand[, !is.na(lat) & !is.na(lon)]
  if (!all(ok)) msg("PCCF+ candidate rows without coordinates dropped: ", sum(!ok))
  cand <- cand[ok]
  cand[, c("x", "y") := to_lambert(lat, lon)]
  setkey(cand, Link_Source, key)
  cand
}

# ---- One province -------------------------------------------------------------------------------

KEEP <- c("F", "C", "D", "5", "4")    # Link_Source 3, 2, 1, 0: partial geography, excluded
key_of <- function(pcode, ls) fcase(ls %chin% c("F", "C", "D"), pcode, ls == "5", substr(pcode, 1, 5),
                                    ls == "4", substr(pcode, 1, 4))

do_province <- function(pr, cand) {
  f_c <- file.path(OUT, paste0("cand_", pr, ".rds")); f_p <- file.path(OUT, paste0("pcode_", pr, ".rds"))
  if (!FORCE && file.exists(f_c) && file.exists(f_p)) { msg(pr, ": done earlier, skipped"); return(invisible()) }
  t0 <- Sys.time()
  a <- readRDS(file.path(LINKED, paste0("linked_", pr, ".rds")))[pccf_coded == TRUE,
         .(PCODE, Link_Source = pccf_Link_Source, pccf_DAUID, pccf_x, pccf_y)]
  p1_counts <- a[, .(n_addr = .N), keyby = Link_Source]
  excl <- p1_counts[!Link_Source %chin% KEEP]
  if (nrow(excl)) msg(pr, ": excluded (partial geography): ", paste(sprintf("Link_Source %s %s", excl$Link_Source,
                                                                     format(excl$n_addr, big.mark = ",")), collapse = "; "))
  a <- a[Link_Source %chin% KEEP]

  # Link_Source is a property of the postal code; stop if it is not
  ls_n <- a[, uniqueN(Link_Source), by = PCODE]
  check(all(ls_n$V1 == 1), sprintf("%s: one Link_Source per postal code (%d postal codes with more)", pr, sum(ls_n$V1 > 1)))
  pc <- a[, .(n_addr = .N), keyby = .(PCODE, Link_Source)]
  pc[, key := key_of(PCODE, Link_Source)]

  # candidates for each postal code (5/4/3: the prefix set, repeated for each postal code)
  cp <- cand[pc[, .(Link_Source, key, PCODE)], on = c("Link_Source", "key"), allow.cartesian = TRUE, nomatch = NULL]
  cp[, cand_id := seq_len(.N), by = PCODE]
  cp[, `:=`(m = .N, w0 = fifelse(is.na(w_raw) | w_raw < 0, 0, w_raw)), by = PCODE]
  cp[, w_fallback := sum(w0) <= 0, by = PCODE]
  cp[, w_pop := if (w_fallback[1]) 1 / .N else w0 / sum(w0), by = PCODE]
  cp[, w_equal := 1 / .N, by = PCODE]

  # postal-code summary
  ps <- cp[, .(n_cand = .N, n_eff = 1 / sum(w_pop^2), n_eff_equal = .N, n_da = uniqueN(DAuid),
               w_fallback = w_fallback[1]), keyby = PCODE]
  ps <- ps[pc, on = "PCODE"]
  ps[is.na(n_cand), n_cand := 0L]

  # realized PCCF+ point: same DA and within TOL_M of a candidate
  r <- a[Link_Source != "F" & !is.na(pccf_x)][, rid := .I]
  hit <- cp[r, on = .(PCODE, DAuid = pccf_DAUID), nomatch = NULL, allow.cartesian = TRUE,
            .(rid, d = sqrt((x - pccf_x)^2 + (y - pccf_y)^2))][, .(hit = min(d) <= TOL_M), by = rid]
  r[, realized_in_set := FALSE][hit[hit == TRUE], realized_in_set := TRUE, on = "rid"]
  rs <- r[, .(n_realized_checked = .N, n_realized_in_set = sum(realized_in_set)), by = PCODE]
  ps[rs, `:=`(n_realized_checked = i.n_realized_checked, n_realized_in_set = i.n_realized_in_set), on = "PCODE"]
  ps[is.na(n_realized_checked), `:=`(n_realized_checked = 0L, n_realized_in_set = 0L)]
  ps[, pr := pr]; setcolorder(ps, c("pr", "PCODE", "Link_Source"))

  cp[, pr := pr]
  cp <- cp[, .(pr, PCODE, Link_Source, cand_id, DAuid, lat, lon, x, y, Rep_Pt_Type, w_raw, w_pop, w_equal, w_fallback)]

  # ---- checks for this province
  wsum <- cp[, .(s_pop = sum(w_pop), s_eq = sum(w_equal)), by = PCODE]
  check(wsum[, all(abs(s_pop - 1) < 1e-9 & abs(s_eq - 1) < 1e-9)], sprintf("%s: weights sum to 1 in all %s postal codes", pr, format(nrow(wsum), big.mark = ",")))
  check(ps[, all(n_cand > 0)], sprintf("%s: every postal code has candidates (%d without)", pr, ps[n_cand == 0, .N]))
  check(ps[Link_Source == "F", all(n_cand == 1)], sprintf("%s: Link_Source F postal codes have one candidate", pr))
  out_counts <- ps[, .(n_addr = sum(n_addr)), keyby = Link_Source]
  cmp <- merge(p1_counts[Link_Source %chin% out_counts$Link_Source], out_counts, by = "Link_Source", suffixes = c("_p1", "_sw1"))
  check(nrow(cmp) == nrow(out_counts) && cmp[, all(n_addr_p1 == n_addr_sw1)], sprintf("%s: address counts by Link_Source match P1 (%s)", pr,
        paste(sprintf("%s %s", cmp$Link_Source, format(cmp$n_addr_sw1, big.mark = ",")), collapse = "; ")))
  share <- r[, mean(realized_in_set)]
  msg(pr, sprintf(": realized PCCF+ point in candidate set for %.4f%% of %s C/D/5/4 addresses", 100 * share, format(nrow(r), big.mark = ",")))
  if (nrow(r)) check(share >= MIN_REALIZED, sprintf("%s: realized point in candidate set >= %.1f%%", pr, 100 * MIN_REALIZED))

  saveRDS(cp, f_c); saveRDS(ps, f_p)
  msg(pr, sprintf(": %s postal codes, %s candidate rows, %.1f s", format(nrow(ps), big.mark = ","),
                  format(nrow(cp), big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---- Run ------------------------------------------------------------------------------------------

msg("SW1 start; provinces: ", paste(PROVS, collapse = " "), if (FORCE) " (FORCE)" else "")
cand <- load_pccf()
msg("PCCF+ candidate rows: ", paste(sprintf("%s %s", cand[, .N, keyby = Link_Source]$Link_Source,
                                            format(cand[, .N, keyby = Link_Source]$N, big.mark = ",")), collapse = "; "))
for (pr in PROVS) do_province(pr, cand)

# ---- Combine and summarize when every province is done ------------------------------------------
all_pr <- sub("^linked_(\\d+)\\.rds$", "\\1", list.files(LINKED, "^linked_\\d+\\.rds$"))
done <- all(file.exists(file.path(OUT, paste0("pcode_", all_pr, ".rds"))))
use_pr <- if (done) all_pr else PROVS
ps <- rbindlist(lapply(use_pr, function(p) readRDS(file.path(OUT, paste0("pcode_", p, ".rds")))))
summ <- ps[, .(postal_codes = .N, addresses = sum(n_addr), candidates = sum(n_cand),
               mean_n_cand = round(mean(n_cand), 2), addr_weighted_n_eff = round(sum(n_addr * n_eff) / sum(n_addr), 2),
               pct_multi_da = round(100 * mean(n_da > 1), 1), fallback_pcodes = sum(w_fallback),
               realized_in_set_pct = round(100 * sum(n_realized_in_set) / pmax(sum(n_realized_checked), 1), 3)),
           keyby = Link_Source]
fwrite(summ, file.path(OUT, if (done) "sw1_summary.csv" else paste0("sw1_summary_", paste(PROVS, collapse = "_"), ".csv")))
print(summ)
if (done) {
  check(ps[, sum(n_realized_in_set) / sum(n_realized_checked)] >= MIN_REALIZED,
        "national: realized point in candidate set >= 99.9% (C/D/5/4)")
  saveRDS(ps, file.path(OUT, "postal_codes.rds"))
  saveRDS(rbindlist(lapply(all_pr, function(p) readRDS(file.path(OUT, paste0("cand_", p, ".rds"))))),
          file.path(OUT, "candidates.rds"))
  msg("All provinces combined: ", format(nrow(ps), big.mark = ","), " postal codes")
}
msg("SW1 done")
