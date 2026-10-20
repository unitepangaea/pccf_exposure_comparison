# =============================================================================
# SW2  Exposures at candidate locations and at weighted centroids (method B)
# Spatial-weight protocol, Section 2.4
#
# Input: data/derived/sw1/candidates.rds (postal code x candidate, weights)
#
# 1. Weighted centroids per postal code (method B), three variants:
#      pop     w_pop, projected coordinates (EPSG:3347)          primary
#      equal   w_equal, projected coordinates                    sensitivity 1 (Link_Source D)
#      latlon  w_pop, weighted mean of NAD83 latitude/longitude  sensitivity 2
#    Lat/long are kept for each centroid (positional error in step 4 uses the haversine
#    distance of the validation study).
# 2. EVI (pixel, 3x3 = b1, 5x5 = b2; 2025, 2026, mean) at the MODIS cell of every candidate
#    and centroid, with the extraction functions of VALIDATION/scripts/p2_variables.R (parsed, not
#    rewritten). Cells already in the validation cache VAL_DERIVED/p2/evi_points.rds are reused;
#    new cells are extracted once and kept in sw2/evi_cells_new.rds (the validation cache is
#    not modified).
# 3. MSDI 2021 (national factor scores and quintiles) joined to candidates by DAuid, and the
#    national quintile cut points recovered as the midpoint between the largest score of
#    quintile k and the smallest of quintile k + 1 (Section 2.4).
# 4. DA of each centroid: a postal code with a single candidate keeps that candidate's DA
#    (method B = method C by definition); otherwise point-in-polygon in the 2021 DA boundary
#    file, and the nearest DA (flagged) for centroids outside every DA.
#
# Outputs (data/derived/sw2/, or sw2/test_<pr>/ when SW_PROV is set):
#   cand_exposures.rds   candidates + cell, evi_<w>_<yr>, evi_<w>_mean, SCOREMAT, SCORESOC,
#                        QuintMat, QuintSoc
#   centroids.rds        pr, PCODE, Link_Source, variant, x, y, lat, lon, n_cand, cell, EVI
#                        columns, DAuid, da_method (candidate | pip | nearest), nearest_m, MSDI
#   msdi_cutpoints.csv   cut points by dimension; sw2_summary.csv
# Checks (stop if any fails):
#   1. MSDI cut points separate the quintiles and reproduce the published national quintile
#      of every DA
#   2. EVI at the drawn candidate equals evi_pccf_* in VAL_DERIVED/p2/p2_analysis.rds (all
#      windows and years) for >= 99.9% of addresses of each Link_Source (F, C, D, 5, 4)
#   3. every centroid has a DA
#
# Run: SW_PROV=11 Rscript scripts/sw2_candidate_exposures.R   (test)
#      Rscript scripts/sw2_candidate_exposures.R             (national)
#      SW_FORCE=1 recomputes outputs (the new-cell EVI cache is always reused)
# Log: data/derived/sw2_candidate_exposures.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(sf); library(ncdf4); library(readxl) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW     <- DERIVED
IN     <- file.path(SW, "sw1")
PROVS  <- if (nzchar(Sys.getenv("SW_PROV"))) strsplit(Sys.getenv("SW_PROV"), ",")[[1]] else NULL
OUT    <- if (is.null(PROVS)) file.path(SW, "sw2") else file.path(SW, "sw2", paste0("test_", paste(PROVS, collapse = "_")))
F_NEWCELLS <- file.path(SW, "sw2", "evi_cells_new.rds")
LOG    <- file.path(SW, "sw2_candidate_exposures.log")
FORCE  <- nzchar(Sys.getenv("SW_FORCE"))
NC_FILE   <- file.path(RAW, "MODIS", "MOD13Q1.061_250m_aid0001.nc")
DA_SHP    <- file.path(RAW, "lda_000b21a_e", "lda_000b21a_e.shp")
MSDI_XLSX <- file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en", "1. EquivalenceTableCanada2021_en.xlsx")
TOL_M <- 1; MIN_EVI_MATCH <- 0.999
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  line <- gsub("  +", " ", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(...)))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK   " else "CHECK FAIL ", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}
fmt <- function(n) format(n, big.mark = ",")

# ---- Reuse the EVI functions and settings of the validation study (P2) ---------------------
# Only these top-level assignments of VALIDATION/scripts/p2_variables.R are evaluated, so the grid,
# quality filter, growing season and windows are identical to the validation study.
P2_OBJECTS <- c("EVI_VAR", "QA_VAR", "GS_START", "GS_END", "YEARS", "BUFFERS", "STRIP",
                "CRS_LAMBERT", "CRS_SINU", "R_EARTH", "haversine", "qa_ok", "to_cell", "extract_evi",
                "nc", "gx", "gy", "res", "nx", "ny", "dates")
for (e in parse(file.path(VALIDATION, "scripts", "p2_variables.R"), keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% P2_OBJECTS) eval(e, globalenv())
stopifnot(all(P2_OBJECTS %in% ls(globalenv())))
CRS_NAD83 <- 4269
WINDOWS <- c("pix", names(BUFFERS))
EVI_COLS <- c(as.vector(outer(WINDOWS, YEARS, paste, sep = "_")), paste0(WINDOWS, "_mean"))

msg("SW2 start", if (is.null(PROVS)) " (national)" else paste0(" (test: ", paste(PROVS, collapse = ","), ")"),
    if (FORCE) " FORCE" else "")
f_cand <- file.path(OUT, "cand_exposures.rds"); f_cen <- file.path(OUT, "centroids.rds")
done <- !FORCE && file.exists(f_cand) && file.exists(f_cen)
if (done) msg("Outputs exist; running checks only (SW_FORCE=1 to recompute)")

# ---- 1. MSDI and cut points ------------------------------------------------------------------
msdi <- as.data.table(read_excel(MSDI_XLSX, sheet = "Data"))[
  , .(DAuid = as.character(DA), SCOREMAT, SCORESOC, QuintMat = as.integer(QuintMat), QuintSoc = as.integer(QuintSoc))]
cuts <- rbindlist(lapply(c(MAT = "Mat", SOC = "Soc"), function(d) {
  s <- if (d == "Mat") "SCOREMAT" else "SCORESOC"; q <- paste0("Quint", d)
  r <- msdi[!is.na(get(q)), .(lo = min(get(s)), hi = max(get(s))), keyby = .(k = get(q))]
  data.table(k = 1:4, max_k = r$hi[1:4], min_k1 = r$lo[2:5], cut = (r$hi[1:4] + r$lo[2:5]) / 2)
}), idcol = "dimension")
fwrite(cuts, file.path(OUT, "msdi_cutpoints.csv"))
quint_of <- function(score, cut) ifelse(is.na(score), NA_integer_, 1L + rowSums(outer(score, cut, ">")))
check(cuts[, all(max_k < min_k1)], "MSDI quintiles are separated by score (max of k < min of k + 1, both dimensions)")
rep_mat <- msdi[!is.na(QuintMat), all(quint_of(SCOREMAT, cuts[dimension == "MAT", cut]) == QuintMat)]
rep_soc <- msdi[!is.na(QuintSoc), all(quint_of(SCORESOC, cuts[dimension == "SOC", cut]) == QuintSoc)]
check(rep_mat && rep_soc, sprintf("cut points reproduce the national quintile of all %s DAs with MSDI (material and social)",
                                  fmt(msdi[!is.na(QuintMat), .N])))

if (!done) {
  # ---- 2. Candidates and centroids ------------------------------------------------------------
  cand <- readRDS(file.path(IN, "candidates.rds"))
  if (!is.null(PROVS)) cand <- cand[pr %chin% PROVS]
  cand[, n_cand := .N, by = .(pr, PCODE)]
  msg("Candidates: ", fmt(nrow(cand)), " rows, ", fmt(uniqueN(cand$PCODE)), " postal codes")

  cen <- rbind(
    cand[, .(variant = "pop",   x = sum(w_pop * x),   y = sum(w_pop * y),   lat = NA_real_, lon = NA_real_), by = .(pr, PCODE, Link_Source, n_cand)],
    cand[, .(variant = "equal", x = sum(w_equal * x), y = sum(w_equal * y), lat = NA_real_, lon = NA_real_), by = .(pr, PCODE, Link_Source, n_cand)],
    cand[, .(variant = "latlon", x = NA_real_, y = NA_real_, lat = sum(w_pop * lat), lon = sum(w_pop * lon)), by = .(pr, PCODE, Link_Source, n_cand)])
  # projected <-> geographic for the centroids
  i <- cen[variant != "latlon", which = TRUE]
  ll <- st_coordinates(st_transform(st_as_sf(cen[i, .(x, y)], coords = c("x", "y"), crs = CRS_LAMBERT), CRS_NAD83))
  cen[i, `:=`(lon = ll[, 1], lat = ll[, 2])]
  i <- cen[variant == "latlon", which = TRUE]
  xy <- st_coordinates(st_transform(st_as_sf(cen[i, .(lon, lat)], coords = c("lon", "lat"), crs = CRS_NAD83), CRS_LAMBERT))
  cen[i, `:=`(x = xy[, 1], y = xy[, 2])]
  rm(ll, xy, i)
  msg("Centroids: ", fmt(nrow(cen)), " (", fmt(cen[variant == "pop", .N]), " postal codes x 3 variants)")

  # ---- 3. EVI at candidate and centroid cells ---------------------------------------------------
  cell_of <- function(px, py) { rc <- to_cell(px, py); fifelse(is.na(rc$col), NA_real_, (rc$row - 1) * nx + rc$col) }
  cand[, cell := cell_of(x, y)]; cen[, cell := cell_of(x, y)]
  cells <- unique(c(cand$cell, cen$cell)); cells <- cells[!is.na(cells)]
  cache <- readRDS(file.path(VAL_DERIVED, "p2", "evi_points.rds"))
  newc <- if (file.exists(F_NEWCELLS)) readRDS(F_NEWCELLS) else cache[0]
  need <- setdiff(cells, c(cache$pid, newc$pid))
  msg("MODIS cells: ", fmt(length(cells)), " needed; ", fmt(sum(cells %in% cache$pid)), " in the validation cache; ",
      fmt(sum(cells %in% newc$pid)), " in the SW2 cache; ", fmt(length(need)), " to extract")
  if (length(need)) {
    pts <- data.table(pid = need, row = as.integer((need - 1) %/% nx) + 1L, col = as.integer((need - 1) %% nx) + 1L)
    t0 <- Sys.time()
    newc <- rbind(newc, extract_evi(pts))
    saveRDS(newc, F_NEWCELLS)
    msg("Extracted EVI at ", fmt(length(need)), " new cells in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
  evi <- rbind(cache[pid %in% cells], newc[pid %in% cells])
  setnames(evi, setdiff(names(evi), "pid"), paste0("evi_", setdiff(names(evi), "pid")))
  for (w in WINDOWS) {               # 2025-26 mean, as in P2
    v <- rowMeans(evi[, paste0("evi_", w, "_", YEARS), with = FALSE], na.rm = TRUE)
    evi[, (paste0("evi_", w, "_mean")) := fifelse(is.nan(v), NA_real_, v)]
  }
  cand <- merge(cand, evi, by.x = "cell", by.y = "pid", all.x = TRUE, sort = FALSE)
  cen  <- merge(cen,  evi, by.x = "cell", by.y = "pid", all.x = TRUE, sort = FALSE)
  rm(cache, newc, evi); invisible(gc())

  # ---- 4. MSDI at candidates; DA and MSDI at centroids ---------------------------------------------
  cand <- merge(cand, msdi, by = "DAuid", all.x = TRUE, sort = FALSE)
  cen[, `:=`(DAuid = NA_character_, da_method = NA_character_, nearest_m = NA_real_)]
  one <- cand[n_cand == 1, .(pr, PCODE, cand_DAuid = DAuid)]
  cen[one, `:=`(DAuid = i.cand_DAuid, da_method = "candidate"), on = .(pr, PCODE)]
  da <- st_read(DA_SHP, quiet = TRUE)[, "DAUID"]
  if (!isTRUE(st_crs(da)$epsg == CRS_LAMBERT)) da <- st_transform(da, CRS_LAMBERT)
  pip <- function(px, py) {
    pts <- st_as_sf(data.table(x = px, y = py), coords = c("x", "y"), crs = CRS_LAMBERT)
    hit <- vapply(st_intersects(pts, da), function(v) if (length(v)) v[1] else NA_integer_, 1L)
    near <- is.na(hit); nd <- rep(NA_real_, length(px))
    if (any(near)) {
      hit[near] <- st_nearest_feature(pts[near, ], da)
      nd[near] <- as.numeric(st_distance(pts[near, ], da[hit[near], ], by_element = TRUE))
    }
    list(DAuid = as.character(da$DAUID[hit]), da_method = fifelse(near, "nearest", "pip"), nearest_m = nd)
  }
  i <- cen[is.na(DAuid), which = TRUE]
  t0 <- Sys.time(); r <- pip(cen$x[i], cen$y[i])
  cen[i, `:=`(DAuid = r$DAuid, da_method = r$da_method, nearest_m = r$nearest_m)]
  msg("Centroid DA by point-in-polygon for ", fmt(length(i)), " centroids in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s; ",
      fmt(cen[da_method == "nearest", .N]), " outside every DA (nearest DA)")
  # diagnostic: for single-candidate postal codes, does point-in-polygon agree with the PCCF+ DA?
  j <- cen[variant == "pop" & da_method == "candidate", which = TRUE]
  if (length(j)) {
    r1 <- pip(cen$x[j], cen$y[j])
    msg(sprintf("Diagnostic: single-candidate postal codes whose point-in-polygon DA equals the PCCF+ DA: %.2f%% of %s",
                100 * mean(r1$DAuid == cen$DAuid[j]), fmt(length(j))))
  }
  cen <- merge(cen, msdi, by = "DAuid", all.x = TRUE, sort = FALSE)
  rm(da); invisible(gc())

  setcolorder(cand, c("pr", "PCODE", "Link_Source", "cand_id", "n_cand", "DAuid", "lat", "lon", "x", "y", "Rep_Pt_Type",
                      "w_raw", "w_pop", "w_equal", "w_fallback", "cell"))
  setorder(cand, pr, PCODE, cand_id)
  setcolorder(cen, c("pr", "PCODE", "Link_Source", "variant", "n_cand", "x", "y", "lat", "lon", "cell"))
  setorder(cen, pr, PCODE, variant)
  saveRDS(cand, f_cand); saveRDS(cen, f_cen)
  msg("Saved ", f_cand, " and ", f_cen)
}

# ---- 5. Checks ---------------------------------------------------------------------------------------
cand <- readRDS(f_cand); cen <- readRDS(f_cen)
check(cen[, !anyNA(DAuid)], sprintf("every centroid has a DA (%s centroids)", fmt(nrow(cen))))

# EVI at the drawn candidate equals the PCCF+ EVI of the validation study
lf <- list.files(file.path(VAL_DERIVED, "linked"), "^linked_[0-9]+[.]rds$", full.names = TRUE)
if (!is.null(PROVS)) lf <- lf[sub("^linked_([0-9]+)[.]rds$", "\\1", basename(lf)) %in% PROVS]
a <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE & pccf_Link_Source %chin% c("F", "C", "D", "5", "4"),
                                                  .(ADDR_GUID, PCODE, Link_Source = pccf_Link_Source, pccf_DAUID, pccf_x, pccf_y)]))
p2 <- readRDS(file.path(VAL_DERIVED, "p2", "p2_analysis.rds"))
p2 <- p2[, c("ADDR_GUID", paste0("evi_pccf_", EVI_COLS)), with = FALSE]
a <- p2[a, on = "ADDR_GUID", nomatch = NULL]; rm(p2); invisible(gc())
a[, aid := .I]
m <- cand[a, on = .(PCODE, DAuid = pccf_DAUID), nomatch = NULL, allow.cartesian = TRUE,
          c(list(aid = i.aid, d = sqrt((x - i.pccf_x)^2 + (y - i.pccf_y)^2)), .SD), .SDcols = paste0("evi_", EVI_COLS)]
m <- m[d <= TOL_M][order(aid, d)][!duplicated(aid)]
a <- a[m, on = "aid", nomatch = NULL]
same <- Reduce(`&`, lapply(EVI_COLS, function(w) {
  u <- a[[paste0("evi_", w)]]; v <- a[[paste0("evi_pccf_", w)]]
  (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & abs(u - v) < 1e-12)
}))
a[, evi_match := same]
em <- a[, .(addresses = .N, evi_match_pct = round(100 * mean(evi_match), 4)), keyby = Link_Source]
for (k in seq_len(nrow(em))) msg(sprintf("EVI at drawn candidate = P2 PCCF+ EVI: Link_Source %s, %.4f%% of %s addresses",
                                         em$Link_Source[k], em$evi_match_pct[k], fmt(em$addresses[k])))
check(all(em$evi_match_pct >= 100 * MIN_EVI_MATCH), "EVI at the drawn candidate equals evi_pccf_* in P2 for >= 99.9% of addresses (each Link_Source)")

# ---- 6. Summary ---------------------------------------------------------------------------------------
summ <- rbind(
  cand[, .(rows = .N, postal_codes = uniqueN(PCODE), pct_evi_b1_missing = round(100 * mean(is.na(evi_b1_mean)), 3),
           pct_outside_grid = round(100 * mean(is.na(cell)), 3), pct_msdi_missing = round(100 * mean(is.na(QuintMat)), 3)),
       keyby = Link_Source][, table := "candidates"],
  cen[, .(rows = .N, postal_codes = uniqueN(PCODE), pct_evi_b1_missing = round(100 * mean(is.na(evi_b1_mean)), 3),
          pct_outside_grid = round(100 * mean(is.na(cell)), 3), pct_msdi_missing = round(100 * mean(is.na(QuintMat)), 3),
          pct_nearest_da = round(100 * mean(da_method == "nearest"), 3)),
      keyby = .(variant, Link_Source)][, table := "centroids"], fill = TRUE)
setcolorder(summ, "table")
fwrite(summ, file.path(OUT, "sw2_summary.csv")); print(summ)
fwrite(em, file.path(OUT, "sw2_evi_check.csv"))
nc_close(nc)
msg("SW2 done")
