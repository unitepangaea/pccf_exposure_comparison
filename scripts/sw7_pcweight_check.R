# =============================================================================
# SW7  Checks of the portable function R/pcweight.R
# Spatial-weight protocol, Section 2.8 (Table 3, last row)
#
# 1. Agreement with the analysis pipeline for every postal code of one province (default
#    Prince Edward Island, 11), for both weight sets (d_weights = "equal" = primary PCCF+
#    selection probabilities, sw3 wset "pccf"; "pop" = DA population share, wset "pop"):
#    number of candidates, method A EVI (9 measures: mean, SD^2 = V_p, share missing), MSDI
#    scores, quintiles, pi and modal quintile (sw5), method B centroid, EVI, DA, score and
#    quintile. The exposure is a raster on the MODIS grid holding the EVI values of the
#    validation cache (VAL_DERIVED/p2/evi_points.rds and sw2/evi_cells_new.rds), so the test also
#    covers the raster path (containing cell, reprojection to the sinusoidal grid).
# 2. Link source inferred by the function (PCCF+ 8D matching sequence) vs the link source
#    reported by PCCF+ for every postal code of a PCCF+-coded NAR address (national).
#
# Outputs (data/derived/sw7/): pcweight_<pr>.rds (function output, equal weights),
#   check_<pr>.csv (agreement by column and weight set), link_source_crosstab.csv,
#   link_source_mismatches.csv
# Checks (stop if any fails):
#   1. every compared column agrees for every postal code of the province (|diff| < 1e-9,
#      identical DA and quintiles), both weight sets
#   2. inferred link source = PCCF+ link source for >= 99.9% of postal codes, and the
#      inferred link source is F, C, D, 5 or 4 for >= 99.9% of postal codes analysed
#
# Run: Rscript scripts/sw7_pcweight_check.R          (PEI + national link source)
#      SW_PROV=12 Rscript scripts/sw7_pcweight_check.R
# Log: data/derived/sw7_pcweight_check.log
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(sf); library(terra); library(ncdf4); library(readxl) })

ROOT <- (function(d) { while (!file.exists(file.path(d, "config.R"))) { if (dirname(d) == d) stop("Run from inside the repository"); d <- dirname(d) }; d })(getwd())
source(file.path(ROOT, "config.R"))   # sets VALIDATION, RAW, PCCF_DIR, VAL_DERIVED and DERIVED
SW   <- DERIVED
OUT  <- file.path(SW, "sw7"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG  <- file.path(SW, "sw7_pcweight_check.log")
PR   <- if (nzchar(Sys.getenv("SW_PROV"))) Sys.getenv("SW_PROV") else "11"
PCCF_R    <- file.path(PCCF_DIR, "DATA_R")
DA_SHP    <- file.path(RAW, "lda_000b21a_e", "lda_000b21a_e.shp")
MSDI_XLSX <- file.path(RAW, "Canada2021_en", "A-MSDIData_Can2021_en", "1. EquivalenceTableCanada2021_en.xlsx")
TOL <- 1e-9

msg <- function(...) {
  line <- gsub("  +", " ", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(...)))
  cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE)
}
check <- function(ok, what) {
  msg(if (ok) "CHECK OK  " else "CHECK FAIL", what)
  if (!ok) stop("Check failed: ", what, call. = FALSE)
}
fmt <- function(n) format(n, big.mark = ",")

source(file.path(ROOT, "R", "pcweight.R"))

# MODIS grid of the validation study (only the grid objects of VALIDATION/scripts/p2_variables.R are evaluated)
P2_OBJECTS <- c("NC_FILE", "CRS_SINU", "nc", "gx", "gy", "res", "nx", "ny", "YEARS")
for (e in parse(file.path(VALIDATION, "scripts", "p2_variables.R"), keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% P2_OBJECTS)
    tryCatch(eval(e, globalenv()), error = function(err) NULL)
if (!exists("gx")) {   # NC_FILE is defined from OUT/ROOT objects of P2; open the file directly
  nc <- nc_open(file.path(RAW, "MODIS", "MOD13Q1.061_250m_aid0001.nc"))
  gx <- nc$dim$xdim$vals; gy <- nc$dim$ydim$vals; res <- abs(gx[2] - gx[1]); nx <- length(gx); ny <- length(gy)
}
stopifnot(exists("CRS_SINU"), exists("gx"))
EVI_M <- c(as.vector(outer(c("pix", "b1", "b2"), c(2025, 2026), paste, sep = "_")), paste0(c("pix", "b1", "b2"), "_mean"))

msg("SW7 start | province", PR)
t0 <- Sys.time(); pccf <- pcw_read_pccf(PCCF_R)
msg("PCCF+ files read in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s |", fmt(nrow(pccf$cand)), "candidate rows")

# ---- 1. Inputs for one province ---------------------------------------------------------------------
ps  <- readRDS(file.path(SW, "sw1", "postal_codes.rds"))[pr == PR, .(PCODE, Link_Source, n_addr)]
pv  <- readRDS(file.path(SW, "sw3", "postal_values.rds"))[pr == PR]
mq  <- readRDS(file.path(SW, "sw5", "modal_quintiles.rds"))[pr == PR]
ce  <- readRDS(file.path(SW, "sw2", "cand_exposures.rds"))[pr == PR, .(cell)]
cen <- readRDS(file.path(SW, "sw2", "centroids.rds"))[pr == PR, .(cell)]
msg("Postal codes:", fmt(nrow(ps)))

# raster of EVI values on the MODIS grid, covering the province's candidate and centroid cells
cells <- unique(na.omit(c(ce$cell, cen$cell)))
rows <- as.integer((cells - 1) %/% nx) + 1L; cols <- as.integer((cells - 1) %% nx) + 1L
r1 <- min(rows) - 2L; r2 <- max(rows) + 2L; c1 <- min(cols) - 2L; c2 <- max(cols) + 2L
evi <- rbind(readRDS(file.path(VAL_DERIVED, "p2", "evi_points.rds")), readRDS(file.path(SW, "sw2", "evi_cells_new.rds")))
evi <- unique(evi, by = "pid")
evi[, `:=`(row = as.integer((pid - 1) %/% nx) + 1L, col = as.integer((pid - 1) %% nx) + 1L)]
evi <- evi[row %between% c(r1, r2) & col %between% c(c1, c2)]
for (w in c("pix", "b1", "b2")) {
  v <- rowMeans(evi[, paste0(w, "_", YEARS), with = FALSE], na.rm = TRUE)
  evi[, (paste0(w, "_mean")) := fifelse(is.nan(v), NA_real_, v)]
}
nr <- r2 - r1 + 1L; ncl <- c2 - c1 + 1L
r <- rast(nrows = nr, ncols = ncl, xmin = gx[c1] - res / 2, xmax = gx[c2] + res / 2,
          ymin = gy[r2] - res / 2, ymax = gy[r1] + res / 2, crs = CRS_SINU, nlyrs = length(EVI_M))
names(r) <- paste0("evi_", EVI_M)
idx <- (evi$row - r1) * ncl + (evi$col - c1) + 1L
for (k in seq_along(EVI_M)) { v <- rep(NA_real_, nr * ncl); v[idx] <- evi[[EVI_M[k]]]; values(r[[k]]) <- v }
msg("Test raster:", nr, "x", ncl, "cells,", length(EVI_M), "layers;", fmt(nrow(evi)), "cells with cached EVI")
rm(evi, ce, cen); invisible(gc())

msdi <- as.data.table(read_excel(MSDI_XLSX, sheet = "Data"))[
  , .(DAuid = as.character(DA), SCOREMAT, SCORESOC, QuintMat = as.integer(QuintMat), QuintSoc = as.integer(QuintSoc))]
da <- st_read(DA_SHP, quiet = TRUE)[, "DAUID"]

# ---- 2. Run the function and compare with steps 3 and 5 -----------------------------------------------
same <- function(u, v) (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & abs(u - v) < TOL)
same_chr <- function(u, v) (is.na(u) & is.na(v)) | (!is.na(u) & !is.na(v) & u == v)
cmp_all <- list()
for (ws in c(pccf = "equal", pop = "pop")) {
  wname <- names(which(c(pccf = "equal", pop = "pop") == ws))
  t0 <- Sys.time()
  res <- pcweight(ps, pccf, pcode = "PCODE", link_source = "Link_Source", exposure = r,
                  da_attr = msdi, da_scores = c(mat = "SCOREMAT", soc = "SCORESOC"),
                  da_quintiles = c(mat = "QuintMat", soc = "QuintSoc"), da_boundaries = da, d_weights = ws)
  msg(sprintf("pcweight (%s weights): %s postal codes in %.1f s", ws, fmt(nrow(res)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (ws == "equal") saveRDS(res, file.path(OUT, paste0("pcweight_", PR, ".rds")))
  ref <- pv[wset == wname][mq[wset == wname, .(PCODE, A_mode_mat, A_mode_soc)], on = "PCODE"]
  m <- merge(res, ref, by = "PCODE", suffixes = c("", ".ref"))
  check(nrow(m) == nrow(ref) && nrow(m) == nrow(res), sprintf("%s weights: same postal codes as step 3 (%s)", wname, fmt(nrow(m))))
  num <- c(n_cand = "n_cand",
           setNames(paste0("A_evi_", EVI_M), paste0("A_evi_", EVI_M)),
           setNames(paste0("A_miss_evi_", EVI_M), paste0("A_miss_evi_", EVI_M)),
           setNames(paste0("B_evi_", EVI_M), paste0("B_evi_", EVI_M)),
           A_S_mat = "A_S_mat", A_S_soc = "A_S_soc", A_miss_S_mat = "A_miss_S_mat", A_miss_S_soc = "A_miss_S_soc",
           setNames(paste0("A_pi_mat_", 1:5), paste0("A_pi_mat_", 1:5)), setNames(paste0("A_pi_soc_", 1:5), paste0("A_pi_soc_", 1:5)),
           A_q_mat = "A_q_mat", A_q_soc = "A_q_soc", A_mode_mat = "A_mode_mat", A_mode_soc = "A_mode_soc",
           B_S_mat = "B_S_mat", B_S_soc = "B_S_soc", B_q_mat = "B_q_mat", B_q_soc = "B_q_soc",
           B_x = "B_x", B_y = "B_y", B_lat = "B_lat", B_lon = "B_lon")
  cmp <- rbindlist(lapply(names(num), function(k) {
    u <- m[[k]]; v <- m[[paste0(num[[k]], ".ref")]]
    data.table(column = k, n = length(u), agree = sum(same(as.numeric(u), as.numeric(v))))
  }))
  # V_p: the function returns the weighted SD
  cmp <- rbind(cmp, rbindlist(lapply(EVI_M, function(e)
    data.table(column = paste0("A_sd_evi_", e, "^2 = A_V"), n = nrow(m),
               agree = sum(same(m[[paste0("A_sd_evi_", e)]]^2, m[[paste0("A_V_evi_", e)]]))))),
    data.table(column = "B_DAuid", n = nrow(m), agree = sum(same_chr(m$B_DAuid, m$B_DAuid.ref))),
    data.table(column = "B_da_method", n = nrow(m), agree = sum(same_chr(m$B_da_method, m$B_da_method.ref))))
  cmp[, wset := wname]
  cmp_all[[wname]] <- cmp
  bad <- cmp[agree < n]
  if (nrow(bad)) print(bad)
  check(nrow(bad) == 0, sprintf("%s weights: function = step 3/5 for all %s postal codes in %d columns", wname, fmt(nrow(m)), nrow(cmp)))
}
fwrite(rbindlist(cmp_all), file.path(OUT, paste0("check_", PR, ".csv")))

# inferred link source for the province
inf <- pcw_link_source(ps$PCODE, pccf)
msg(sprintf("Province %s: inferred link source = PCCF+ link source for %.3f%% of %s postal codes", PR, 100 * mean(inf == ps$Link_Source), fmt(nrow(ps))))

# ---- 3. Inferred vs PCCF+ link source, national ---------------------------------------------------------
lf <- list.files(file.path(VAL_DERIVED, "linked"), "^linked_[0-9]+[.]rds$", full.names = TRUE)
a <- rbindlist(lapply(lf, function(f) readRDS(f)[pccf_coded == TRUE & !is.na(PCODE), .(n_addr = .N), by = .(PCODE, ls_pccf = pccf_Link_Source)]))
n_multi <- a[, uniqueN(ls_pccf), by = PCODE][V1 > 1, .N]
a <- a[, .(n_addr = sum(n_addr), ls_pccf = ls_pccf[which.max(n_addr)]), by = PCODE]
t0 <- Sys.time(); a[, ls_inf := pcw_link_source(PCODE, pccf)]
msg(sprintf("National: %s postal codes (%s with more than one PCCF+ link source); inference %.1f s", fmt(nrow(a)), fmt(n_multi),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
xt <- dcast(a[, .(postal_codes = .N, addresses = sum(n_addr)), by = .(ls_pccf, ls_inf)], ls_pccf ~ ls_inf, value.var = "postal_codes", fill = 0)
fwrite(xt, file.path(OUT, "link_source_crosstab.csv")); print(xt)
fwrite(a[ls_pccf != ls_inf][order(ls_pccf, ls_inf, PCODE)][, .(PCODE, ls_pccf, ls_inf, n_addr)], file.path(OUT, "link_source_mismatches.csv"))
agree <- a[, mean(ls_pccf == ls_inf)]; agree_addr <- a[, sum(n_addr * (ls_pccf == ls_inf)) / sum(n_addr)]
an <- a[ls_pccf %chin% PCW_KEEP]
msg(sprintf("Inferred = PCCF+ link source: %.3f%% of postal codes, %.3f%% of addresses; analysed link sources (F, C, D, 5, 4): %.3f%%",
            100 * agree, 100 * agree_addr, 100 * an[, mean(ls_pccf == ls_inf)]))
check(agree >= 0.999 && an[, mean(ls_pccf == ls_inf)] >= 0.999, "inferred link source = PCCF+ link source for >= 99.9% of postal codes (all and analysed)")
nc_close(nc)
msg("SW7 done")
