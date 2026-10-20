# =============================================================================
# pcweight.R  Population-weighted exposure assignment from residential postal codes
# Spatial-weight protocol, Section 2.8 (methods A and B of Sections 2.3-2.4)
#
# For each six-character postal code, lists the PCCF+ 8D records from which the residential
# routine (codeVersion 0) draws a location, weights them by the PCCF+ selection probabilities
# and returns, without random allocation:
#   method A  the weighted mean of the exposure over the candidates (with the weighted SD and
#             the share of weight with missing values), the weighted mean of DA attributes
#             such as MSDI factor scores, its quintile, the weighted quintile distribution
#             pi(1..5) and the modal quintile;
#   method B  the weighted centroid of the candidates, with the exposure and the DA (and DA
#             attributes) at that point.
#
# Candidates and weights (Table 1 of the protocol):
#   F  pccfuniq   one record                                   weight 1
#   C  wc6dups    DA representative points, by PCODE           PC6DAWT
#   D  pccfdups   blockface/block/DA points, by PCODE          equal (PCCF+ 8D selection);
#                                                              d_weights = "pop": 2021 DA
#                                                              population / rows in that DA
#   5  wc5dups    DA points, by first 5 characters             PC5DAWT
#   4  wc4dups    DA points, by first 4 characters             PC4DAWT
#   3, 2, 1, 0    no candidate location retained by PCCF+ 8D   returned without values
# Missing or negative raw weights count as 0; a postal code whose raw weights are all 0
# gets equal weights (w_fallback), as in the patched PCCF+ run of the validation study.
#
# Link source: taken from a PCCF+ run when the input has a link-source column; otherwise
# inferred with the residential matching sequence of PCCF+ 8D (engine steps 10-24):
#   wc6point -> C; pccfuniq (DA not 0000) -> F; pointdup -> D; wc5point -> 5; wc4point -> 4;
#   F, D and unmatched postal codes with delivery mode H, J, K, M, R, T or X (or Z with one
#   of these as DMTDIFF) are recoded from the FSA: wc3point -> 3; wc2point -> 2; valid first
#   letter -> 1; otherwise 0. A pccfuniq record with DA 0000 is not coded (0).
# pccf_nonres flags F and D postal codes that PCCF+ codes as non-residential (cpc_egmres
# EGMRES "-"); PCCF+ residential coding returns no coordinates for them, whereas this
# function still returns values.
#
# Exposures: a terra SpatRaster (value of the containing cell, in the raster's CRS; one
# output per layer) or an sf POINT layer (value of the nearest point, in EPSG:3347; one
# output per column named in exposure_cols). Weights are renormalized over candidates with a
# non-missing value. For moving-window means (e.g. 3 x 3 cells), pass a raster of window means
# (terra::focal).
# DA attributes: a table with a DA identifier (8 characters, PR + CD + DA) and score columns;
# quintile columns, when given, define the national cut points (midpoint between the largest
# score of quintile k and the smallest of quintile k + 1, checked to reproduce every DA).
# Method B DA: a postal code with one candidate keeps that candidate's DA; otherwise the DA
# containing the centroid (da_boundaries, an sf polygon layer with DAUID), or the nearest DA
# (flagged) for centroids outside every DA. Without da_boundaries, method B DA attributes are
# returned only for single-candidate postal codes.
#
# PCCF+ files are licensed by Statistics Canada and are not distributed with this code; point
# pcw_read_pccf() to the DATA_R folder of a licensed PCCF+ 8D installation.
#
# Main functions
#   pccf <- pcw_read_pccf("<PCCF+ 8D>/DATA_R")
#   res  <- pcweight(data, pccf, pcode = "PCODE", exposure = r, da_attr = msdi,
#                    da_scores = c(mat = "SCOREMAT", soc = "SCORESOC"),
#                    da_quintiles = c(mat = "QuintMat", soc = "QuintSoc"), da_boundaries = da)
#   pcw_link_source(pcodes, pccf)        inferred link source
#   pcw_candidates(pcodes, ls, pccf)     candidate table with weights
# Requires data.table and sf; terra for raster exposures.
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(sf) })

PCW_CRS <- 3347                 # Statistics Canada Lambert (NAD83), as in the validation study
PCW_RECODE_DMT <- c("H", "J", "K", "M", "R", "T", "X")
PCW_KEEP <- c("F", "C", "D", "5", "4")

pcw_clean_pcode <- function(x) toupper(gsub("\\s+", "", as.character(x)))

# ---- PCCF+ files -----------------------------------------------------------------------------
pcw_read_pccf <- function(dir) {
  rd <- function(f) {
    p <- file.path(dir, paste0(f, ".rds"))
    if (!file.exists(p)) stop("PCCF+ file not found: ", p, call. = FALSE)
    data.table::as.data.table(readRDS(p))
  }
  std <- function(d, key, w) {
    d[, list(key = get(key), DAuid = paste0(PR, CD, DA), lat = as.numeric(LAT), lon = as.numeric(LONG),
             Rep_Pt_Type = if ("Rep_Pt_Type" %in% names(d)) as.character(Rep_Pt_Type) else NA_character_,
             w_raw = if (is.null(w)) NA_real_ else as.numeric(get(w)))]
  }
  uniq <- rd("pccf_pccfuniq"); dups <- rd("pccf_pccfdups")
  rec <- function(d) d$DMT %in% PCW_RECODE_DMT | (d$DMT %in% "Z" & d$DMTDIFF %in% PCW_RECODE_DMT)
  cand <- data.table::rbindlist(list(
    F = std(uniq[DA != "0000"], "PCODE", NULL),
    C = std(rd("pccf_wc6dups"), "PCODE", "PC6DAWT"),
    D = std(dups, "PCODE", NULL),
    `5` = std(rd("pccf_wc5dups"), "PCODE5", "PC5DAWT"),
    `4` = std(rd("pccf_wc4dups"), "PCODE4", "PC4DAWT")), idcol = "Link_Source")
  cand <- cand[!is.na(lat) & !is.na(lon)]
  data.table::setkey(cand, Link_Source, key)
  dapop <- unique(rd("geo_dablkpnt")[, list(DAuid, DA21pop = as.numeric(DA21pop))], by = "DAuid")
  egm <- rd("cpc_egmres")[order(PCODE)][!duplicated(PCODE)]
  structure(list(
    cand  = cand,
    dapop = dapop,
    wc6 = unique(rd("pccf_wc6point")$PCODE),
    uniq = data.table::data.table(PCODE = uniq$PCODE, da0 = uniq$DA == "0000", recode = rec(uniq)),
    pointdup = unique(rd("pccf_pointdup")$PCODE),
    dup_recode = dups[, list(recode = all(rec(.SD))), by = PCODE, .SDcols = c("DMT", "DMTDIFF")],
    wc5 = unique(rd("pccf_wc5point")$PCODE5), wc4 = unique(rd("pccf_wc4point")$PCODE4),
    wc3 = unique(rd("pccf_wc3point")$PCODE3), wc2 = unique(rd("pccf_wc2point")$PCODE2),
    nonres = egm[EGMRES == "-", PCODE]), class = "pcw_pccf")
}

# ---- Link source -------------------------------------------------------------------------------
pcw_link_source <- function(pcodes, pccf) {
  p <- pcw_clean_pcode(pcodes)
  u <- pccf$uniq[match(p, pccf$uniq$PCODE)]
  dr <- pccf$dup_recode[match(p, pccf$dup_recode$PCODE), recode]
  ls <- data.table::fcase(
    p %chin% pccf$wc6, "C",
    !is.na(u$PCODE) & u$da0, "0",
    !is.na(u$PCODE) & u$recode, "fsa",
    !is.na(u$PCODE), "F",
    p %chin% pccf$pointdup & dr %in% TRUE, "fsa",
    p %chin% pccf$pointdup, "D",
    substr(p, 1, 5) %chin% pccf$wc5, "5",
    substr(p, 1, 4) %chin% pccf$wc4, "4",
    default = "fsa")
  f <- ls == "fsa"
  ls[f] <- data.table::fcase(
    substr(p[f], 1, 3) %chin% pccf$wc3, "3",
    substr(p[f], 1, 2) %chin% pccf$wc2, "2",
    substr(p[f], 1, 1) %chin% c("A", "B", "C", "E", "G", "H", "J", "K", "L", "M", "N", "P", "R", "S", "T", "V", "X", "Y"), "1",
    default = "0")
  ls
}

# ---- Candidates and weights -----------------------------------------------------------------------
pcw_candidates <- function(pcodes, ls, pccf, d_weights = c("equal", "pop")) {
  d_weights <- match.arg(d_weights)
  pc <- unique(data.table::data.table(PCODE = pcw_clean_pcode(pcodes), Link_Source = as.character(ls)))[Link_Source %chin% PCW_KEEP]
  pc[, key := data.table::fcase(Link_Source %chin% c("F", "C", "D"), PCODE, Link_Source == "5", substr(PCODE, 1, 5),
                                Link_Source == "4", substr(PCODE, 1, 4))]
  cp <- pccf$cand[pc, on = c("Link_Source", "key"), allow.cartesian = TRUE, nomatch = NULL]
  if (d_weights == "pop") {
    cp[pccf$dapop, DA21pop := i.DA21pop, on = "DAuid"]
    cp[Link_Source == "D", w_raw := DA21pop / .N, by = list(PCODE, DAuid)]
    cp[, DA21pop := NULL]
  } else cp[Link_Source == "D", w_raw := 1]
  cp[Link_Source == "F", w_raw := 1]
  cp[, cand_id := seq_len(.N), by = PCODE]
  cp[, w0 := data.table::fifelse(is.na(w_raw) | w_raw < 0, 0, w_raw)]
  cp[, w_fallback := sum(w0) <= 0, by = PCODE]
  cp[, w := if (w_fallback[1]) 1 / .N else w0 / sum(w0), by = PCODE]
  cp[, n_cand := .N, by = PCODE]
  xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(cp[, list(lon, lat)], coords = c("lon", "lat"), crs = 4269), PCW_CRS))
  cp[, `:=`(x = xy[, 1], y = xy[, 2], w0 = NULL, key = NULL)]
  data.table::setcolorder(cp, c("PCODE", "Link_Source", "cand_id", "n_cand", "DAuid", "lat", "lon", "x", "y", "w", "w_raw", "w_fallback"))
  cp[]
}

# ---- Exposures, DA attributes and helpers ------------------------------------------------------------
pcw_extract <- function(x, y, exposure, exposure_cols = NULL) {
  pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = PCW_CRS)
  if (inherits(exposure, "SpatRaster")) {
    xy <- sf::st_coordinates(sf::st_transform(pts, sf::st_crs(terra::crs(exposure))))
    cell <- terra::cellFromXY(exposure, xy)
    v <- data.table::as.data.table(stats::setNames(rep(list(rep(NA_real_, length(cell))), terra::nlyr(exposure)), names(exposure)))
    ok <- which(!is.na(cell))
    if (length(ok)) { vals <- exposure[cell[ok]]; for (m in names(v)) data.table::set(v, ok, m, as.numeric(vals[[m]])) }
    v
  } else if (inherits(exposure, "sf")) {
    if (is.null(exposure_cols)) stop("exposure_cols is required for point exposures", call. = FALSE)
    e <- sf::st_transform(exposure, PCW_CRS)
    i <- sf::st_nearest_feature(pts, e)
    v <- data.table::as.data.table(sf::st_drop_geometry(e))
    v[i, exposure_cols, with = FALSE]
  } else stop("exposure must be a terra SpatRaster or an sf POINT layer", call. = FALSE)
}

pcw_quintile_cuts <- function(da_attr, score, quint) {
  d <- data.table::data.table(s = da_attr[[score]], q = as.integer(da_attr[[quint]]))[!is.na(s) & !is.na(q)]
  r <- d[, list(lo = min(s), hi = max(s)), keyby = q]
  if (!identical(r$q, 1:5) || any(r$hi[1:4] >= r$lo[2:5]))
    stop("quintiles of ", quint, " are not separated by ", score, call. = FALSE)
  cut <- (r$hi[1:4] + r$lo[2:5]) / 2
  if (!all(pcw_quint_of(d$s, cut) == d$q)) stop("cut points do not reproduce ", quint, call. = FALSE)
  cut
}
pcw_quint_of <- function(score, cut) ifelse(is.na(score), NA_integer_, 1L + rowSums(outer(score, cut, ">")))

# modal quintile: argmax pi; ties go to the quintile closest to the mean-score quintile, then the lower
pcw_mode_q <- function(P, qa) {
  ok <- !is.na(P[, 1])
  mx <- apply(P, 1, max)
  S <- ifelse(P >= mx - 1e-9, abs(col(P) - ifelse(is.na(qa), 3L, qa)) * 10 + col(P), Inf)
  out <- max.col(-S, ties.method = "first"); out[!ok] <- NA_integer_; out
}

# weighted mean, variance and share of weight missing, by postal code (weights renormalized)
pcw_wmean <- function(cp, e) {
  ok <- !is.na(e)
  d <- data.table::data.table(PCODE = cp$PCODE, wo = cp$w * ok, e0 = data.table::fifelse(ok, e, 0))
  s <- d[, list(sw = sum(wo), swe = sum(wo * e0)), by = PCODE]
  s[, mu := data.table::fifelse(sw > 0, swe / sw, NA_real_)]
  d[s, mu := i.mu, on = "PCODE"]
  v <- d[, list(sv = sum(wo * (e0 - data.table::fcoalesce(mu, 0))^2)), by = PCODE]
  s[v, V := data.table::fifelse(sw > 0, i.sv / sw, NA_real_), on = "PCODE"]
  s[, list(PCODE, mu, V, miss = 1 - sw)]
}

pcw_da_of <- function(x, y, da_boundaries) {
  da <- sf::st_transform(da_boundaries[, "DAUID"], PCW_CRS)
  pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = PCW_CRS)
  hit <- vapply(sf::st_intersects(pts, da), function(v) if (length(v)) v[1] else NA_integer_, 1L)
  near <- is.na(hit); nd <- rep(NA_real_, length(x))
  if (any(near)) {
    hit[near] <- sf::st_nearest_feature(pts[near, ], da)
    nd[near] <- as.numeric(sf::st_distance(pts[near, ], da[hit[near], ], by_element = TRUE))
  }
  list(DAuid = as.character(da$DAUID[hit]), da_method = data.table::fifelse(near, "nearest", "pip"), nearest_m = nd)
}

# ---- Main function --------------------------------------------------------------------------------
pcweight <- function(data, pccf, pcode = "PCODE", link_source = NULL,
                     exposure = NULL, exposure_cols = NULL,
                     da_attr = NULL, da_id = "DAuid", da_scores = NULL, da_quintiles = NULL,
                     da_boundaries = NULL, d_weights = c("equal", "pop"), return_candidates = FALSE) {
  d_weights <- match.arg(d_weights)
  if (is.character(pccf)) pccf <- pcw_read_pccf(pccf)
  stopifnot(inherits(pccf, "pcw_pccf"))
  data <- data.table::as.data.table(data)

  # postal codes and link sources
  pc <- data.table::data.table(PCODE = pcw_clean_pcode(data[[pcode]]))
  if (!is.null(link_source)) {
    pc[, Link_Source := as.character(data[[link_source]])]
    pc <- pc[, list(Link_Source = names(which.max(table(Link_Source))), n_ls = data.table::uniqueN(Link_Source), n_rec = .N), by = PCODE]
    if (any(pc$n_ls > 1)) warning(sum(pc$n_ls > 1), " postal codes have more than one link source; the most frequent is used", call. = FALSE)
    pc[, `:=`(ls_origin = "pccf", n_ls = NULL)]
  } else {
    pc <- pc[, list(n_rec = .N), by = PCODE]
    pc[, `:=`(Link_Source = pcw_link_source(PCODE, pccf), ls_origin = "inferred")]
  }
  pc[, pccf_nonres := Link_Source %chin% c("F", "D") & PCODE %chin% pccf$nonres]

  # candidates and summary
  cp <- pcw_candidates(pc$PCODE, pc$Link_Source, pccf, d_weights)
  ps <- cp[, list(n_cand = .N, n_eff = 1 / sum(w^2), n_da = data.table::uniqueN(DAuid), w_fallback = w_fallback[1]), by = PCODE]
  out <- merge(pc, ps, by = "PCODE", all.x = TRUE, sort = FALSE)
  out[is.na(n_cand), n_cand := 0L]

  # method B centroid
  cen <- cp[, list(B_x = sum(w * x), B_y = sum(w * y)), by = PCODE]
  ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(cen[, list(B_x, B_y)], coords = c("B_x", "B_y"), crs = PCW_CRS), 4269))
  cen[, `:=`(B_lon = ll[, 1], B_lat = ll[, 2], B_DAuid = NA_character_, B_da_method = NA_character_, B_nearest_m = NA_real_)]
  one <- cp[n_cand == 1, list(PCODE, cand_DAuid = DAuid)]
  cen[one, `:=`(B_DAuid = i.cand_DAuid, B_da_method = "candidate"), on = "PCODE"]
  if (!is.null(da_boundaries)) {
    i <- cen[is.na(B_DAuid), which = TRUE]
    if (length(i)) {
      r <- pcw_da_of(cen$B_x[i], cen$B_y[i], da_boundaries)
      cen[i, `:=`(B_DAuid = r$DAuid, B_da_method = r$da_method, B_nearest_m = r$nearest_m)]
    }
  }

  # exposures: method A over candidates, method B at the centroid
  if (!is.null(exposure)) {
    ec <- pcw_extract(cp$x, cp$y, exposure, exposure_cols)
    eb <- pcw_extract(cen$B_x, cen$B_y, exposure, exposure_cols)
    for (m in names(ec)) {
      s <- pcw_wmean(cp, ec[[m]])
      out[s, `:=`(A = i.mu, A_sd = sqrt(i.V), A_miss = i.miss), on = "PCODE"]
      data.table::setnames(out, c("A", "A_sd", "A_miss"), paste0(c("A_", "A_sd_", "A_miss_"), m))
      cen[, (paste0("B_", m)) := eb[[m]]]
    }
  }

  # DA attributes: weighted mean score, quintile, pi, modal quintile; method B at the centroid DA
  if (!is.null(da_attr) && length(da_scores)) {
    da_attr <- data.table::as.data.table(da_attr)
    ida <- match(cp$DAuid, as.character(da_attr[[da_id]]))
    idb <- match(cen$B_DAuid, as.character(da_attr[[da_id]]))
    for (k in names(da_scores)) {
      sc <- da_scores[[k]]
      s <- pcw_wmean(cp, da_attr[[sc]][ida])
      out[s, `:=`(A_S = i.mu, A_sd_S = sqrt(i.V), A_miss_S = i.miss), on = "PCODE"]
      data.table::setnames(out, c("A_S", "A_sd_S", "A_miss_S"), paste0(c("A_S_", "A_sd_S_", "A_miss_S_"), k))
      cen[, (paste0("B_S_", k)) := da_attr[[sc]][idb]]
      if (!is.null(da_quintiles) && k %in% names(da_quintiles)) {
        qc <- da_quintiles[[k]]
        cut <- pcw_quintile_cuts(da_attr, sc, qc)
        out[, (paste0("A_q_", k)) := pcw_quint_of(get(paste0("A_S_", k)), cut)]
        q <- as.integer(da_attr[[qc]][ida]); okq <- !is.na(q)
        d <- data.table::data.table(PCODE = cp$PCODE, wq = cp$w * okq, qq = data.table::fifelse(okq, q, 0L))
        p <- d[, list(sw = sum(wq), p1 = sum(wq * (qq == 1L)), p2 = sum(wq * (qq == 2L)), p3 = sum(wq * (qq == 3L)),
                      p4 = sum(wq * (qq == 4L)), p5 = sum(wq * (qq == 5L))), by = PCODE]
        for (j in 1:5) p[, (paste0("p", j)) := data.table::fifelse(sw > 0, get(paste0("p", j)) / sw, NA_real_)]
        data.table::setnames(p, paste0("p", 1:5), paste0("A_pi_", k, "_", 1:5))
        out <- merge(out, p[, !"sw"], by = "PCODE", all.x = TRUE, sort = FALSE)
        P <- as.matrix(out[, paste0("A_pi_", k, "_", 1:5), with = FALSE])
        out[, (paste0("A_mode_", k)) := pcw_mode_q(P, get(paste0("A_q_", k)))]
        cen[, (paste0("B_q_", k)) := as.integer(da_attr[[qc]][idb])]
      }
    }
  }

  out <- merge(out, cen, by = "PCODE", all.x = TRUE, sort = FALSE)
  data.table::setcolorder(out, c("PCODE", "Link_Source", "ls_origin", "n_rec", "n_cand", "n_eff", "n_da", "w_fallback", "pccf_nonres"))
  data.table::setorder(out, PCODE)
  if (return_candidates) attr(out, "candidates") <- cp
  out[]
}
