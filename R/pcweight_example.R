# =============================================================================
# pcweight_example.R  Example use of pcweight.R
#
# Assigns greenness (EVI) and the 2021 Material and Social Deprivation Index (MSDI) to a
# file of records with six-character postal codes, by method A (population-weighted mean
# over the PCCF+ candidate locations) and method B (population-weighted centroid).
#
# You need:
#   1. a licensed PCCF+ 8D installation (Statistics Canada); its DATA_R folder holds the
#      .rds files read by pcw_read_pccf(). PCCF+ files are not distributed with this code.
#   2. a CSV of records with an identifier and a postal code (and, optionally, the
#      Link_Source column of a PCCF+ run; otherwise the link source is inferred).
#   3. an exposure raster (any CRS; values are taken from the containing cell), e.g. a
#      growing-season EVI GeoTIFF; or an sf point layer (values from the nearest point).
#   4. optionally, the MSDI 2021 equivalence table (INSPQ) and the 2021 DA boundary file
#      (Statistics Canada, lda_000b21a_e) for method B DAs.
# Edit the paths below, then: Rscript R/pcweight_example.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(sf); library(terra); library(readxl) })
source("R/pcweight.R")

PCCF_DIR  <- "path/to/PCCF8D/DATA_R"
RECORDS   <- "path/to/records.csv"                      # columns ID, PCODE
EVI_TIF   <- "path/to/evi_growing_season.tif"
MSDI_XLSX <- "path/to/1. EquivalenceTableCanada2021_en.xlsx"
DA_SHP    <- "path/to/lda_000b21a_e.shp"
OUT_CSV   <- "pcweight_postal_codes.csv"

pccf <- pcw_read_pccf(PCCF_DIR)                         # about 20 s; reuse for several runs
rec  <- fread(RECORDS, colClasses = "character")

# exposure layers: the pixel value and the 3 x 3 window mean of non-missing cells, as in the
# validation study (a cell with a missing value still gets the mean of its neighbours)
evi <- rast(EVI_TIF)
evi <- c(evi, focal(evi, w = 3, fun = "mean", na.rm = TRUE, na.policy = "all"))
names(evi) <- c("evi_pix", "evi_3x3")

msdi <- as.data.table(read_excel(MSDI_XLSX, sheet = "Data"))[
  , .(DAuid = as.character(DA), SCOREMAT, SCORESOC, QuintMat = as.integer(QuintMat), QuintSoc = as.integer(QuintSoc))]
da <- st_read(DA_SHP, quiet = TRUE)[, "DAUID"]

pc <- pcweight(rec, pccf, pcode = "PCODE", exposure = evi,
               da_attr = msdi, da_scores = c(mat = "SCOREMAT", soc = "SCORESOC"),
               da_quintiles = c(mat = "QuintMat", soc = "QuintSoc"),
               da_boundaries = da, d_weights = "equal")

# one row per postal code:
#   Link_Source, n_cand, n_eff, n_da        candidates and their spread
#   A_evi_3x3, A_sd_evi_3x3, A_miss_evi_3x3 method A greenness, SD over candidates, share missing
#   A_S_mat, A_q_mat, A_mode_mat            method A material score, its quintile, modal quintile
#   A_pi_mat_1 .. A_pi_mat_5                probability of each quintile under PCCF+ allocation
#   B_x, B_y, B_lat, B_lon, B_DAuid         method B centroid and its DA
#   B_evi_3x3, B_S_mat, B_q_mat             method B values
# postal codes resolved by PCCF+ from three or fewer characters (Link_Source 3, 2, 1, 0) have
# no values.
fwrite(pc, OUT_CSV)
print(pc[, .N, keyby = Link_Source])

# attach to the records
rec[, PCODE_clean := pcw_clean_pcode(PCODE)]
rec <- pc[, .(PCODE, Link_Source, evi_A = A_evi_3x3, evi_B = B_evi_3x3, mat_q_mode = A_mode_mat, soc_q_mode = A_mode_soc)][
  rec, on = .(PCODE = PCODE_clean)]
