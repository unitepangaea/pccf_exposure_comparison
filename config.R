# =============================================================================
# Paths used by every script. Each script finds the repository root (the folder
# containing this file) and sources it, so no script contains a machine-specific path.
#
# This repository extends the PCCF+/NAR validation study. Its scripts read the outputs of
# that pipeline and reuse functions from three of its scripts, so a checkout of
# https://github.com/unitepangaea/pccf_nar_validation, with its pipeline run, is required.
#
#   VALIDATION   the pccf_nar_validation checkout (default: a sibling folder)
#   RAW          public input data          } taken from the configuration of the
#   PCCF_DIR     licensed PCCF+ 8D files    } validation checkout unless set in
#   VAL_DERIVED  validation pipeline output } config_local.R
#   DERIVED      outputs of this repository (data/derived)
#
# To change any path, copy config_local_example.R to config_local.R and edit it.
# config_local.R is not tracked by git.
# =============================================================================

VALIDATION <- file.path(dirname(ROOT), "pccf_nar_validation")
DERIVED    <- file.path(ROOT, "data", "derived")
if (file.exists(file.path(ROOT, "config_local.R"))) source(file.path(ROOT, "config_local.R"))

VALIDATION <- normalizePath(path.expand(VALIDATION), mustWork = FALSE)
if (!file.exists(file.path(VALIDATION, "config.R")))
  stop("pccf_nar_validation checkout not found at ", VALIDATION,
       "; clone it there or set VALIDATION in config_local.R", call. = FALSE)
local({
  v <- new.env(); v$ROOT <- VALIDATION
  sys.source(file.path(VALIDATION, "config.R"), envir = v)
  if (!exists("RAW", envir = globalenv(), inherits = FALSE))         assign("RAW", v$RAW, envir = globalenv())
  if (!exists("PCCF_DIR", envir = globalenv(), inherits = FALSE))    assign("PCCF_DIR", v$PCCF_DIR, envir = globalenv())
  if (!exists("VAL_DERIVED", envir = globalenv(), inherits = FALSE)) assign("VAL_DERIVED", v$DERIVED, envir = globalenv())
})
RAW <- path.expand(RAW); PCCF_DIR <- path.expand(PCCF_DIR); VAL_DERIVED <- path.expand(VAL_DERIVED)
DERIVED <- path.expand(DERIVED)
dir.create(DERIVED, recursive = TRUE, showWarnings = FALSE)
