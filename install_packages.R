# Installs the R packages used by this repository (developed with R 4.6.1).
# The portable function R/pcweight.R needs only data.table and sf (and terra for raster
# exposures). The pccf_nar_validation pipeline has its own install_packages.R.
pkgs <- c("data.table", "sf", "terra", "ncdf4", "readxl", "writexl", "ggplot2", "patchwork", "scales")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
message("All packages installed: ", all(pkgs %in% rownames(installed.packages())))
