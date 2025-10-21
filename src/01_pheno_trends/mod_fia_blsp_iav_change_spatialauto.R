# ******************************************************************************
# Spatial autocorrelation of changes in IAV
# 
# Author: Xiaojie Gao
# Date: 2025-09-21
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(terra)
library(tmap)
library(spdep)


if (!exists("var_dt")) {
    source("src/01_pheno_trends/mod_fia_blsp_iav_change.R")
}


# Get lat and lon
fia_site_loc <- fread("data/fia_blsp_fit.csv")
fia_site_loc <- unique(fia_site_loc[, .(ID, x, y)])
var_dt <- merge(var_dt, fia_site_loc, by = "ID")
var_dt_vec <- vect(var_dt, geom = c("x", "y"), crs = "epsg:4326")

specnames <- unique(var_dt$spec) %>%
    sort()

# Moran's I test for spatial autocorrelation
moranI_dt <- lapply(specnames, function(sp) {
    sp_dt <- var_dt[spec == sp]

    dists <- as.matrix(dist(sp_dt[, .(x, y)]))
    inv_dists <- 1 / dists
    diag(inv_dists) <- 0
    lw <- mat2listw(inv_dists, style = "W")
    
    # SOS
    moran_result <- moran.test(sp_dt$sos_iav_slp, lw)
    sos_moran <- moran_result$statistic
    sos_moran_pval <- moran_result$p.value
    # EOS
    moran_result <- moran.test(sp_dt$eos_iav_slp, lw)
    eos_moran <- moran_result$statistic
    eos_moran_pval <- moran_result$p.value
    # GSL
    moran_result <- moran.test(sp_dt$gsl_iav_slp, lw)
    gsl_moran <- moran_result$statistic
    gsl_moran_pval <- moran_result$p.value

    arow <- data.table(
        spec = sp, 
        sos_moran = sos_moran, 
        sos_moran_pval = sos_moran_pval, 
        eos_moran = eos_moran, 
        eos_moran_pval = eos_moran_pval, 
        gsl_moran = gsl_moran,
        gsl_moran_pval = gsl_moran_pval
    )

    return(arow)
})
moranI_dt <- rbindlist(moranI_dt)

# out: 
fwrite(moranI_dt, "out/iav_change_moran_test.csv")



