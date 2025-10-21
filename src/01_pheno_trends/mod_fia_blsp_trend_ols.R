# ******************************************************************************
# Trend analysis on FIA BLSP time series. Same as the `mod_fia_blsp_trend.R` but
# with ordinary linear regression instead of Mann-Kendall and Theil-Sen slope.
# 
# Author: Xiaojie Gao
# Date: 2024-10-03
# ******************************************************************************
rm(list = ls())

source("src/base.R")
library(data.table)
library(lubridate)




# ~ Site-specific trend ####
# ~ ----------------------------------------------------------------------------
SiteSpecTrend <- function(specname) {
    # Read the species file
    spec_fit <- readRDS(file.path("pipe/blsp_fit", paste0(specname, ".Rds")))

    # Use the EVI2 dt to get the latitude for each site
    evi2_dt <- rbind(
        fread("data/fia_dom_landsat.csv"), 
        fread("data/fia_dom_landsat_more.csv")
    )
    # B/c I updated the major species site table, here I need to filter out the
    # sites that were removed and add sites that were added!
    majorspec_dt <- fread("data/fia_domspec_truexy.csv")
    major_ids <- majorspec_dt[, paste0("X", X)]
    evi2_dt <- evi2_dt[ID %in% major_ids]

    yrs <- spec_fit[[1]]$phenos[, as.numeric(Year)]

    # Iterate all the fit, use a linear regression to get a slope
    trend_dt <- lapply(spec_fit, function(thefit) {
        phenos <- thefit$phenos
        if (nrow(phenos) < length(yrs) || ncol(phenos) < 22) {
            return(NULL)
        }
        setorder(phenos, Year)

        # ~ Spring ----------------------------------------
        sos <- phenos[
            MidGreenup_upr - MidGreenup_lwr < 30, 
            .(Year = as.numeric(Year), MidGreenup)
        ]
        sos_noNA <- na.omit(sos)
        if (nrow(sos_noNA) > 10) {
            lmfit <- base$FitLm(as.numeric(sos_noNA$Year), sos_noNA$MidGreenup)
            spr_pval <- lmfit$pval
            spr_slp <- lmfit$cf[2]
        } else {
            spr_pval <- NA
            spr_slp <- NA
        }
        

        # ~ Autumn ----------------------------------------
        eos <- phenos[
            MidGreendown_upr - MidGreendown_lwr < 30,
            .(Year = as.numeric(Year), MidGreendown)
        ]
        eos_noNA <- na.omit(eos)
        if (nrow(eos_noNA) > 10) {
            lmfit <- base$FitLm(as.numeric(eos_noNA$Year), eos_noNA$MidGreendown)
            atm_pval <- lmfit$pval
            atm_slp <- lmfit$cf[2]
        } else {
            atm_pval <- NA
            atm_slp <- NA
        }


        # ~ Growing season length ----------------------------------------
        gsl <- merge(sos_noNA, eos_noNA, by = "Year")
        gsl <- gsl[, gsl := MidGreendown - MidGreenup]
        gsl_noNA <- na.omit(gsl)
        if (nrow(gsl_noNA) > 10) {
            lmfit <- base$FitLm(as.numeric(gsl_noNA$Year), gsl_noNA$gsl)
            gsl_pval <- lmfit$pval
            gsl_slp <- lmfit$cf[2]
        } else {
            gsl_pval <- NA
            gsl_slp <- NA
        }

        # Get latitude
        lat <- evi2_dt[ID == thefit$ID, Latitude][1]

        arow <- data.table(
            spec = specname, ID = thefit$ID, lat = lat,
            spr_slp = spr_slp, spr_pval = spr_pval,
            atm_slp = atm_slp, atm_pval = atm_pval,
            gsl_slp = gsl_slp, gsl_pval = gsl_pval
        )
        return(arow)
    })
    trend_dt <- do.call(rbind, trend_dt)
    
    return(trend_dt)
}



evi2_dt <- rbind(
    fread("data/fia_dom_landsat.csv"), 
    fread("data/fia_dom_landsat_more.csv")
)
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
majorspec_dt <- fread("data/fia_domspec_truexy.csv")
major_ids <- majorspec_dt[, paste0("X", X)]
evi2_dt <- evi2_dt[ID %in% major_ids]

specnames <- unique(evi2_dt$Category)
specnames <- specnames[specnames != "PITA"]
site_trend_dt <- NULL
for (specname in specnames) {
    site_trend_dt <- rbind(site_trend_dt, SiteSpecTrend(specname))
}
site_trend_dt <- merge(site_trend_dt, 
    unique(evi2_dt[, .(ID, lon = Longitude)]), 
    by = "ID"
)
# Get deci/ever info
fia_domspec_sites <- fread("data/fia_domspec_truexy.csv")
fia_domspec_sites <- fia_domspec_sites[, .(spec, CommonName, Type)] %>%
    unique()
site_trend_dt <- merge(site_trend_dt, fia_domspec_sites, by = "spec")

fwrite(site_trend_dt, "pipe/fia_dom_blsp_trend_ols.csv")


