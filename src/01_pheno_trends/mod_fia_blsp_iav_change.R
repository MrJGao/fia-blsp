# ******************************************************************************
# Investigate whether spring/fall/gsl IAV is increasing/decreasing w/ climate
# change.
# 
# Author: Xiaojie Gao
# Date: 2025-06-20
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(lubridate)



AnalyzeSiteSpecData <- function(specname) {
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

    var_dt <- lapply(spec_fit, function(thefit) {
        phenos <- thefit$phenos
        if (nrow(phenos) < length(yrs) || ncol(phenos) < 22) {
            return(NULL)
        }
        setorder(phenos, Year)

        phenos[MidGreenup_upr - MidGreenup_lwr > 30, MidGreenup := NA]
        phenos[MidGreendown_upr - MidGreendown_lwr > 30, MidGreenup := NA]

        phenos <- na.omit(phenos[, .(MidGreenup, MidGreendown, yr = Year)])
        phenos[, yr := as.numeric(yr)]

        if (nrow(phenos) < 15) {
            return(NULL)
        }

        sos_de <- DecomposeSignal(phenos[, .(pheno = MidGreenup, yr)])
        sos_iav_trend <- CalIAVTrend(sos_de$detrend, phenos$yr)

        eos_de <- DecomposeSignal(phenos[, .(pheno = MidGreendown, yr)])
        eos_iav_trend <- CalIAVTrend(eos_de$detrend, phenos$yr)

        gsl_de <- DecomposeSignal(
            phenos[, .(pheno = MidGreendown - MidGreenup, yr)]
        )
        gsl_iav_trend <- CalIAVTrend(gsl_de$detrend, phenos$yr)


        # Get latitude
        lat <- evi2_dt[ID == thefit$ID, Latitude][1]


        arow <- data.table(
            spec = specname, ID = thefit$ID, lat = lat,
            sos_iav_slp = sos_iav_trend$slp,
            sos_iav_pval = sos_iav_trend$pval,
            eos_iav_slp = eos_iav_trend$slp,
            eos_iav_pval = eos_iav_trend$pval,
            gsl_iav_slp = gsl_iav_trend$slp,
            gsl_iav_pval = gsl_iav_trend$pval
        )

        return(arow)
    })
    var_dt <- rbindlist(var_dt)

    return(var_dt)
}




# ~ Main -----------------------------------------------------------------------
evi2_dt <- rbind(
    fread("data/fia_dom_landsat.csv"),
    fread("data/fia_dom_landsat_more.csv")
)
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
majorspec_dt <- fread("data/fia_domspec_truexy.csv")
major_ids <- majorspec_dt[, paste0("X", X)]
evi2_dt <- evi2_dt[ID %in% major_ids]

specnames <- majorspec_dt[Type == "deci", spec] %>%
    unique()


var_dt <- lapply(specnames, AnalyzeSiteSpecData)
var_dt <- rbindlist(var_dt)





