# ******************************************************************************
# Check the robustness of the spring/fall/gsl IAV stationarity.
# 
# Author: Xiaojie Gao
# Date: 2026-06-20
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(lubridate)



# Check if the IAV stationarity was correlated with BLSP CI
CalStationarity_WithCI <- function(specname) {
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

        phenos[MidGreenup_upr - MidGreenup_lwr > 30, MidGreenup := NA]
        phenos[MidGreendown_upr - MidGreendown_lwr > 30, MidGreendown := NA]
        phenos[, yr := as.numeric(Year)]
        setorder(phenos, yr)

        phenos <- na.omit(phenos)
        if (nrow(phenos) < 20) {
            return(NULL)
        }

        sos_de <- DecomposeSignal(phenos[, .(pheno = MidGreenup, yr)])
        sos_iav_trend <- CalIAVTrend(sos_de$detrend, phenos$yr)
        phenos[, sos_ci := MidGreenup_upr - MidGreenup_lwr]
        sos_ci_trend <- CalIAVTrend(phenos$sos_ci, phenos$yr)

        eos_de <- DecomposeSignal(phenos[, .(pheno = MidGreendown, yr)])
        eos_iav_trend <- CalIAVTrend(eos_de$detrend, phenos$yr)
        phenos[, eos_ci := MidGreendown_upr - MidGreendown_lwr]
        eos_ci_trend <- CalIAVTrend(phenos$eos_ci, phenos$yr)

        arow <- data.table(
            spec = specname, ID = thefit$ID,
            sos_iav_slp = sos_iav_trend$slp,
            sos_iav_pval = sos_iav_trend$pval,
            sos_ci_slp = sos_ci_trend$slp,
            sos_ci_pval = sos_ci_trend$pval,

            eos_iav_slp = eos_iav_trend$slp,
            eos_iav_pval = eos_iav_trend$pval,
            eos_ci_slp = eos_ci_trend$slp,
            eos_ci_pval = eos_ci_trend$pval
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


iav_ci_var_dt <- NULL
for (spec in specnames) {
    spec_var_dt <- CalStationarity_WithCI(spec)
    iav_ci_var_dt <- rbind(iav_ci_var_dt, spec_var_dt)
}



{
    png(
        "out/species_phenology/blsp_iav_ci_corr.png",
        width = 1200, height = 600, res = 200
    )
    par(mfrow = c(1, 2), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 1, 1))

    plot(
        iav_ci_var_dt[, .(sos_ci_slp, sos_iav_slp)],
        pch = 16, cex = 0.3, col = "seagreen",
        ylim = c(-10, 10),
        xlab = "ΔCI SOS", ylab = "ΔIAV SOS"
    )
    abline(v = 0, lty = 2, col = "grey")
    abline(h = 0, lty = 2, col = "grey")
    cortest <- cor.test(
        iav_ci_var_dt[sos_iav_pval < 0.05]$sos_ci_slp,
        iav_ci_var_dt[sos_iav_pval < 0.05]$sos_iav_slp
    )
    legend("bottomright", bty = "n", legend = c(
        bquote(r == .(round(cortest$estimate, 3))),
        bquote(p == .(round(cortest$p.value, 3)))
    ))

    plot(
        iav_ci_var_dt[, .(eos_ci_slp, eos_iav_slp)],
        pch = 16, cex = 0.3, col = "orange",
        ylim = c(-10, 10),
        xlab = "ΔCI EOS", ylab = "ΔIAV EOS"
    )
    abline(v = 0, lty = 2, col = "grey")
    abline(h = 0, lty = 2, col = "grey")
    cortest <- cor.test(
        iav_ci_var_dt[eos_iav_pval < 0.05]$eos_ci_slp,
        iav_ci_var_dt[eos_iav_pval < 0.05]$eos_iav_slp
    )
    legend("bottomright", bty = "n", legend = c(
        bquote(r == .(round(cortest$estimate, 3))),
        bquote(p == .(round(cortest$p.value, 3)))
    ))

    dev.off()
}

