# ******************************************************************************
# Bootstrap BLSP based on posterior distributions and then re-run the trend
# analysis to see how much uncertain the result can be.
# 
# Author: Xiaojie Gao
# Date: 2026-06-22
# ******************************************************************************
rm(list = ls())

source("src/base.R")
library(data.table)
library(lubridate)
library(trend)
library(parallel)



cachedir <- "pipe/blsp_trend_bootstrap"
dir.create(cachedir, showWarnings = FALSE, recursive = TRUE)


# ~ Site-specific trend ####
# ~ ----------------------------------------------------------------------------
SamplePosterior <- function(ph, lwr, upr) {
    if (is.na(upr) || is.na(lwr)) {
        stop("something is wrong...")
    }
    if (upr - lwr > 30) {
        return(NA)
    }
    z <- qnorm(1 - (1 - 0.95) / 2)
    sd <- (upr - lwr) / (2 * z)

    return(rnorm(1, ph, sd))
}

SamplePhenoFromPosterior <- function(phenos) {
    # Sample phenos from the posterior distribution
    sos <- sapply(1:nrow(phenos), \(i) {
        if (is.na(phenos[i, MidGreenup])) {
            return(NA)
        } else {
            return(SamplePosterior(
                phenos[i, MidGreenup],
                lwr = phenos[i, MidGreenup_lwr],
                upr = phenos[i, MidGreenup_upr]
            ))
        }
    })
    eos <- sapply(1:nrow(phenos), \(i) {
        if (is.na(phenos[i, MidGreendown])) {
            return(NA)
        } else {
            return(SamplePosterior(
                phenos[i, MidGreendown],
                lwr = phenos[i, MidGreendown_lwr],
                upr = phenos[i, MidGreendown_upr]
            ))
        }
    })

    gsl <- eos - sos
    
    sampled_phenos <- data.table(
        Year = phenos$Year,
        SOS = sos, SOS_lwr = phenos$MidGreenup_lwr, SOS_upr = phenos$MidGreenup_upr,
        EOS = eos, EOS_lwr = phenos$MidGreendown_lwr, EOS_upr = phenos$MidGreendown_upr,
        GSL = gsl
    )
    
    return(sampled_phenos)
}

CalTrend <- function(pheno) {
    # ~ Spring ----------------------------------------
    # Sample phenos from the posterior distribution
    sos_noNA <- na.omit(pheno$SOS)
    if (length(sos_noNA) > 10) {
        spr_pval <- mk.test(sos_noNA)$p.value
        spr_slp <- sens.slope(sos_noNA)$estimates
    } else {
        spr_pval <- NA
        spr_slp <- NA
    }


    # ~ Autumn ----------------------------------------
    eos_noNA <- na.omit(pheno$EOS)
    if (length(eos_noNA) > 10) {
        atm_pval <- mk.test(eos_noNA)$p.value
        atm_slp <- sens.slope(eos_noNA)$estimates
    } else {
        atm_pval <- NA
        atm_slp <- NA
    }


    # ~ Growing season length ----------------------------------------
    gsl_noNA <- na.omit(pheno$GSL)
    if (length(gsl_noNA) > 10) {
        gsl_pval <- mk.test(gsl_noNA)$p.value
        gsl_slp <- sens.slope(gsl_noNA)$estimates
    } else {
        gsl_pval <- NA
        gsl_slp <- NA
    }

    res <- data.table(
        spr_slp = spr_slp, spr_pval = spr_pval,
        atm_slp = atm_slp, atm_pval = atm_pval,
        gsl_slp = gsl_slp, gsl_pval = gsl_pval
    )

    return(res)
}

DecomposeSignalCorrected <- function(pheno_dt, metric) {
    # Fit a linear model to get trend
    model <- lm(pheno_dt[[metric]] ~ pheno_dt$yr)
    # Trend
    slp <- coef(model)[2]

    # Isolate IAV from trend and calculate variance
    trend_estimated <- predict(model)
    # var_trend <- var(trend_estimated)
    detrend <- pheno_dt[[metric]] - trend_estimated
    
    # Uncertainty corrected variance
    if (metric %in% c("SOS", "EOS")) {
        sigma2 <- mean(
            (pheno_dt[[paste0(metric, "_upr")]] - 
                pheno_dt[[paste0(metric, "_lwr")]]) / 3.92, 
            na.rm = TRUE
        )^2
        var_iav <- max(var(detrend) - sigma2, 0)
    } else if (metric == "GSL") {
        # Here we have to assume SOS and EOS are independent b/c we don't have
        # their true covariance that is not affected by their estimate uncertainty
        SOS_sigma2 <- mean(
            (pheno_dt[["SOS_upr"]] - pheno_dt[["SOS_lwr"]]) / 3.92,
            na.rm = TRUE
        )^2
        EOS_sigma2 <- mean(
            (pheno_dt[["EOS_upr"]] - pheno_dt[["EOS_lwr"]]) / 3.92,
            na.rm = TRUE
        )^2
        GSL_sigma2 <- SOS_sigma2 + EOS_sigma2
        var_iav <- max(var(detrend) - GSL_sigma2)
    }

    return(list(
        stat = data.table(
            slp, var_iav
        ),
        detrend = detrend
    ))
}

CalContribution <- function(pheno) {
    pheno <- na.omit(pheno)
    pheno[, yr := as.numeric(Year)]

    if (nrow(pheno) < 20) {
        return(NULL)
    }

    sos_de <- DecomposeSignalCorrected(pheno, "SOS")
    eos_de <- DecomposeSignalCorrected(pheno, "EOS")
    gsl_de <- DecomposeSignalCorrected(pheno, "GSL")

    # Calculate trend cotribution
    sos_trend_pct <- -sos_de$stat$slp / gsl_de$stat$slp * 100
    eos_trend_pct <- eos_de$stat$slp / gsl_de$stat$slp * 100

    # Calculate IAV contribution
    rho <- cor(sos_de$detrend, eos_de$detrend)
    covar <- rho * sqrt(sos_de$stat$var_iav) * sqrt(eos_de$stat$var_iav)
    sos_iav_pct <- (sos_de$stat$var_iav - covar) / gsl_de$stat$var_iav * 100
    eos_iav_pct <- (eos_de$stat$var_iav - covar) / gsl_de$stat$var_iav * 100

    # Correlation contribution
    covar_pct <- covar / gsl_de$stat$var_iav * 100

    res <- data.table(
        rho, covar_pct,
        sos_ols_slp = sos_de$stat$slp,
        sos_ols_iav = sos_de$stat$var_iav,
        eos_ols_slp = eos_de$stat$slp,
        eos_ols_iav = eos_de$stat$var_iav,
        gsl_ols_slp = gsl_de$stat$slp,
        gsl_ols_iav = gsl_de$stat$var_iav,
        
        sos_trend_pct = sos_trend_pct,
        eos_trend_pct = eos_trend_pct,
        sos_iav_pct = sos_iav_pct,
        eos_iav_pct = eos_iav_pct
    )

    return(res)
}


BootstrapAnalysis <- function(specname) {
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
    bs_dt <- lapply(spec_fit, function(thefit) {
        phenos <- thefit$phenos
        if (nrow(phenos) < length(yrs) || ncol(phenos) < 20) {
            return(NULL)
        }
        setorder(phenos, Year)
        
        samp_phenos <- SamplePhenoFromPosterior(phenos)

        trend_res <- CalTrend(samp_phenos)
        contri_res <- CalContribution(samp_phenos)

        if (is.null(trend_res) || is.null(contri_res)) {
            arow <- NULL
        } else {
            arow <- cbind(spec = specname, ID = thefit$ID, trend_res, contri_res)
        }
        
        return(arow)
    })
    bs_dt <- do.call(rbind, bs_dt)
    
    return(bs_dt)
}


cl <- makeCluster(20)
calls <- clusterCall(cl, function() {
    suppressPackageStartupMessages({
        source("src/base.R")
        library(data.table)
        library(lubridate)
        library(trend)
    })
})

clusterExport(cl, c(
    "BootstrapAnalysis", "cachedir", 
    "SamplePosterior", "SamplePhenoFromPosterior", 
    "CalTrend", "DecomposeSignalCorrected", "CalContribution"
))
clusterApplyLB(cl, x = 1:1000, fun = function(i) {
    # Get deci/ever info
    fia_domspec_sites <- fread("data/fia_domspec_truexy.csv")
    fia_domspec_sites <- unique(fia_domspec_sites[, .(spec, CommonName, Type)])
    deci_sites <- fia_domspec_sites[Type == "deci", ]

    site_dt <- NULL
    for (specname in deci_sites[, spec]) {
        site_dt <- rbind(site_dt, BootstrapAnalysis(specname))
    }
    site_dt <- merge(site_dt, fia_domspec_sites, by = "spec")

    # out:
    fwrite(site_dt, file.path(cachedir, paste0("deci_blsp_trend_", i, ".csv")))
})
stopCluster(cl)


