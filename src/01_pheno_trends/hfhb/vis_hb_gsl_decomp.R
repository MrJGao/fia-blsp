# ******************************************************************************
# Analyze contributions of spring and autumn phenology on GSL in terms of trend
# and interannual variability for HB data.
# 
# Functions were adapted from the BLSP version.
# 
# Author: Xiaojie Gao
# Date: 2025-03-05
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(lubridate)
library(ggplot2)
library(ggridges)



DecomposeSignal <- function(pheno_dt) {
    # Fit a linear model to get trend
    model <- lm(pheno_dt$pheno ~ pheno_dt$yr)
    # Trend
    slp <- coef(model)[2]

    # Isolate IAV from trend and calculate variance
    trend_estimated <- predict(model)
    # var_trend <- var(trend_estimated)
    detrend <- pheno_dt$pheno - trend_estimated
    var_iav <- var(detrend)

    var_total <- var(pheno_dt$pheno)

    return(list(
        stat = data.table(
            slp, var_iav, var_total
        ),
        detrend = detrend
    ))
}


AnalyzeSiteSpecData <- function(specname) {
    spec_dt <- hb_pheno[SPECIES == specname]

    var_dt <- lapply(unique(spec_dt$SITE), function(ind) {
        ind_dt <- spec_dt[SITE == ind]
        setorder(ind_dt, YEAR)
        ind_dt <- na.omit(ind_dt)

        sos_de <- DecomposeSignal(ind_dt[, .(pheno = spring, yr = YEAR)])
        eos_de <- DecomposeSignal(ind_dt[, .(pheno = fall, yr = YEAR)])
        gsl_de <- DecomposeSignal(ind_dt[, .(pheno = gsl, yr = YEAR)])
        
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

        arow <- data.table(
            spec = specname, SITE = ind,
            rho, covar_pct,
            sos_trend_pct = sos_trend_pct,
            eos_trend_pct = eos_trend_pct,
            sos_iav_pct = sos_iav_pct,
            eos_iav_pct = eos_iav_pct
        )

        return(arow)
    }) %>%
        rbindlist()

    return(var_dt)
}



# ~ Load data ####
# ~ ----------------------------------------------------------------------------
hb_pheno <- fread("data/raw/hb_pheno.csv")

# Linear intepolation
hb_pheno <- by(hb_pheno, list(hb_pheno$SITE, hb_pheno$SPECIES), function(sp) {
    alldays <- seq(min(sp$DATE), max(sp$DATE), by = "day")
    interpo <- approx(sp$DATE, sp$Phenology_Stage, xout = alldays)

    interpo_dt <- data.table(DATE = interpo$x, Phenology_Stage = interpo$y)
    interpo_dt[, DAY := yday(as_date(DATE))]
    interpo_dt[, YEAR := year(DATE)]
    
    spr <- lapply(unique(interpo_dt$YEAR), function(yr) {
        yr_dt <- interpo_dt[YEAR == yr]
        spr <- yr_dt[DAY < 185][which.min(abs(Phenology_Stage - 3)), ]
        return(spr)
    })
    spr <- data.table(do.call(rbind, spr))[, .(
        YEAR, spring = DAY, spring_date = DATE
    )]

    atm <- lapply(unique(interpo_dt$YEAR), function(yr) {
        yr_dt <- interpo_dt[YEAR == yr]
        atm <- yr_dt[DAY > 185][which.min(abs(Phenology_Stage - 1)), ]
        return(atm)
    })
    atm <- data.table(do.call(rbind, atm))[, .(
        YEAR, fall = DAY, fall_date = DATE
    )]
    
    arow <- cbind(spr, atm)

    arow[, gsl := fall - spring]

    arow$SITE <- sp$SITE[1]
    arow$SPECIES <- sp$SPECIES[1]

    
    return(arow)
})
hb_pheno <- do.call(rbind, hb_pheno)

specnames <- unique(hb_pheno$SPECIES)

var_dt <- lapply(specnames, AnalyzeSiteSpecData)
var_dt <- rbindlist(var_dt)

# Get trend significance
# HERE needs to run `vis_hb_pheno_trends.R` to get `gsl_trend_ols`.
var_dt <- merge(
    var_dt, 
    gsl_trend_ols[, .(spec, SITE = site, gsl_slp = slp, gsl_pval = pval)], 
    by = c("spec", "SITE")
)

var_dt[, ":="(
    iav_diff = sos_iav_pct - eos_iav_pct,
    trend_diff = sos_trend_pct - eos_trend_pct
)]

# Calculate some proportions
c_eos_iav_pct <- var_dt[iav_diff < 0, .N / nrow(var_dt)] * 100
c_sos_iav_pct <- 100 - c_eos_iav_pct
c_eos_iav_sig_pct <- var_dt[gsl_pval < 0.05 & iav_diff < 0, .N / nrow(var_dt)] * 100
c_sos_iav_sig_pct <- var_dt[gsl_pval < 0.05 & iav_diff > 0, .N / nrow(var_dt)] * 100

c_eos_trend_pct <- var_dt[trend_diff < 0, .N / nrow(var_dt)] * 100
c_sos_trend_pct <- 100 - c_eos_trend_pct
c_eos_trend_sig_pct <- var_dt[gsl_pval < 0 & trend_diff < 0, .N / nrow(var_dt)] * 100
c_sos_trend_sig_pct <- var_dt[gsl_pval < 0 & trend_diff > 0, .N / nrow(var_dt)] * 100


# Here report the number of trees
var_dt[iav_diff > 0, .N / nrow(var_dt)]
var_dt[
    gsl_pval < 0.05 & iav_diff < 0, 
    .(.N, .N / nrow(var_dt[gsl_pval < 0.05]))
]
var_dt[
    gsl_pval < 0.05 & iav_diff > 0, 
    .(.N, .N / nrow(var_dt[gsl_pval < 0.05]))
]

var_dt[trend_diff < 0, .N / nrow(var_dt)]
var_dt[trend_diff > 0, .N / nrow(var_dt)]
var_dt[
    gsl_pval < 0.05 & trend_diff < 0, 
    .(.N, .N / nrow(var_dt[gsl_pval < 0.05]))
]
var_dt[
    gsl_pval < 0.05 & trend_diff > 0, 
    .(.N, .N / nrow(var_dt[gsl_pval < 0.05]))
]







