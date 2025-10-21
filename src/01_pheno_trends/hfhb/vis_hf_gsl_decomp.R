# ******************************************************************************
# Analyze contributions of spring and autumn phenology on GSL in terms of trend
# and interannual variability for HF data.
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
    spec_dt <- com_dt[species == specname]

    var_dt <- lapply(unique(spec_dt$tree.id), function(ind) {
        ind_dt <- spec_dt[tree.id == ind]
        setorder(ind_dt, year)
        ind_dt <- na.omit(ind_dt)

        sos_de <- DecomposeSignal(ind_dt[, .(pheno = l75.doy, yr = year)])
        eos_de <- DecomposeSignal(ind_dt[, .(pheno = lf.doy, yr = year)])
        gsl_de <- DecomposeSignal(ind_dt[, .(pheno = gsl, yr = year)])

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
            spec = specname, tree.id = ind,
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
spr_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-05-spring-mean-ind.csv")
atm_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-07-fall-mean-ind.csv")
com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))

# Select only trees w/ longer than 10 years of data
spr_tree_records <- spr_dt[, .(num_records = .N), by = tree.id]
spr_dt <- spr_dt[tree.id %in% spr_tree_records[num_records > 10, tree.id]]

atm_tree_records <- atm_dt[, .(num_records = .N), by = tree.id]
atm_dt <- atm_dt[tree.id %in% atm_tree_records[num_records > 10, tree.id]]

com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))
# Remove the trees that only have fall phenology
com_dt <- na.omit(com_dt[, .(year, tree.id, species, l75.doy, lf.doy)])
com_dt[, gsl := lf.doy - l75.doy]

# Select only deciduous tree species
deci_spec <- c(
    "ACPE", "ACRU", "ACSA", "AMSP", "BEAL", "BELE", "BEPA", "BEPO", "CADE",
    "FAGR", "FRAM", "NYSY", "POTR", "PRSE",
    "QUAL", "QURU", "QUVE"
)
com_dt <- com_dt[species %in% deci_spec, ]

uniqueN(com_dt$species)
uniqueN(com_dt$tree.id)

specnames <- unique(com_dt$species)
var_dt <- lapply(specnames, AnalyzeSiteSpecData)
var_dt <- rbindlist(var_dt)

# Get trend significance 
# HACK: HERE needs to run `DoTrendOLS()` in `vis_hf_pheno_trends.R` to get
# `gsl_trend_ols`.
var_dt <- merge(
    var_dt,
    gsl_trend_ols[, .(spec, tree.id, gsl_slp = slp, gsl_pval = pval)],
    by = c("spec", "tree.id")
)

var_dt[, ":="(
    iav_diff = sos_iav_pct - eos_iav_pct,
    trend_diff = sos_trend_pct - eos_trend_pct
)]

# Calculate some proportions
var_dt[iav_diff < 0, .N / nrow(var_dt)]
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
