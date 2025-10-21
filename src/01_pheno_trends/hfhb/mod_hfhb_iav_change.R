# ******************************************************************************
# Investigate whether spring/fall/gsl IAV is increasing/decreasing w/ climate
# change for HF and HB data.
# 
# Author: Xiaojie Gao
# Date: 2025-05-21
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(lubridate)



# ~ HF ####
# ~ ----------------------------------------------------------------------------
spr_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-05-spring-mean-ind.csv")
atm_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-07-fall-mean-ind.csv")
# com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))

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


hf_var_dt <- lapply(unique(com_dt$tree.id), function(ind) {
    ind_dt <- com_dt[tree.id == ind]
    setorder(ind_dt, year)
    ind_dt <- na.omit(ind_dt)

    if (nrow(ind_dt) < 20) {
        return(NULL)
    }

    sos_de <- DecomposeSignal(ind_dt[, .(pheno = l75.doy, yr = year)])
    sos_iav_trend <- CalIAVTrend(sos_de$detrend, ind_dt$year)

    eos_de <- DecomposeSignal(ind_dt[, .(pheno = lf.doy, yr = year)])
    eos_iav_trend <- CalIAVTrend(eos_de$detrend, ind_dt$year)

    gsl_de <- DecomposeSignal(ind_dt[, .(pheno = gsl, yr = year)])
    gsl_iav_trend <- CalIAVTrend(gsl_de$detrend, ind_dt$year)

    arow <- data.table(
        spec = ind_dt$species[1], tree.id = ind,
        sos_iav_slp = sos_iav_trend$slp,
        sos_iav_pval = sos_iav_trend$pval,
        eos_iav_slp = eos_iav_trend$slp,
        eos_iav_pval = eos_iav_trend$pval,
        gsl_iav_slp = gsl_iav_trend$slp,
        gsl_iav_pval = gsl_iav_trend$pval
    )

    return(arow)
}) %>%
    rbindlist()



# ~ HB ####
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
        YEAR,
        spring = DAY, spring_date = DATE
    )]

    atm <- lapply(unique(interpo_dt$YEAR), function(yr) {
        yr_dt <- interpo_dt[YEAR == yr]
        atm <- yr_dt[DAY > 185][which.min(abs(Phenology_Stage - 1)), ]
        return(atm)
    })
    atm <- data.table(do.call(rbind, atm))[, .(
        YEAR,
        fall = DAY, fall_date = DATE
    )]

    arow <- cbind(spr, atm)

    arow[, gsl := fall - spring]

    arow$SITE <- sp$SITE[1]
    arow$SPECIES <- sp$SPECIES[1]


    return(arow)
})
hb_pheno <- do.call(rbind, hb_pheno)

specnames <- unique(hb_pheno$SPECIES)

hb_var_dt <- lapply(specnames, function(specname) {
    spec_dt <- hb_pheno[SPECIES == specname]

    var_dt <- lapply(unique(spec_dt$SITE), function(ind) {
        ind_dt <- spec_dt[SITE == ind]
        setorder(ind_dt, YEAR)
        ind_dt <- na.omit(ind_dt)

        if (nrow(ind_dt) < 20) {
            return(NULL)
        }

        sos_de <- DecomposeSignal(ind_dt[, .(pheno = spring, yr = YEAR)])
        sos_iav_trend <- CalIAVTrend(sos_de$detrend, ind_dt$YEAR)

        eos_de <- DecomposeSignal(ind_dt[, .(pheno = fall, yr = YEAR)])
        eos_iav_trend <- CalIAVTrend(eos_de$detrend, ind_dt$YEAR)

        gsl_de <- DecomposeSignal(ind_dt[, .(pheno = gsl, yr = YEAR)])
        gsl_iav_trend <- CalIAVTrend(gsl_de$detrend, ind_dt$YEAR)


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
            sos_iav_slp = sos_iav_trend$slp,
            sos_iav_pval = sos_iav_trend$pval,
            eos_iav_slp = eos_iav_trend$slp,
            eos_iav_pval = eos_iav_trend$pval,
            gsl_iav_slp = gsl_iav_trend$slp,
            gsl_iav_pval = gsl_iav_trend$pval
        )

        return(arow)
    }) %>%
        rbindlist()

    return(var_dt)
})
hb_var_dt <- rbindlist(hb_var_dt)



# ~ Combine HF & HB ####
# ~ ----------------------------------------------------------------------------
setnames(hf_var_dt, "tree.id", "ID")
hf_var_dt[, site := "HF"]

setnames(hb_var_dt, "SITE", "ID")
hb_var_dt[, ID := paste0(spec, "-", ID)]
hb_var_dt[, site := "HB"]

hfhb_var_dt <- rbind(hf_var_dt, hb_var_dt)

# out:
fwrite(hfhb_var_dt, "pipe/species_pheno/hfhb_var_dt.csv")
