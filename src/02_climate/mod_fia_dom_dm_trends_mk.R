# ******************************************************************************
# Trends in climate from Daymet data for the FIA dominant species sites
# 
# Author: Xiaojie Gao
# Date: 2024-07-19
# ******************************************************************************
rm(list = ls())
library(data.table)
library(lubridate)
library(trend)



# Calculate trends
SiteTrend <- function(spr_dm, tavg_name) {
    ids <- unique(spr_dm$ID)

    # Iterate all the fit, use Theil-Sens and Mann-Kendall to analyze trend
    trend_dt <- lapply(ids, function(id) {
        thefit <- spr_dm[ID == id, ]
        yrs <- thefit$year

        # ~ Spring ----------------------------------------
        tavg <- thefit[, get(tavg_name)]

        pval <- mk.test(tavg)$p.value
        slp <- sens.slope(tavg)$estimates

        arow <- data.table(
            ID = id, slp = slp, pval = pval
        )
        return(arow)
    })
    trend_dt <- do.call(rbind, trend_dt)

    return(trend_dt)
}



dm_dt <- fread("data/fia_domspec_daymet.csv")
dm_dt[, ID := paste0("X", site)]
dm_dt[, date := as_date(paste0(year, "-01-01")) + yday - 1]
dm_dt[, mth := month(date)]


# ~ Temperature trends ####
# ~ ----------------------------------------------------------------------------

# ~ Spring ----------------------------------------
spr_dm <- dm_dt[mth %in% c(3, 4, 5), .(
    spr_Tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2)
), by = c("ID", "year")]

spr_dm_trend <- SiteTrend(spr_dm, tavg_name = "spr_Tavg")


# ~ Fall ----------------------------------------
atm_dm <- dm_dt[mth %in% c(9, 10, 11), .(
    atm_Tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2)
), by = c("ID", "year")]

atm_dm_trend <- SiteTrend(atm_dm, tavg_name = "atm_Tavg")


# ~ Annual ----------------------------------------
ann_dm <- dm_dt[, .(
    ann_Tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2)
), by = c("ID", "year")]

ann_dm_trend <- SiteTrend(ann_dm, tavg_name = "ann_Tavg")


temp_trend <- rbind(
    spr_dm_trend[, season := "spr"],
    atm_dm_trend[, season := "atm"],
    ann_dm_trend[, season := "ann"]
)
# out: fia_domspec_daymet_trend.csv
fwrite(temp_trend, "pipe/fia_domspec_daymet_temp_trends.csv")



# ~ Prcp trends ####
# ~ ----------------------------------------------------------------------------

# ~ Spring ----------------------------------------
spr_dm <- dm_dt[mth %in% c(3, 4, 5), .(
    spr_Tavg = mean(`prcp (mm/day)`)
), by = c("ID", "year")]

spr_dm_trend <- SiteTrend(spr_dm, tavg_name = "spr_Tavg")


# ~ Fall ----------------------------------------
atm_dm <- dm_dt[mth %in% c(9, 10, 11), .(
    atm_Tavg = mean(`prcp (mm/day)`)
), by = c("ID", "year")]

atm_dm_trend <- SiteTrend(atm_dm, tavg_name = "atm_Tavg")


# ~ Annual ----------------------------------------
ann_dm <- dm_dt[, .(
    ann_Tavg = mean(`prcp (mm/day)`)
), by = c("ID", "year")]

ann_dm_trend <- SiteTrend(ann_dm, tavg_name = "ann_Tavg")


prcp_trend <- rbind(
    spr_dm_trend[, season := "spr"],
    atm_dm_trend[, season := "atm"],
    ann_dm_trend[, season := "ann"]
)
# out: fia_domspec_daymet_trend.csv
fwrite(prcp_trend, "pipe/fia_domspec_daymet_prcp_trends.csv")


# ~ Srad trend ####
# ~ ----------------------------------------------------------------------------

# ~ Spring ----------------------------------------
spr_dm <- dm_dt[mth %in% c(3, 4, 5), .(
    spr_Tavg = mean(`srad (W/m^2)`)
), by = c("ID", "year")]

spr_dm_trend <- SiteTrend(spr_dm, tavg_name = "spr_Tavg")


# ~ Fall ----------------------------------------
atm_dm <- dm_dt[mth %in% c(9, 10, 11), .(
    atm_Tavg = mean(`srad (W/m^2)`)
), by = c("ID", "year")]

atm_dm_trend <- SiteTrend(atm_dm, tavg_name = "atm_Tavg")


# ~ Annual ----------------------------------------
ann_dm <- dm_dt[, .(
    ann_Tavg = mean(`srad (W/m^2)`)
), by = c("ID", "year")]

ann_dm_trend <- SiteTrend(ann_dm, tavg_name = "ann_Tavg")


srad_trend <- rbind(
    spr_dm_trend[, season := "spr"],
    atm_dm_trend[, season := "atm"],
    ann_dm_trend[, season := "ann"]
)
# out: fia_domspec_daymet_trend.csv
fwrite(srad_trend, "pipe/fia_domspec_daymet_srad_trends.csv")



# ~ VP trend ####
# ~ ----------------------------------------------------------------------------

# ~ Spring ----------------------------------------
spr_dm <- dm_dt[mth %in% c(3, 4, 5), .(
    spr_Tavg = mean(`vp (Pa)`)
), by = c("ID", "year")]

spr_dm_trend <- SiteTrend(spr_dm, tavg_name = "spr_Tavg")


# ~ Fall ----------------------------------------
atm_dm <- dm_dt[mth %in% c(9, 10, 11), .(
    atm_Tavg = mean(`vp (Pa)`)
), by = c("ID", "year")]

atm_dm_trend <- SiteTrend(atm_dm, tavg_name = "atm_Tavg")


# ~ Annual ----------------------------------------
ann_dm <- dm_dt[, .(
    ann_Tavg = mean(`vp (Pa)`)
), by = c("ID", "year")]

ann_dm_trend <- SiteTrend(ann_dm, tavg_name = "ann_Tavg")


vp_trend <- rbind(
    spr_dm_trend[, season := "spr"],
    atm_dm_trend[, season := "atm"],
    ann_dm_trend[, season := "ann"]
)
# out: fia_domspec_daymet_trend.csv
fwrite(vp_trend, "pipe/fia_domspec_daymet_vp_trends.csv")

