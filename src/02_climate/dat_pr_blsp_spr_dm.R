# ******************************************************************************
# Format BLSP and Daymet climate data together
# 
# Author: Xiaojie Gao
# Date: 2024-06-26
# ******************************************************************************
rm(list = ls())

library(data.table)
library(parallel)
library(lubridate)



# Get and format a site-year's Daymet temperature
MergePhenoDm <- function(i) {
    cur_sy <- blsp_dt[i, ]
    end_date <- as_date(
        255 - 1,
        origin = paste0(cur_sy[["PhenoYear"]], "-01-01")
    )
    start_date <- end_date - 364
    if (leap_year(as_date(paste0(as.numeric(cur_sy[["PhenoYear"]]) - 1, 
        "-01-01")))) {
        start_date <- start_date - 1
    }

    arow <- dm_dt[ID == cur_sy$siteID & between(date, start_date, end_date), ]
    arow$SOS <- cur_sy$SOS
    arow$PhenoYear <- cur_sy$PhenoYear
    arow$spec <- cur_sy$spec

    return(arow)
}


blsp <- fread("data/fia_blsp_fit.csv")
# Remove large uncertainty observations
blsp[MidGreenup_upr - MidGreenup_lwr > 30, MidGreenup := NA]
blsp[MidGreendown_upr - MidGreendown_lwr > 30, MidGreendown := NA]


dm_dt <- fread("data/fia_domspec_daymet.csv")
dm_dt[, ID := paste0("X", site)]
dm_dt[, date := as_date(paste0(year, "-01-01")) + yday - 1]


blsp_dt <- blsp[Year > 1984, .(siteID = ID, spec, PhenoYear = Year, SOS = MidGreenup)]


cl <- makeCluster(60, outfile = "")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        library(data.table)
        library(lubridate)
    })
})
clusterExport(cl, c("dm_dt", "blsp_dt"))
blsp_dm_dt <- clusterApplyLB(cl, x = 1:nrow(blsp_dt), fun = MergePhenoDm)
stopCluster(cl)

blsp_dm_dt <- data.table(do.call(rbind, blsp_dm_dt))


# out: data/blsp_dm_dt.csv
fwrite(blsp_dm_dt, file = "data/blsp_dm_dt.csv")


