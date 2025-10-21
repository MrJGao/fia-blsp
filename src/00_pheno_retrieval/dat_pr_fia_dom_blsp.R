# ******************************************************************************
# Process BLSP for FIA sites.
# 
# Author: Xiaojie Gao
# Date: 2024-05-06
# ******************************************************************************
library(data.table)
library(blsp)
library(parallel)
library(lubridate)



# evi2_dt <- fread("data/fia_dom_landsat.csv")
evi2_dt <- fread("data/fia_dom_landsat_more.csv")
length(unique(evi2_dt$ID))


# Export a no xy version to run on the cluster
evi2_dt[, ":=" (Latitude = NULL, Longitude = NULL)]
fwrite(evi2_dt, file = "data/fia_dom_landsat_noxy.csv")



cl <- makeCluster(50)
calls <- clusterCall(cl, function() {
    suppressWarnings({
        library(data.table)
        library(blsp)
    })
})
clusterExport(cl, c("evi2_dt"))

ids <- unique(evi2_dt$ID)
blsp_fit <- clusterApplyLB(cl, x = ids, fun = function(id) {
    thedt <- evi2_dt[ID == id,]
    avgfit <- blsp::FitAvgModel(
        date_vec = thedt$Date,
        vi_vec = thedt$evi2
    )
    blsp_fit <- blsp::FitBLSP(
        date_vec = thedt$Date,
        vi_vec = thedt$evi2,
        init_values = avgfit,
        start_yr = 1984,
        end_yr = 2023,
        opt = list(method = "threshold"),
        verbose = TRUE
    )
    blsp_fit$ID <- id

    return(blsp_fit)
})
stopCluster(cl)

saveRDS(blsp_fit, file = "pipe/blsp_fit_more.Rds")

# NOTE: b/c I processed the previous BLSP data using parallel processing of the
# Harvard cluster, the previous results were saved in inidividual Rds files.
# Here, to make the rest of the analyses consistent, I export individual Rds
# files as well.
for (i in 1:length(blsp_fit)) {
    thefit <- blsp_fit[[i]]
    saveRDS(thefit, file = paste0("pipe/fia_dom_blsp/", "more_", i, ".Rds"))
}


