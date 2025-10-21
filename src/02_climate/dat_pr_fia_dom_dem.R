# ******************************************************************************
# Extract DEM for the sites of FIA dominant species 
# 
# Author: Xiaojie Gao
# Date: 2024-04-17
# ******************************************************************************
library(data.table)
library(magrittr)


pheno_dm_dt <- fread("data/sp_60_mslsp_dm.csv")

dem_dt <- lapply(
    list.files("data/raw/fia_srtm", ".csv", full.names = TRUE), 
    fread
) %>%
    do.call(rbind, .)

pheno_dm_dt <- merge(pheno_dm_dt, 
    dem_dt[, .(ID, DEM = SRTMGL1_NC_003_SRTMGL1_DEM)], 
    by.x = "X", by.y = "ID"
)


# out: "data/sp_60_mslsp_dm.csv"
fwrite(pheno_dm_dt, file = "data/sp_60_mslsp_dm.csv")


