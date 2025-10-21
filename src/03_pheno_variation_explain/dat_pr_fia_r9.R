# ******************************************************************************
# Make a table for FIA attribute data in ecoregion 9
# 
# Author: Xiaojie Gao
# Date: 2024-12-17
# ******************************************************************************
source("src/base.R")
library(terra)


# Original FIA location data in 2018
fia_true_vect <- vect("Y:/FIA/GEE/SpatialData/mostRecentFICEpointscw.shp")
# Project it to WGS84, although NAD83 and WGS84 are almost identical in Eastern
# US
fia_true_vect <- terra::project(fia_true_vect, "epsg:4326")

# According to Danelle, only use the X column.
fia_true_vect <- fia_true_vect[, c("X")]

# Table from Danelle
crosswalk <- fread("data/raw/xConcatPlotcw.csv")
# Merge crosswalk table to get the concatPlot column
fia_true_vect <- merge(fia_true_vect, crosswalk, by = "X")

# FIA plot attributes w/ fuzzed corrdinates
fia_attr <- fread(base$fia_fullwide_file)[, PLT_CN := paste0("X", PLT_CN)]

# Merge point locations with measurements
gm <- data.table(geom(fia_true_vect))[, .(x, y)]
fia_attr <- merge(
    cbind(as.data.table(fia_true_vect)[, .(X, concatPlot)], gm),
    fia_attr,
    by = "concatPlot"
)
fia_attr[, .N, keyby = "MEASYEAR"]


# Filter out sites that are less than 100% forested
fia_attr_dt <- fia_attr[nsubplotsForest == 4, ]



# out: data/fia_attr_dt_r9.csv
fwrite(fia_attr_dt, "data/fia_attr_dt_r9.csv")



