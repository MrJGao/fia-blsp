# ******************************************************************************
# Select dominant species for each FIA site.
# 
# The sites and their dominant species are determined by:
#   1. The site must have been measured at least twice and the latest
#      measurement year is within 2013-2023, i.e., during the last 10 years.
#   2. We use baseal area 60% as a threshold to determine the dominant species
#      for a site, and the dominant species must agree across multi-year
#      measurements. Sites with none or multiple species having a basal area
#      exceeding the threshold were removed. 
# 
# Author: Xiaojie Gao
# Date: 2024-04-01
# ******************************************************************************
rm(list=ls())

source("src/base.R")
library(terra)
library(magrittr)
library(data.table)
library(readxl)
library(parallel)



# Get a major species for a site
GetSiteMajorSpec <- function(site_id) {
    # site_id = fia_attr$concatPlot[1]

    theplot <- fia_attr[concatPlot == site_id, ]
    
    spec_cols <- grep(
        paste0("spBAm2ha_*.*_live"),
        fia_attrnames,
        value = TRUE
    )

    setorder(theplot, MEASYEAR)
    # The earliest record
    theplot_early <- theplot[1, ]
    pct_early <- theplot_early[, spec_cols, with = FALSE] / 
        theplot_early$plotGroupBAm2ha_live
    dom_early <- spec_cols[which(pct_early > 0.6)]
    # The newest record
    theplot_now <- theplot[nrow(theplot), ]
    pct_now <- theplot_now[, spec_cols, with = FALSE] / 
        theplot_now$plotGroupBAm2ha_live
    dom_now <- spec_cols[which(pct_now > 0.6)]

    if (length(dom_early) == length(dom_now) && length(dom_early) == 1) {
        dom_spec <- substr(dom_early, 10, 13)
        arow <- cbind(
            site_id, dom_spec, 
            theplot_early$MEASYEAR, theplot_now$MEASYEAR
        )
        return(arow)
    } else if (length(dom_early) == length(dom_now) && length(dom_early) > 1) {
        dom_early <- spec_cols[which.max(pct_early)]
        dom_now <- spec_cols[which.max(pct_now)]
        if (dom_early == dom_now) {
            dom_spec <- substr(dom_early, 10, 13)
            arow <- cbind(
                site_id, dom_spec,
                theplot_early$MEASYEAR, theplot_now$MEASYEAR
            )
            return(arow)
        } else {
            return(NULL)
        }
    } else {
        return(NULL)
    }
}



# ~ Select major species sites ####
# ~ ----------------------------------------------------------------------------
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
fia_attr <- fia_attr[nsubplotsForest == 4, ]

fia_attrnames <- names(fia_attr)

uni_concatplot <- unique(fia_attr$concatPlot)


# Make it parallel and do the work
cl <- makeCluster(30, outfile = "")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/base.R")
        library(terra)
        library(magrittr)
        library(data.table)
        library(readxl)
    })
})
clusterExport(cl, c("fia_attr", "fia_attrnames"))
majorspec <- clusterApplyLB(cl, x = uni_concatplot, fun = GetSiteMajorSpec)
stopCluster(cl)


majorspec_dt <- do.call(rbind, majorspec)
majorspec_dt <- data.table(majorspec_dt)
colnames(majorspec_dt) <- c("concatPlot", "domspec", "earliest", "latest")

# Check time span and select only the sites within the last 10 years
majorspec_dt[, timespan := as.numeric(latest) - as.numeric(earliest)]
majorspec_dt <- majorspec_dt[latest > 2013, ]

# Merge w/ the FIA table to get coordinates
majorspec_dt <- merge(
    fia_attr[MEASYEAR %in% 2013:2023, .(concatPlot, X, x, y)], 
    majorspec_dt, 
    by = "concatPlot"
)
majorspec_dt <- unique(majorspec_dt)
colnames(majorspec_dt)[5] <- "spec"

# Remove sites that have too few records
# majorspec_dt[, .N, by = spec]
majorspec_dt[, N := .N, by = spec]
majorspec_dt <- majorspec_dt[N > 100, ]



# ~ Pre analysis ####
# ~ ----------------------------------------------------------------------------

# ~ What are the species? ----------------------------------------

domspec <- unique(majorspec_dt$spec)

# Look them up in the FIA code table
fia_spec_codes <- read_xlsx(base$fia_species_codes_xlxs, skip = 6)
colnames(fia_spec_codes) <- gsub("\r\n", "", colnames(fia_spec_codes))
setDT(fia_spec_codes)

idx <- sapply(domspec, function(x) {
    grep(x, fia_spec_codes$PLANTSCode)
})
domspec_dt <- lapply(domspec, function(x) {
    fia_spec_codes[grep(x, fia_spec_codes$PLANTSCode), .(
        PLANTSCode, CommonName, Genus, Species
    )]
}) %>%
    do.call(rbind, .)


# Match species w/ their common names
majorspec_dt[, c("CommonName", "Genus", "Species") := ""]
for (sp in domspec) {
    set(
        majorspec_dt, which(majorspec_dt$spec == sp),
        c("CommonName", "Genus", "Species"),
        fia_spec_codes[grep(sp, fia_spec_codes$PLANTSCode), .(
            CommonName, Genus, Species
        )][1, ]
    )
}


# Categorize them as deciduous and evergreen
deci_codes <- c(
    "LITU", "ACRU", "LIST", "QUAL", "CAGL", "QUVE", "FRPE", "QUPR", "FRAM",
    "ROPS", "QUMA", "TIAM", "QURU", "PRSE", "ULAM", "ACSA", "FAGR", "CAOV",
    "BEAL", "LALA", "BEPO", "POTR", "PRPE", "BEPA", "QUCO", "BELE", "ULRU",
    "FRNI", "QUST", "CAAL", "QUFA", "ULAL", "ACPE", "NYSY", "COFL", "OXAR"
)
ever_codes <- c(
    "JUVI", "JUVIS", "PITA", "PIEC", "PIST", "PIRE", "PIMA", "THOC", "PIRU",
    "ABBA", "TSCA", "PIGL", "PIAB", "PIRI", "PIBA", "PIVI"
)


majorspec_dt[spec %in% deci_codes, Type := "deci"]
majorspec_dt[spec %in% ever_codes, Type := "ever"]

majorspec_dt <- majorspec_dt[Type != "", ]

majorspec_dt[, .(.N, Type), by = Type]


# out: fia_domspec_truexy.csv
fwrite(majorspec_dt, "data/fia_domspec_truexy.csv")
# out: domspec_sites.csv
domspec_sites <- majorspec_dt[, .(
    ID = X, Category = Type, Latitude = y, Longitude = x
)]
fwrite(domspec_sites, "data/domspec_sites.csv")


# ~ Where are those plots? ----------------------------------------

majorspec_dt <- fread("data/fia_domspec_truexy.csv") 


majorspec_vect <- vect(majorspec_dt, geom = c("x", "y"), crs = "EPSG:4326")

majorspec_sf <- sf::st_as_sf(majorspec_vect)
mapview::mapview(majorspec_sf, 
    map.types = c("OpenStreetMap.Mapnik", "Esri.WorldImagery"),
    zcol = "Type",
    cex = 1
)


