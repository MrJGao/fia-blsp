# ******************************************************************************
# Parse the result of BLSP fit for FIA dominant species
# 
# Author: Xiaojie Gao
# Date: 2024-05-26
# ******************************************************************************
rm(list = ls())

library(data.table)
library(lubridate)
library(parallel)



# Parse BLSP fit to species-specific obj files
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
ParseBlspFit <- function() {
    allfitfiles <- list.files("pipe/fia_dom_blsp", ".Rds$", full.names = TRUE)

    evi2_dt <- rbind(
        fread("data/fia_dom_landsat.csv"),
        fread("data/fia_dom_landsat_more.csv")
    )
    
    specnames <- unique(evi2_dt$Category)

    outdir <- "pipe/blsp_fit"
    dir.create("pipe/blsp_fit", showWarnings = FALSE)


    cl <- makeCluster(length(specnames))
    calls <- clusterCall(cl, function() {
        suppressWarnings({
            library(data.table)
            library(lubridate)
            library(parallel)
        })
    })
    # clusterExport(cl, c("evi2_dt"))
    null <- clusterApplyLB(cl, specnames, function(specname) {
        evi2_dt <- rbind(
            fread("data/fia_dom_landsat.csv"), 
            fread("data/fia_dom_landsat_more.csv")
        )
        ids <- evi2_dt[Category == specname, unique(ID)]
        majorspec_dt <- fread("data/fia_domspec_truexy.csv")
        major_ids <- majorspec_dt[, paste0("X", X)]
        ids <- ids[ids %in% major_ids]

        spec_fit <- list()
        for (i in 1:length(allfitfiles)) {
            thefit <- readRDS(allfitfiles[i])
            if (!thefit$ID %in% ids) {
                next
            }

            spec_fit <- append(spec_fit, list(thefit))
        }

        saveRDS(spec_fit, file.path(outdir, paste0(specname, ".Rds")))
    })
    stopCluster(cl)
}



ParseBlspFit()



# ~ Save the entire data to a data.table ####
# ~ ----------------------------------------------------------------------------

# Select the only the deciduous species
fia_domspec_sites <- fread("data/fia_domspec_truexy.csv")
fia_domspec_spec <- fia_domspec_sites[, .(spec, CommonName, Type)] %>%
    unique() %>%
    .[Type == "deci"]

specnames <- fia_domspec_spec$spec %>%
    sort()


# Read all data and store in a data.table
fia_domspec_sites[, ID := paste0("X", X)]

blsp_fit_dt <- NULL
for (specname in specnames) {
    spec_blsp <- readRDS(paste0(blsp_dir, specname, ".Rds"))

    spec_blsp_fit <- lapply(seq_along(spec_blsp), function(i) {
        thefit <- spec_blsp[[i]]
        
        pheno_dt <- thefit$phenos
        if (ncol(pheno_dt) < 22) {
            return(NULL)
        }
        id <- thefit$ID
        # Location
        xy <- fia_domspec_sites[ID == id, .(x, y)]

        pheno_dt[, spec := specname]
        pheno_dt[, ID := id]
        pheno_dt[, x := xy$x]
        pheno_dt[, y := xy$y]

        return(pheno_dt)
    }) %>%
        do.call(rbind, .)

    blsp_fit_dt <- rbind(blsp_fit_dt, spec_blsp_fit)
}

# Reorder columns
setcolorder(blsp_fit_dt, c("ID", "spec", "x", "y"))

# out: data/fia_blsp_fit
fwrite(blsp_fit_dt, "data/fia_blsp_fit.csv")
