# ******************************************************************************
# Trend analysis on FIA BLSP time series
# 
# Author: Xiaojie Gao
# Date: 2024-10-03
# ******************************************************************************
rm(list = ls())

source("src/base.R")
library(data.table)
library(lubridate)
library(trend)
library(parallel)




# ~ Organize BLSP fit ####
# ~ ----------------------------------------------------------------------------
# Determine a species name, find all the IDs, then, check each blsp fit files to
# get the corresponding data
DoSpec <- function(specname, allfitfiles, outdir) {
    # specname <- "PITA"
    ids <- evi2_dt[Category == specname, unique(ID)]
    # B/c I updated the major species site table, here I need to filter out the
    # sites that were removed and add sites that were added!
    majorspec_dt <- fread("data/fia_domspec_truexy.csv")
    major_ids <- majorspec_dt[, paste0("X", X)]
    ids <- ids[ids %in% major_ids]
    
    message(paste("Processing", specname, "---------"))
    
    outpdf <- file.path(outdir, paste0(specname, ".pdf"))
    pdf(outpdf, width = 12, height = 4)

    # Store all found fits
    blsp_fit <- list()

    for (i in 1:length(allfitfiles)) {
        thefit <- readRDS(allfitfiles[i])
        if (!thefit$ID %in% ids) {
            next
        }

        if (ncol(thefit$phenos) < 22) {
            next
        }

        # Found 
        blsp_fit <- append(blsp_fit, list(thefit))

        yrs <- as.numeric(thefit$phenos$Year)
        
        if (length(yrs) < 39) {
            next
        }
        
        ptcols <- RColorBrewer::brewer.pal(6, "Set1")

        thedt <- thefit$data[, .(Date = date, evi2 = vi)]
        if (nrow(thedt) == 0 | all(is.na(thedt$evi2))) {
            next
        }
        
        plot(
            thedt[, .(Date, evi2)],
            cex = 0,
            xlab = "Date", ylab = "EVI2",
            bty = "L",
            mgp = c(1.5, 0.5, 0)
        )

        # Plot fitted lines
        mean_param <- lapply(thefit$params, function(ll) {
            apply(ll, 2, median)
        })
        for (j in 1:length(yrs)) {
            # The fitted line
            param_list <- lapply(mean_param, "[[", j)
            names(param_list) <- paste0("m", 1:7)
            dates <- seq(
                as_date(paste0(yrs[j], "-01-01")),
                as_date(paste0(yrs[j], "-12-31")),
                by = "day"
            )
            days <- yday(dates)
            param_list$t <- days
            pred <- eval(str2expression(thefit$model$model_str), envir = param_list)

            lines(dates, pred, col = "seagreen", lwd = 2)

            null <- lapply(seq(2, 22, by = 3), function(k) {
                pheno <- unlist(thefit$phenos[j, k, with = FALSE])
                pheno_dates <- as_date(paste0(yrs[j], "-01-01")) + pheno
                pheno_vi <- pred[pheno]
                points(pheno_dates, pheno_vi, pch = "x", col = "red")
            })
        }
        points(thedt[, .(Date, evi2)], cex = 0.5, col = "grey", pch = 16)
    }


    # SOS and EOS trends
    yrs <- 1984:2023
    par(mfrow = c(1, 2))
    all_sos <- lapply(blsp_fit, function(thefit) {
        if (nrow(thefit$phenos) < length(yrs)) {
            return(NULL)
        }
        sos <- thefit$phenos[, MidGreenup]
        sos_lwr <- thefit$phenos[, MidGreenup_lwr]
        sos_upr <- thefit$phenos[, MidGreenup_upr]
        sos[sos_upr - sos_lwr > 30] <- NA

        if (length(sos) < 39) {
            return(NULL)
        }

        return(sos)
    })
    qt <- quantile(unlist(all_sos), c(0.025, 0.975), na.rm = TRUE)
    min_sos <- qt[1]
    max_sos <- qt[2]

    if (!is.na(min_sos) & !is.na(max_sos)) {
        plot(NA,
            xlim = c(yrs[1], yrs[length(yrs)]), ylim = c(min_sos, max_sos),
            bty = "l", mgp = c(1.5, 0.5, 0),
            xlab = "Year", ylab = "MidGreenup"
        )
        null <- lapply(all_sos, function(sos) {
            if (length(sos) != length(yrs)) {
                return(NULL)
            }
            lines(yrs, sos, col = adjustcolor("grey", 0.3))
        })
        median_sos <- apply(do.call(rbind, all_sos), 2, median, na.rm = TRUE)
        lines(yrs, median_sos, lwd = 2, col = "seagreen")
    }


    all_eos <- lapply(blsp_fit, function(thefit) {
        if (nrow(thefit$phenos) < length(yrs)) {
            return(NULL)
        }
        eos <- thefit$phenos[, MidGreendown]
        eos_lwr <- thefit$phenos[, MidGreendown_lwr]
        eos_upr <- thefit$phenos[, MidGreendown_upr]
        eos[eos_upr - eos_lwr > 30] <- NA
        return(eos)
    })
    qt <- quantile(unlist(all_eos), c(0.025, 0.975), na.rm = TRUE)
    min_eos <- qt[1]
    max_eos <- qt[2]

    if (!is.na(min_eos) & !is.na(max_eos)) {
        plot(NA,
            xlim = c(yrs[1], yrs[length(yrs)]), ylim = c(min_eos, max_eos),
            bty = "l", mgp = c(1.5, 0.5, 0),
            xlab = "Year", ylab = "MidGreendown"
        )
        null <- lapply(all_eos, function(eos) {
            lines(yrs, eos, col = adjustcolor("grey", 0.3))
        })
        median_eos <- apply(do.call(rbind, all_eos), 2, median, na.rm = TRUE)
        lines(yrs, median_eos, lwd = 2, col = "darkorange")
    }

    dev.off()
}



allfitfiles <- list.files("pipe/fia_dom_blsp", ".Rds$", full.names = TRUE)

evi2_dt <- rbind(
    fread("data/fia_dom_landsat.csv"), 
    fread("data/fia_dom_landsat_more.csv")
)
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
majorspec_dt <- fread("data/fia_domspec_truexy.csv")
major_ids <- majorspec_dt[, paste0("X", X)]
evi2_dt <- evi2_dt[ID %in% major_ids]

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
clusterExport(cl, c("evi2_dt"))
null <- clusterApplyLB(
    cl, 
    specnames, DoSpec, 
    allfitfiles = allfitfiles, outdir = outdir
)
stopCluster(cl)





# ~ Site-specific trend ####
# ~ ----------------------------------------------------------------------------
SiteSpecTrend <- function(specname) {
    # Read the species file
    spec_fit <- readRDS(file.path("pipe/blsp_fit", paste0(specname, ".Rds")))

    # Use the EVI2 dt to get the latitude for each site
    evi2_dt <- rbind(
        fread("data/fia_dom_landsat.csv"), 
        fread("data/fia_dom_landsat_more.csv")
    )
    # B/c I updated the major species site table, here I need to filter out the
    # sites that were removed and add sites that were added!
    majorspec_dt <- fread("data/fia_domspec_truexy.csv")
    major_ids <- majorspec_dt[, paste0("X", X)]
    evi2_dt <- evi2_dt[ID %in% major_ids]

    yrs <- spec_fit[[1]]$phenos[, as.numeric(Year)]

    # Iterate all the fit, use a linear regression to get a slope
    trend_dt <- lapply(spec_fit, function(thefit) {
        phenos <- thefit$phenos
        if (nrow(phenos) < length(yrs) || ncol(phenos) < 22) {
            return(NULL)
        }
        setorder(phenos, Year)

        # ~ Spring ----------------------------------------
        sos <- phenos[, MidGreenup]
        sos_lwr <- phenos[, MidGreenup_lwr]
        sos_upr <- phenos[, MidGreenup_upr]
        sos[sos_upr - sos_lwr > 30] <- NA
        
        sos_noNA <- na.omit(sos)
        if (length(sos_noNA) > 10) {
            spr_pval <- mk.test(sos_noNA)$p.value
            spr_slp <- sens.slope(sos_noNA)$estimates
        } else {
            spr_pval <- NA
            spr_slp <- NA
        }
        

        # ~ Autumn ----------------------------------------
        eos <- phenos[, MidGreendown]
        eos_lwr <- phenos[, MidGreendown_lwr]
        eos_upr <- phenos[, MidGreendown_upr]
        eos[eos_upr - eos_lwr > 30] <- NA

        eos_noNA <- na.omit(eos)
        if (length(eos_noNA) > 10) {
            atm_pval <- mk.test(eos_noNA)$p.value
            atm_slp <- sens.slope(eos_noNA)$estimates
        } else {
            atm_pval <- NA
            atm_slp <- NA
        }


        # ~ Growing season length ----------------------------------------
        gsl <- eos - sos
        gsl_noNA <- na.omit(gsl)
        if (length(gsl_noNA) > 10) {
            gsl_pval <- mk.test(gsl_noNA)$p.value
            gsl_slp <- sens.slope(gsl_noNA)$estimates
        } else {
            gsl_pval <- NA
            gsl_slp <- NA
        }

        # Get latitude
        lat <- evi2_dt[ID == thefit$ID, Latitude][1]

        arow <- data.table(
            spec = specname, ID = thefit$ID, lat = lat,
            spr_slp = spr_slp, spr_pval = spr_pval,
            atm_slp = atm_slp, atm_pval = atm_pval,
            gsl_slp = gsl_slp, gsl_pval = gsl_pval
        )
        return(arow)
    })
    trend_dt <- do.call(rbind, trend_dt)
    
    return(trend_dt)
}



evi2_dt <- rbind(
    fread("data/fia_dom_landsat.csv"), 
    fread("data/fia_dom_landsat_more.csv")
)
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
majorspec_dt <- fread("data/fia_domspec_truexy.csv")
major_ids <- majorspec_dt[, paste0("X", X)]
evi2_dt <- evi2_dt[ID %in% major_ids]

specnames <- unique(evi2_dt$Category)
specnames <- specnames[specnames != "PITA"]
site_trend_dt <- NULL
for (specname in specnames) {
    site_trend_dt <- rbind(site_trend_dt, SiteSpecTrend(specname))
}
site_trend_dt <- merge(site_trend_dt, 
    unique(evi2_dt[, .(ID, lon = Longitude)]), 
    by = "ID"
)
# Get deci/ever info
fia_domspec_sites <- fread("data/fia_domspec_truexy.csv")
fia_domspec_sites <- fia_domspec_sites[, .(spec, CommonName, Type)] %>%
    unique()
site_trend_dt <- merge(site_trend_dt, fia_domspec_sites, by = "spec")

fwrite(site_trend_dt, "pipe/fia_dom_blsp_trend.csv")


